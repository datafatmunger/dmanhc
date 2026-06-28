from __future__ import annotations

"""Fit per-timestep readout slopes and plot the inferred x expectation.

The saved notebook arrays have shape (beta_index, subcircuit_step).  This
script fits the saved measurement values against the notebook imBeta override
axis.  The expected overlays simulate the same fixed-reBeta readout sweep,
fit that simulated readout in the same way, and compare the resulting
slope-derived trace to the data.
"""

import argparse
import csv
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cvdv-matplotlib-cache"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulator import (
    FOCK_CUTOFF,
    apply_cosine_step_from_parameters,
    apply_harmonic_frame,
    apply_prep,
    compiled_gate_x_trace,
    exact_hsim_x_trace,
    hsim_hamiltonian,
    initial_state,
    jaqal_timestep_us,
    load_program_for_hsim_trace,
    postselected_motional_vector,
)
from measure import characteristic_function, motional_density_matrix


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = REPO_ROOT / "results" / "20260625_Data"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "build" / "readout_x_trace"
DEFAULT_EXPECTED_JAQAL = REPO_ROOT / "build" / "dmanh.jaqal"
DEFAULT_IM_BETAS = np.linspace(-0.4, 0.4, 5)
DEFAULT_DT_MS = 0.8 / 5.09628e3 * 1.0e3
DEFAULT_SELECTED_STEPS = [0, 13, 26, 49]
DEFAULT_EXPECTED_READOUT_RE_BETA = 0.0
DEFAULT_EXPECTED_BETA_MAPPING = "direct"
BETA_MAPPING_CHOICES = ("sandia-local", "direct", "direct-double")


def notebook_beta_to_mcgarry_beta(raw_beta: complex, beta_mapping: str) -> complex:
    if beta_mapping == "sandia-local":
        # Legacy pre-2026-06-28 compiler convention: xCD(s) was interpreted as
        # D(-i s), so s = i beta / 2.
        return -2.0j * raw_beta
    if beta_mapping == "direct":
        return raw_beta
    if beta_mapping == "direct-double":
        return 2.0 * raw_beta
    raise ValueError(f"unsupported beta mapping {beta_mapping!r}")


@dataclass(frozen=True)
class DatasetFit:
    path: Path
    label: str
    repeats: int | None
    measurements: np.ndarray
    times_ms: np.ndarray
    slopes: np.ndarray
    intercepts: np.ndarray
    slope_stderr: np.ndarray
    raw_slope_over_sqrt2: np.ndarray
    x_estimates: np.ndarray
    x_stderr: np.ndarray
    r_squared: np.ndarray
    fitted: np.ndarray
    measurement_im_chi_sign: float


@dataclass(frozen=True)
class ExpectedReadoutTrace:
    beta_mapping: str
    readout_re_beta: float
    mapped_betas: np.ndarray
    im_chi_values: np.ndarray
    slopes: np.ndarray
    intercepts: np.ndarray
    corrected_x: np.ndarray
    raw_slope_over_sqrt2: np.ndarray


def repeat_count(path: Path) -> int | None:
    match = re.search(r"_(\d+)\s+repeats$", path.name)
    if match:
        return int(match.group(1))
    return None


def discover_repeat_dirs(data_root: Path) -> list[Path]:
    dirs = sorted(path.parent for path in data_root.glob("*/*expZ_imMeas.npy"))
    if not dirs:
        raise FileNotFoundError(f"no expZ_imMeas.npy files found below {data_root}")
    return dirs


def fit_dataset(
    repeat_dir: Path,
    im_betas: np.ndarray,
    dt_ms: float,
    measurement_im_chi_sign: float,
) -> DatasetFit:
    if math.isclose(measurement_im_chi_sign, 0.0, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("--measurement-im-chi-sign must be nonzero")

    measurements = np.load(repeat_dir / "expZ_imMeas.npy")
    if measurements.ndim != 2:
        raise ValueError(f"{repeat_dir / 'expZ_imMeas.npy'} is not a 2-D array: {measurements.shape}")
    if measurements.shape[0] != im_betas.size:
        raise ValueError(
            f"{repeat_dir / 'expZ_imMeas.npy'} has {measurements.shape[0]} beta rows; "
            f"expected {im_betas.size}"
        )

    design = np.column_stack([im_betas, np.ones_like(im_betas)])
    coeffs, *_ = np.linalg.lstsq(design, measurements, rcond=None)
    slopes = coeffs[0]
    intercepts = coeffs[1]
    fitted = design @ coeffs

    residuals = measurements - fitted
    dof = max(im_betas.size - 2, 1)
    residual_variance = np.sum(residuals**2, axis=0) / dof
    sxx = float(np.sum((im_betas - im_betas.mean()) ** 2))
    slope_stderr = np.sqrt(residual_variance / sxx)

    centered = measurements - measurements.mean(axis=0, keepdims=True)
    total_variance = np.sum(centered**2, axis=0)
    r_squared = np.full(measurements.shape[1], np.nan, dtype=np.float64)
    nonzero = total_variance > np.finfo(float).eps
    r_squared[nonzero] = 1.0 - np.sum(residuals[:, nonzero] ** 2, axis=0) / total_variance[nonzero]

    times_ms = np.arange(measurements.shape[1], dtype=np.float64) * dt_ms
    return DatasetFit(
        path=repeat_dir,
        label=repeat_dir.name,
        repeats=repeat_count(repeat_dir),
        measurements=measurements,
        times_ms=times_ms,
        slopes=slopes,
        intercepts=intercepts,
        slope_stderr=slope_stderr,
        raw_slope_over_sqrt2=slopes / math.sqrt(2.0),
        x_estimates=slopes / (measurement_im_chi_sign * math.sqrt(2.0)),
        x_stderr=np.abs(slope_stderr / (measurement_im_chi_sign * math.sqrt(2.0))),
        r_squared=r_squared,
        fitted=fitted,
        measurement_im_chi_sign=measurement_im_chi_sign,
    )


def label_for(fit: DatasetFit) -> str:
    if fit.repeats is None:
        return fit.label
    return f"{fit.repeats} repeats"


def plot_x_trace(fits: list[DatasetFit], output: Path) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
    colors = ["#1f5aa6", "#b3433f", "#2d7f5e", "#6f4aa8"]

    for index, fit in enumerate(fits):
        axis.errorbar(
            fit.times_ms,
            fit.x_estimates,
            yerr=fit.x_stderr,
            color=colors[index % len(colors)],
            marker="o",
            markersize=3.0,
            linewidth=1.4,
            elinewidth=0.8,
            capsize=1.8,
            label=label_for(fit),
        )

    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.3)
    axis.set_xlabel("time (ms)")
    axis.set_ylabel(r"readout-corrected $\langle x\rangle$")
    axis.set_title(r"Per-timestep readout-slope estimate of $\langle x\rangle(t)$")
    axis.grid(alpha=0.16)
    axis.legend(frameon=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_fit_examples(
    fits: list[DatasetFit],
    im_betas: np.ndarray,
    selected_steps: list[int],
    output: Path,
) -> None:
    figure, axes = plt.subplots(
        len(fits),
        len(selected_steps),
        figsize=(3.15 * len(selected_steps), 2.45 * len(fits)),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes.reshape((len(fits), len(selected_steps)))

    beta_dense = np.linspace(float(np.min(im_betas)), float(np.max(im_betas)), 101)
    for row, fit in enumerate(fits):
        for column, step in enumerate(selected_steps):
            if step < 0 or step >= fit.measurements.shape[1]:
                raise ValueError(f"selected step {step} is outside 0..{fit.measurements.shape[1] - 1}")
            axis = axes[row, column]
            line = fit.slopes[step] * beta_dense + fit.intercepts[step]
            axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.2)
            axis.scatter(im_betas, fit.measurements[:, step], color="#1f5aa6", s=22, zorder=3)
            axis.plot(beta_dense, line, color="#b3433f", linewidth=1.5)
            axis.grid(alpha=0.15)
            if row == 0:
                axis.set_title(f"step {step}\nt={fit.times_ms[step]:.3f} ms", fontsize=9)
            if column == 0:
                axis.set_ylabel(f"{label_for(fit)}\nsaved Z")
            if row == len(fits) - 1:
                axis.set_xlabel(r"Im[$\beta$]")
            axis.text(
                0.03,
                0.94,
                rf"$\langle x\rangle={fit.x_estimates[step]:+.3f}$",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
            )

    figure.suptitle("Representative per-timestep line fits", fontsize=11)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def select_fit_for_repeats(fits: list[DatasetFit], repeats: int) -> DatasetFit:
    for fit in fits:
        if fit.repeats == repeats:
            return fit
    labels = ", ".join(label_for(fit) for fit in fits)
    raise ValueError(f"could not find {repeats}-repeat dataset among: {labels}")


def comparison_fits(fits: list[DatasetFit], requested_repeats: list[int] | None) -> list[DatasetFit]:
    if requested_repeats is None or not requested_repeats:
        return fits
    return [select_fit_for_repeats(fits, repeats) for repeats in requested_repeats]


def comparison_slug(fit: DatasetFit) -> str:
    if fit.repeats is not None:
        return f"{fit.repeats}_repeats"
    slug = re.sub(r"[^A-Za-z0-9]+", "_", fit.label).strip("_").lower()
    return slug or "dataset"


def mapped_readout_betas(
    readout_re_beta: float,
    im_betas: np.ndarray,
    beta_mapping: str,
) -> np.ndarray:
    return np.array(
        [
            notebook_beta_to_mcgarry_beta(complex(readout_re_beta, im_beta), beta_mapping)
            for im_beta in im_betas
        ],
        dtype=np.complex128,
    )


def fit_expected_readout_trace(
    rhos: list[np.ndarray],
    im_betas: np.ndarray,
    mapped_betas: np.ndarray,
    cutoff: int,
    measurement_im_chi_sign: float,
    beta_mapping: str,
    readout_re_beta: float,
) -> ExpectedReadoutTrace:
    values = np.empty((im_betas.size, len(rhos)), dtype=np.float64)
    for row, beta in enumerate(mapped_betas):
        for column, rho in enumerate(rhos):
            values[row, column] = characteristic_function(beta, rho, cutoff).imag

    design = np.column_stack([im_betas, np.ones_like(im_betas)])
    coeffs, *_ = np.linalg.lstsq(design, values, rcond=None)
    slopes = coeffs[0]
    intercepts = coeffs[1]
    corrected_x = slopes / math.sqrt(2.0)
    return ExpectedReadoutTrace(
        beta_mapping=beta_mapping,
        readout_re_beta=readout_re_beta,
        mapped_betas=mapped_betas,
        im_chi_values=values,
        slopes=slopes,
        intercepts=intercepts,
        corrected_x=corrected_x,
        raw_slope_over_sqrt2=measurement_im_chi_sign * corrected_x,
    )


def compiled_gate_rhos(
    prep,
    model: dict[str, float],
    times_ms: np.ndarray,
    cutoff: int,
    readout_state: str = "down",
) -> list[np.ndarray]:
    max_time_s = float(times_ms[-1]) * 1e-3
    dt = model["dt"]
    steps = round(max_time_s / dt)
    if not math.isclose(steps * dt, max_time_s, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("expected times must end on an integer multiple of the inferred timestep")

    expected_times = np.arange(steps + 1, dtype=np.float64) * dt * 1e3
    if expected_times.shape != times_ms.shape or not np.allclose(expected_times, times_ms, atol=1e-9):
        raise ValueError("readout-simulated expected traces require the saved data timestep grid")

    state = apply_prep(initial_state(cutoff), prep, cutoff)
    rhos: list[np.ndarray] = [motional_density_matrix(state, cutoff, readout_state)[0]]

    for step_index in range(steps):
        phase = model["alpha_phase_offset"] + model["delta"] * step_index * dt
        alpha = model["alpha0"] * complex(math.cos(phase), math.sin(phase))
        state = apply_cosine_step_from_parameters(
            state,
            alpha,
            model["dt"] * model["amplitude"],
            model["varphi"],
            cutoff,
        )
        elapsed = (step_index + 1) * dt
        snapshot = apply_harmonic_frame(state, model["delta"], elapsed, cutoff)
        rhos.append(motional_density_matrix(snapshot, cutoff, readout_state)[0])

    return rhos


def exact_hsim_rhos(
    prep,
    model: dict[str, float],
    times_ms: np.ndarray,
    cutoff: int,
) -> list[np.ndarray]:
    prepared_state = apply_prep(initial_state(cutoff), prep, cutoff)
    motional_state = postselected_motional_vector(prepared_state, cutoff, "down")

    hamiltonian = hsim_hamiltonian(model, cutoff)
    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    coefficients = eigenvectors.conj().T @ motional_state

    rhos: list[np.ndarray] = []
    for time_s in times_ms * 1e-3:
        evolved = eigenvectors @ (np.exp(-1.0j * energies * time_s) * coefficients)
        rhos.append(np.outer(evolved, evolved.conj()))
    return rhos


def expected_traces(
    jaqal: Path,
    times_ms: np.ndarray,
    cutoff: int,
) -> tuple[np.ndarray, np.ndarray]:
    dt_us = jaqal_timestep_us(jaqal)
    model, prep = load_program_for_hsim_trace(jaqal, dt_us)

    compiled_times, compiled_x = compiled_gate_x_trace(
        prep,
        model,
        float(times_ms[-1]),
        cutoff,
    )
    if compiled_times.shape != times_ms.shape or not np.allclose(compiled_times, times_ms, atol=1e-9):
        compiled_x = np.interp(times_ms, compiled_times, compiled_x)

    exact_x = exact_hsim_x_trace(prep, model, times_ms, cutoff)
    return compiled_x, exact_x


def expected_readout_traces(
    jaqal: Path,
    times_ms: np.ndarray,
    im_betas: np.ndarray,
    readout_re_beta: float,
    beta_mapping: str,
    cutoff: int,
    measurement_im_chi_sign: float,
) -> tuple[ExpectedReadoutTrace, ExpectedReadoutTrace, np.ndarray, np.ndarray]:
    dt_us = jaqal_timestep_us(jaqal)
    model, prep = load_program_for_hsim_trace(jaqal, dt_us)
    mapped_betas = mapped_readout_betas(readout_re_beta, im_betas, beta_mapping)

    compiled_rhos = compiled_gate_rhos(prep, model, times_ms, cutoff)
    exact_rhos = exact_hsim_rhos(prep, model, times_ms, cutoff)
    compiled_readout = fit_expected_readout_trace(
        compiled_rhos,
        im_betas,
        mapped_betas,
        cutoff,
        measurement_im_chi_sign,
        beta_mapping,
        readout_re_beta,
    )
    exact_readout = fit_expected_readout_trace(
        exact_rhos,
        im_betas,
        mapped_betas,
        cutoff,
        measurement_im_chi_sign,
        beta_mapping,
        readout_re_beta,
    )

    compiled_times, compiled_direct_x = compiled_gate_x_trace(
        prep,
        model,
        float(times_ms[-1]),
        cutoff,
    )
    if compiled_times.shape != times_ms.shape or not np.allclose(compiled_times, times_ms, atol=1e-9):
        compiled_direct_x = np.interp(times_ms, compiled_times, compiled_direct_x)
    exact_direct_x = exact_hsim_x_trace(prep, model, times_ms, cutoff)

    return compiled_readout, exact_readout, compiled_direct_x, exact_direct_x


def plot_expected_overlay(
    fit: DatasetFit,
    compiled_readout: ExpectedReadoutTrace,
    exact_readout: ExpectedReadoutTrace,
    output: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(7.3, 4.4), constrained_layout=True)

    axis.plot(
        fit.times_ms,
        exact_readout.corrected_x,
        color="#1f5aa6",
        linewidth=2.0,
        label=r"exact $H_{\mathrm{sim}}$ readout fit",
    )
    axis.plot(
        fit.times_ms,
        compiled_readout.corrected_x,
        color="#b3433f",
        linewidth=1.6,
        linestyle="--",
        label="compiled gate readout fit",
    )
    axis.errorbar(
        fit.times_ms,
        fit.x_estimates,
        yerr=fit.x_stderr,
        color="#2d7f5e",
        marker="o",
        markersize=3.0,
        linewidth=1.2,
        elinewidth=0.8,
        capsize=1.8,
        label=f"{label_for(fit)} slope extraction",
    )
    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.28)
    axis.set_xlabel("time (ms)")
    axis.set_ylabel(r"readout-slope $\langle x\rangle$")
    axis.set_title("Sign-fixed readout slope trace", fontsize=10)
    axis.grid(alpha=0.16)
    axis.legend(loc="best", fontsize=8, frameon=True)
    figure.suptitle(
        (
            r"Expected fixed-reBeta readout-slope trace versus saved data"
            "\n"
            f"raw reBeta={compiled_readout.readout_re_beta:+.3f}, "
            f"beta mapping={compiled_readout.beta_mapping}; mapped beta values are in the CSV"
        ),
        fontsize=11,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def write_expected_csv(
    fit: DatasetFit,
    compiled_readout: ExpectedReadoutTrace,
    exact_readout: ExpectedReadoutTrace,
    compiled_direct_x: np.ndarray,
    exact_direct_x: np.ndarray,
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        fieldnames = [
            "step",
            "time_ms",
            "data_x",
            "data_x_stderr",
            "raw_slope_over_sqrt2",
            "compiled_gate_readout_x",
            "exact_hsim_readout_x",
            "compiled_gate_raw_slope_over_sqrt2",
            "exact_hsim_raw_slope_over_sqrt2",
            "compiled_gate_direct_x",
            "exact_hsim_direct_x",
            "data_minus_compiled_readout",
            "raw_slope_over_sqrt2_minus_compiled_readout",
            "data_minus_exact_readout",
            "raw_slope_over_sqrt2_minus_exact_readout",
            "measurement_im_chi_sign",
            "readout_re_beta",
            "beta_mapping",
            "mapped_betas",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        mapped = " ".join(
            f"{beta.real:+.17g}{beta.imag:+.17g}i" for beta in compiled_readout.mapped_betas
        )
        for step in range(fit.times_ms.size):
            data_x = fit.x_estimates[step]
            writer.writerow(
                {
                    "step": step,
                    "time_ms": f"{fit.times_ms[step]:.12g}",
                    "data_x": f"{data_x:.17g}",
                    "data_x_stderr": f"{fit.x_stderr[step]:.17g}",
                    "raw_slope_over_sqrt2": f"{fit.raw_slope_over_sqrt2[step]:.17g}",
                    "compiled_gate_readout_x": f"{compiled_readout.corrected_x[step]:.17g}",
                    "exact_hsim_readout_x": f"{exact_readout.corrected_x[step]:.17g}",
                    "compiled_gate_raw_slope_over_sqrt2": f"{compiled_readout.raw_slope_over_sqrt2[step]:.17g}",
                    "exact_hsim_raw_slope_over_sqrt2": f"{exact_readout.raw_slope_over_sqrt2[step]:.17g}",
                    "compiled_gate_direct_x": f"{compiled_direct_x[step]:.17g}",
                    "exact_hsim_direct_x": f"{exact_direct_x[step]:.17g}",
                    "data_minus_compiled_readout": f"{data_x - compiled_readout.corrected_x[step]:.17g}",
                    "raw_slope_over_sqrt2_minus_compiled_readout": f"{fit.raw_slope_over_sqrt2[step] - compiled_readout.raw_slope_over_sqrt2[step]:.17g}",
                    "data_minus_exact_readout": f"{data_x - exact_readout.corrected_x[step]:.17g}",
                    "raw_slope_over_sqrt2_minus_exact_readout": f"{fit.raw_slope_over_sqrt2[step] - exact_readout.raw_slope_over_sqrt2[step]:.17g}",
                    "measurement_im_chi_sign": f"{fit.measurement_im_chi_sign:.17g}",
                    "readout_re_beta": f"{compiled_readout.readout_re_beta:.17g}",
                    "beta_mapping": compiled_readout.beta_mapping,
                    "mapped_betas": mapped,
                }
            )


def write_csv(fits: list[DatasetFit], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        fieldnames = [
            "dataset",
            "repeats",
            "step",
            "time_ms",
            "slope",
            "slope_stderr",
            "raw_slope_over_sqrt2",
            "x_estimate",
            "x_stderr",
            "intercept",
            "r_squared",
            "measurement_im_chi_sign",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for fit in fits:
            for step in range(fit.measurements.shape[1]):
                writer.writerow(
                    {
                        "dataset": fit.label,
                        "repeats": "" if fit.repeats is None else fit.repeats,
                        "step": step,
                        "time_ms": f"{fit.times_ms[step]:.12g}",
                        "slope": f"{fit.slopes[step]:.17g}",
                        "slope_stderr": f"{fit.slope_stderr[step]:.17g}",
                        "raw_slope_over_sqrt2": f"{fit.raw_slope_over_sqrt2[step]:.17g}",
                        "x_estimate": f"{fit.x_estimates[step]:.17g}",
                        "x_stderr": f"{fit.x_stderr[step]:.17g}",
                        "intercept": f"{fit.intercepts[step]:.17g}",
                        "r_squared": f"{fit.r_squared[step]:.17g}",
                        "measurement_im_chi_sign": f"{fit.measurement_im_chi_sign:.17g}",
                    }
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit saved readout versus Im[beta] at each timestep and plot <x>(t)."
    )
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--repeat-dir", type=Path, action="append", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--dt-ms", type=float, default=DEFAULT_DT_MS)
    parser.add_argument("--im-betas", type=float, nargs="+", default=DEFAULT_IM_BETAS.tolist())
    parser.add_argument("--selected-steps", type=int, nargs="+", default=DEFAULT_SELECTED_STEPS)
    parser.add_argument("--expected-jaqal", type=Path, default=DEFAULT_EXPECTED_JAQAL)
    parser.add_argument("--expected-readout-re-beta", type=float, default=DEFAULT_EXPECTED_READOUT_RE_BETA)
    parser.add_argument(
        "--expected-beta-mapping",
        choices=BETA_MAPPING_CHOICES,
        default=DEFAULT_EXPECTED_BETA_MAPPING,
        help=(
            "How to map the raw reBeta+i*imBeta notebook readout coordinate "
            "to McGarry beta before simulating the expected readout trace."
        ),
    )
    parser.add_argument(
        "--comparison-repeats",
        type=int,
        nargs="*",
        default=None,
        help=(
            "Repeat counts to overlay against expected dynamics. Defaults to "
            "all discovered datasets."
        ),
    )
    parser.add_argument("--cutoff", type=int, default=FOCK_CUTOFF)
    parser.add_argument(
        "--measurement-im-chi-sign",
        type=float,
        default=-1.0,
        help=(
            "Sign relating saved measurement values to Im[chi]. Use -1 for "
            "the current saved notebook arrays, where the probe readout gives -Im[chi]."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    im_betas = np.asarray(args.im_betas, dtype=np.float64)
    if im_betas.ndim != 1 or im_betas.size < 2:
        raise ValueError("--im-betas must provide at least two values")

    repeat_dirs = args.repeat_dir if args.repeat_dir else discover_repeat_dirs(args.data_root)
    fits = [
        fit_dataset(repeat_dir, im_betas, args.dt_ms, args.measurement_im_chi_sign)
        for repeat_dir in repeat_dirs
    ]

    x_trace_output = args.output_dir / "x_trace_from_slopes.png"
    fit_examples_output = args.output_dir / "selected_step_slope_fits.png"
    csv_output = args.output_dir / "x_trace_from_slopes.csv"

    plot_x_trace(fits, x_trace_output)
    plot_fit_examples(fits, im_betas, args.selected_steps, fit_examples_output)
    write_csv(fits, csv_output)

    print(x_trace_output)
    print(fit_examples_output)
    print(csv_output)

    for fit in comparison_fits(fits, args.comparison_repeats):
        slug = comparison_slug(fit)
        overlay_output = args.output_dir / f"actual_vs_expected_{slug}.png"
        overlay_csv_output = args.output_dir / f"actual_vs_expected_{slug}.csv"
        compiled_readout, exact_readout, compiled_direct_x, exact_direct_x = expected_readout_traces(
            args.expected_jaqal,
            fit.times_ms,
            im_betas,
            args.expected_readout_re_beta,
            args.expected_beta_mapping,
            args.cutoff,
            fit.measurement_im_chi_sign,
        )
        plot_expected_overlay(fit, compiled_readout, exact_readout, overlay_output)
        write_expected_csv(
            fit,
            compiled_readout,
            exact_readout,
            compiled_direct_x,
            exact_direct_x,
            overlay_csv_output,
        )
        print(overlay_output)
        print(overlay_csv_output)


if __name__ == "__main__":
    main()
