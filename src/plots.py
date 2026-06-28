from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

from experiment import add_experiment_arg, apply_toml_defaults

REPO_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cvdv-matplotlib-cache"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulator import (
    DEFAULT_DT_US,
    DEFAULT_JAQAL,
    DEFAULT_TIMES_MS,
    FOCK_CUTOFF,
    X_MAX,
    X_MIN,
    X_SAMPLES,
    compiled_gate_x_trace,
    exact_hsim_x_trace,
    ho_basis_functions,
    jaqal_timestep_us,
    load_program_for_hsim_trace,
    load_program_model_and_snapshots,
    oscillator_density,
    potential_curve,
)


DEFAULT_DIRECT_OUTPUT = REPO_ROOT / "build" / "sandia_double_well_oracle_panels.png"
DEFAULT_HSIM_TRACE_OUTPUT = REPO_ROOT / "build" / "sandia_double_well_hsim_vs_jaqal_trace.png"
DEFAULT_DIRECT_TITLE = "Symmetric double well from ideal truncated-Fock simulator"
DEFAULT_HSIM_TRACE_TITLE = r"Exact $H_{\mathrm{sim}}$ versus compiled-gate dynamics"
DEFAULT_CHI_TITLE = "McGarry Fig. 7-style characteristic-function slice"
DEFAULT_COMPARISON_TITLE = "Direct state access versus McGarry characteristic-function readout"


def time_title(time_ms: float) -> str:
    return f"t = {time_ms:.2f} ms"


def plot_panel(
    axis: plt.Axes,
    xs: np.ndarray,
    potential: np.ndarray,
    density: np.ndarray,
    time_ms: float,
) -> None:
    potential_min = np.min(potential)
    potential_max = np.max(potential)
    density_max = np.max(density)
    scale = 0.18 * (potential_max - potential_min) / density_max if density_max > 0 else 0.0
    overlay = potential_min + scale * density

    axis.plot(xs, potential, color="#d97b2d", linewidth=2.0)
    axis.fill_between(xs, potential_min, overlay, color="#2f5aa6", alpha=0.35)
    axis.plot(xs, overlay, color="#2f5aa6", linewidth=1.6)
    axis.set_title(time_title(time_ms))
    axis.set_xlabel("x")
    axis.grid(alpha=0.15)


def plot_density_panels(
    densities: dict[float, np.ndarray],
    model: dict[str, float],
    times_ms: list[float],
    output: Path,
    title: str,
    ylabel: str = "Potential / shifted density",
) -> None:
    xs = np.linspace(X_MIN, X_MAX, X_SAMPLES)
    potential = potential_curve(xs, model)

    figure, axes = plt.subplots(
        1,
        len(times_ms),
        figsize=(12.0, 4.0),
        sharey=True,
        constrained_layout=True,
    )
    if len(times_ms) == 1:
        axes = [axes]

    for axis, time_ms in zip(axes, times_ms):
        plot_panel(axis, xs, potential, densities[time_ms], time_ms)

    axes[0].set_ylabel(ylabel)
    figure.suptitle(title, fontsize=13)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_direct_state_panels(
    snapshots: dict[float, np.ndarray],
    model: dict[str, float],
    cutoff: int,
    times_ms: list[float],
    output: Path,
    title: str = DEFAULT_DIRECT_TITLE,
) -> None:
    xs = np.linspace(X_MIN, X_MAX, X_SAMPLES)
    basis = ho_basis_functions(cutoff, xs)
    densities = {
        time_ms: oscillator_density(snapshots[time_ms], cutoff, basis)
        for time_ms in times_ms
    }
    plot_density_panels(
        densities,
        model,
        times_ms,
        output,
        title=title,
    )


def plot_measurement_comparison(
    results: dict[float, dict[str, object]],
    model: dict[str, float],
    times_ms: list[float],
    output: Path,
    title: str = DEFAULT_COMPARISON_TITLE,
) -> None:
    xs = np.linspace(X_MIN, X_MAX, X_SAMPLES)
    potential = potential_curve(xs, model)

    figure, axes = plt.subplots(
        2,
        len(times_ms),
        figsize=(12.0, 7.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    if len(times_ms) == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    for column, time_ms in enumerate(times_ms):
        plot_panel(axes[0, column], xs, potential, results[time_ms]["direct_density"], time_ms)
        plot_panel(axes[1, column], xs, potential, results[time_ms]["measured_density"], time_ms)
        axes[0, column].set_xlabel("")

    axes[0, 0].set_ylabel("Direct postselected density")
    axes[1, 0].set_ylabel("Eq. 33 from simulated readout")
    figure.suptitle(title, fontsize=13)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_characteristic_slice_panels(
    results: dict[float, dict[str, object]],
    times_ms: list[float],
    output: Path,
    title: str = DEFAULT_CHI_TITLE,
) -> None:
    figure, axes = plt.subplots(
        2,
        len(times_ms),
        figsize=(12.0, 6.0),
        sharex=True,
        sharey="row",
        constrained_layout=True,
    )
    if len(times_ms) == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    for column, time_ms in enumerate(times_ms):
        beta_imag = np.asarray(results[time_ms]["beta_imag"])
        chi = np.asarray(results[time_ms]["chi"])
        positive = beta_imag >= -1e-14
        beta_plot = beta_imag[positive]
        chi_plot = chi[positive]

        axes[0, column].plot(beta_plot, chi_plot.real, color="#2f5aa6", linewidth=1.8)
        axes[1, column].plot(beta_plot, chi_plot.imag, color="#b3453c", linewidth=1.8)
        axes[0, column].set_title(time_title(time_ms))
        axes[1, column].set_xlabel(r"Im[$\beta$]")
        axes[0, column].grid(alpha=0.15)
        axes[1, column].grid(alpha=0.15)
        axes[0, column].axhline(0.0, color="black", linewidth=0.8, alpha=0.25)
        axes[1, column].axhline(0.0, color="black", linewidth=0.8, alpha=0.25)

    axes[0, 0].set_ylabel(r"Re[$\chi(i y)$]")
    axes[1, 0].set_ylabel(r"Im[$\chi(i y)$]")
    figure.suptitle(title, fontsize=13)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_hsim_vs_compiled_gate_x_trace(
    model: dict[str, float],
    prep,
    max_time_ms: float,
    points: int,
    cutoff: int,
    output: Path,
    title: str = DEFAULT_HSIM_TRACE_TITLE,
) -> None:
    if points < 2:
        raise ValueError("--trace-points must be at least 2")

    times_ms = np.linspace(0.0, max_time_ms, points)
    x_expectation = exact_hsim_x_trace(prep, model, times_ms, cutoff)
    gate_times_ms, gate_x_expectation = compiled_gate_x_trace(
        prep,
        model,
        max_time_ms,
        cutoff,
    )
    xsdf_gate_times_ms, xsdf_gate_x_expectation = compiled_gate_x_trace(
        prep,
        model,
        max_time_ms,
        cutoff,
        displacement_gate="xSDF",
    )

    figure, axis = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    axis.plot(
        times_ms,
        x_expectation,
        color="#2f5aa6",
        linewidth=2.0,
        label=r"exact $H_{\mathrm{sim}}$",
    )
    axis.plot(
        gate_times_ms,
        gate_x_expectation,
        color="#b3453c",
        linewidth=1.5,
        marker=".",
        markersize=3.0,
        alpha=0.9,
        label="compiled xCD/Rz sequence",
    )
    axis.plot(
        xsdf_gate_times_ms,
        xsdf_gate_x_expectation,
        color="#3f7f4a",
        linewidth=1.2,
        linestyle="--",
        marker="x",
        markersize=3.0,
        alpha=0.85,
        label="xSDF/Rz semantic sequence",
    )
    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.3)
    axis.set_xlabel("time (ms)")
    axis.set_ylabel(r"$\langle x \rangle$")
    axis.set_title(title)
    axis.legend(frameon=False)
    axis.grid(alpha=0.15)

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


PLOTS_TOML_MAP = {
    "output.jaqal": "jaqal",
    "plots.direct_png": "output",
    "plots.hsim_png": "hsim_output",
    "plots.direct_title": "title",
    "plots.hsim_title": "hsim_title",
    "plots.times_ms": "times_ms",
    "plots.hsim_max_time_ms": "hsim_max_time_ms",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the direct ideal-simulator state density for the McGarry double well."
    )
    add_experiment_arg(parser)
    parser.add_argument("--jaqal", type=Path, default=DEFAULT_JAQAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_DIRECT_OUTPUT)
    parser.add_argument("--hsim-output", type=Path, default=DEFAULT_HSIM_TRACE_OUTPUT)
    parser.add_argument("--title", default=DEFAULT_DIRECT_TITLE)
    parser.add_argument("--hsim-title", default=DEFAULT_HSIM_TRACE_TITLE)
    parser.add_argument(
        "--no-hsim-output",
        action="store_true",
        help="Do not emit the McGarry Fig. 6-style x-expectation comparison.",
    )
    parser.add_argument("--cutoff", type=int, default=FOCK_CUTOFF)
    parser.add_argument("--times-ms", type=float, nargs="+", default=DEFAULT_TIMES_MS)
    parser.add_argument(
        "--hsim-max-time-ms",
        type=float,
        default=16.0,
        help="Maximum time for the Fig. 6-style x-expectation comparison.",
    )
    parser.add_argument("--trace-points", type=int, default=201)
    parser.add_argument(
        "--dt-us",
        type=float,
        default=None,
        help=(
            "Simulation timestep in microseconds. Defaults to the generated "
            f"Jaqal timestep comment, or {DEFAULT_DT_US:g} if absent."
        ),
    )
    parser.add_argument(
        "--interaction-frame",
        action="store_true",
        help="plot raw interaction-frame gate state instead of applying the inferred harmonic frame",
    )
    apply_toml_defaults(parser, PLOTS_TOML_MAP, path_dests={"jaqal", "output", "hsim_output"})
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dt_us = jaqal_timestep_us(args.jaqal, args.dt_us)
    model, snapshots = load_program_model_and_snapshots(
        args.jaqal,
        args.times_ms,
        args.cutoff,
        dt_us,
        use_schrodinger_frame=not args.interaction_frame,
    )
    plot_direct_state_panels(
        snapshots,
        model,
        args.cutoff,
        args.times_ms,
        args.output,
        args.title,
    )
    print(args.output)
    if not args.no_hsim_output:
        hsim_model, prep = load_program_for_hsim_trace(args.jaqal, dt_us)
        plot_hsim_vs_compiled_gate_x_trace(
            hsim_model,
            prep,
            args.hsim_max_time_ms,
            args.trace_points,
            args.cutoff,
            args.hsim_output,
            args.hsim_title,
        )
        print(args.hsim_output)


if __name__ == "__main__":
    main()
