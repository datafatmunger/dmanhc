from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from dataclasses import dataclass
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
    FOCK_CUTOFF,
    apply_prep,
    displacement_matrix,
    hsim_hamiltonian,
    initial_state,
    jaqal_timestep_us,
    load_program_for_hsim_trace,
    position_matrix,
    postselected_motional_vector,
)


DEFAULT_JAQAL = REPO_ROOT / "build" / "dmanh.jaqal"
DEFAULT_CSV_OUTPUT = REPO_ROOT / "build" / "dmanh_frequency_sweep.csv"
DEFAULT_SPECTRA_OUTPUT = REPO_ROOT / "build" / "dmanh_frequency_spectra.png"
DEFAULT_TRACE_OUTPUT = REPO_ROOT / "build" / "dmanh_frequency_trace_examples.png"

Z_TO_X = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128) / math.sqrt(2.0)


@dataclass(frozen=True)
class Spectrum:
    frequencies_hz: np.ndarray
    amplitudes: np.ndarray
    peak_hz: float
    peak_amplitude: float


@dataclass(frozen=True)
class SweepRow:
    vartheta: float
    dt_us: float
    steps: int
    actual_time_ms: float
    xcd_gates: int
    rz_gates: int
    total_evolution_gates: int
    exact_peak_hz: float
    compiled_peak_hz: float
    compiled_global_peak_hz: float
    frequency_shift_hz: float
    frequency_shift_rad_s: float
    relative_frequency_shift: float
    rms_x_error: float
    max_abs_x_error: float


class ExactEvolution:
    def __init__(self, prep, model: dict[str, float], cutoff: int) -> None:
        prepared_state = apply_prep(initial_state(cutoff), prep, cutoff)
        motional_state = postselected_motional_vector(prepared_state, cutoff, "down")
        hamiltonian = hsim_hamiltonian(model, cutoff)
        energies, eigenvectors = np.linalg.eigh(hamiltonian)

        self.energies = energies
        self.coefficients = eigenvectors.conj().T @ motional_state
        self.x_in_energy_basis = eigenvectors.conj().T @ position_matrix(cutoff) @ eigenvectors

    def x_trace(self, times_s: np.ndarray) -> np.ndarray:
        evolved = np.exp(-1.0j * np.outer(times_s, self.energies)) * self.coefficients
        values = np.einsum(
            "ti,ij,tj->t",
            evolved.conj(),
            self.x_in_energy_basis,
            evolved,
            optimize=True,
        )
        return values.real.astype(np.float64)


class DisplacementApplicator:
    """Apply the same truncated-Fock displacement matrices used by simulator.py.

    The sweep needs thousands of displacements with a fixed magnitude and
    changing phase.  In the truncated Fock basis,
    D(r exp(i phi)) = exp(i phi N) D(r) exp(-i phi N), so one matrix
    exponential per magnitude is enough.
    """

    def __init__(self, cutoff: int) -> None:
        self.cutoff = cutoff
        self.levels = np.arange(cutoff, dtype=np.float64)
        self.base_by_magnitude: dict[float, np.ndarray] = {}

    def base_matrix(self, magnitude: float) -> np.ndarray:
        if math.isclose(magnitude, 0.0, rel_tol=0.0, abs_tol=1e-15):
            return np.eye(self.cutoff, dtype=np.complex128)
        key = float(magnitude)
        matrix = self.base_by_magnitude.get(key)
        if matrix is None:
            matrix = displacement_matrix(complex(magnitude, 0.0), self.cutoff)
            self.base_by_magnitude[key] = matrix
        return matrix

    def apply(self, vector: np.ndarray, alpha: complex) -> np.ndarray:
        magnitude = abs(alpha)
        if math.isclose(magnitude, 0.0, rel_tol=0.0, abs_tol=1e-15):
            return vector.copy()

        phase = math.atan2(alpha.imag, alpha.real)
        rotation = np.exp(1.0j * self.levels * phase)
        base = self.base_matrix(magnitude)
        return rotation * (base @ (rotation.conj() * vector))


def sqr_matrix(vartheta: float, varphi: float) -> np.ndarray:
    axis = np.array(
        [
            [math.cos(varphi), -1.0j * math.sin(varphi)],
            [1.0j * math.sin(varphi), -math.cos(varphi)],
        ],
        dtype=np.complex128,
    )
    return math.cos(0.5 * vartheta) * np.eye(2, dtype=np.complex128) - 1.0j * math.sin(
        0.5 * vartheta
    ) * axis


def apply_xcd_fast(
    amplitudes: np.ndarray,
    alpha: complex,
    displacement: DisplacementApplicator,
) -> np.ndarray:
    x_amplitudes = Z_TO_X.conj().T @ amplitudes
    updated = x_amplitudes.copy()
    updated[0] = displacement.apply(x_amplitudes[0], alpha)
    updated[1] = displacement.apply(x_amplitudes[1], -alpha)
    return Z_TO_X @ updated


def schrodinger_down_x_expectation(
    amplitudes: np.ndarray,
    model: dict[str, float],
    elapsed_s: float,
    x_operator: np.ndarray,
) -> float:
    cutoff = x_operator.shape[0]
    levels = np.arange(cutoff, dtype=np.float64)
    oscillator_phases = np.exp(-1.0j * model["delta"] * (levels + 0.5) * elapsed_s)
    vector = amplitudes[1] * oscillator_phases
    norm_squared = float(np.vdot(vector, vector).real)
    if norm_squared <= 0.0:
        raise ValueError("postselected down-branch motional vector vanished")
    return float((np.vdot(vector, x_operator @ vector).real) / norm_squared)


def compiled_gate_x_trace_fast(
    prep,
    model: dict[str, float],
    steps: int,
    cutoff: int,
) -> tuple[np.ndarray, np.ndarray]:
    dt = model["dt"]
    vartheta = model["dt"] * model["amplitude"]
    varphi = model["varphi"]
    displacement = DisplacementApplicator(cutoff)
    x_operator = position_matrix(cutoff)

    amplitudes = apply_prep(initial_state(cutoff), prep, cutoff).reshape((2, cutoff))
    times_s = np.arange(steps + 1, dtype=np.float64) * dt
    values = np.empty(steps + 1, dtype=np.float64)
    values[0] = schrodinger_down_x_expectation(amplitudes, model, 0.0, x_operator)

    sqr_left = sqr_matrix(-vartheta, -varphi)
    sqr_right = sqr_matrix(-vartheta, varphi)

    for step_index in range(steps):
        phase = model["alpha_phase_offset"] + model["delta"] * step_index * dt
        alpha = model["alpha0"] * complex(math.cos(phase), math.sin(phase))

        amplitudes = apply_xcd_fast(amplitudes, -alpha, displacement)
        amplitudes = sqr_left @ amplitudes
        amplitudes = apply_xcd_fast(amplitudes, alpha, displacement)
        amplitudes = apply_xcd_fast(amplitudes, alpha, displacement)
        amplitudes = sqr_right @ amplitudes
        amplitudes = apply_xcd_fast(amplitudes, -alpha, displacement)

        values[step_index + 1] = schrodinger_down_x_expectation(
            amplitudes,
            model,
            times_s[step_index + 1],
            x_operator,
        )

    return times_s, values


def interpolated_peak(
    frequencies_hz: np.ndarray,
    amplitudes: np.ndarray,
    index: int,
) -> tuple[float, float]:
    if index <= 0 or index >= len(amplitudes) - 1:
        return float(frequencies_hz[index]), float(amplitudes[index])

    y0, y1, y2 = np.log(np.maximum(amplitudes[index - 1 : index + 2], np.finfo(float).tiny))
    denominator = y0 - 2.0 * y1 + y2
    if math.isclose(float(denominator), 0.0, rel_tol=0.0, abs_tol=1e-14):
        return float(frequencies_hz[index]), float(amplitudes[index])

    bin_offset = float(np.clip(0.5 * (y0 - y2) / denominator, -1.0, 1.0))
    frequency_step = frequencies_hz[1] - frequencies_hz[0]
    peak_hz = frequencies_hz[index] + bin_offset * frequency_step
    peak_log_amplitude = y1 - 0.25 * (y0 - y2) * bin_offset
    return float(peak_hz), float(math.exp(peak_log_amplitude))


def fft_spectrum(
    times_s: np.ndarray,
    values: np.ndarray,
    *,
    pad_factor: int,
    min_frequency_hz: float,
    max_frequency_hz: float | None = None,
) -> Spectrum:
    if len(times_s) < 3:
        raise ValueError("at least three samples are required for an FFT peak estimate")

    dt = float(times_s[1] - times_s[0])
    if not np.allclose(np.diff(times_s), dt, rtol=1e-8, atol=1e-14):
        raise ValueError("FFT peak estimate requires uniformly sampled traces")

    centered = np.asarray(values, dtype=np.float64) - float(np.mean(values))
    window = np.hanning(centered.size)
    pad_length = 1 << math.ceil(math.log2(max(centered.size * pad_factor, centered.size)))
    amplitudes = np.abs(np.fft.rfft(centered * window, n=pad_length))
    frequencies_hz = np.fft.rfftfreq(pad_length, dt)

    mask = frequencies_hz >= min_frequency_hz
    if max_frequency_hz is not None:
        mask &= frequencies_hz <= max_frequency_hz
    candidate_indices = np.flatnonzero(mask)
    if candidate_indices.size == 0:
        raise ValueError("frequency search window contains no FFT bins")

    peak_index = int(candidate_indices[np.argmax(amplitudes[candidate_indices])])
    peak_hz, peak_amplitude = interpolated_peak(frequencies_hz, amplitudes, peak_index)
    return Spectrum(
        frequencies_hz=frequencies_hz,
        amplitudes=amplitudes,
        peak_hz=peak_hz,
        peak_amplitude=peak_amplitude,
    )


def vartheta_grid(args: argparse.Namespace) -> list[float]:
    if args.vartheta_values is not None:
        values = [float(value) for value in args.vartheta_values]
    else:
        if args.vartheta_step <= 0.0:
            raise ValueError("--vartheta-step must be positive")
        if args.vartheta_min <= 0.0 or args.vartheta_max <= 0.0:
            raise ValueError("--vartheta-min and --vartheta-max must be positive")
        if args.vartheta_max < args.vartheta_min:
            raise ValueError("--vartheta-max must be at least --vartheta-min")
        count = int(math.floor((args.vartheta_max - args.vartheta_min) / args.vartheta_step + 1e-9))
        values = [args.vartheta_min + index * args.vartheta_step for index in range(count + 1)]
        if values[-1] < args.vartheta_max - 1e-9:
            values.append(args.vartheta_max)

    if any(value <= 0.0 for value in values):
        raise ValueError("all vartheta values must be positive")
    return values


def nearest_selected_indices(values: list[float], selected: list[float]) -> set[int]:
    indices: set[int] = set()
    for target in selected:
        index = min(range(len(values)), key=lambda candidate: abs(values[candidate] - target))
        indices.add(index)
    return indices


def run_sweep(args: argparse.Namespace) -> tuple[list[SweepRow], list[dict[str, object]]]:
    dt_us = jaqal_timestep_us(args.jaqal, args.dt_us)
    base_model, prep = load_program_for_hsim_trace(args.jaqal, dt_us)
    exact = ExactEvolution(prep, base_model, args.cutoff)
    b_rad_s = base_model["amplitude"]
    max_time_s = args.max_time_ms * 1e-3
    values = vartheta_grid(args)
    selected = nearest_selected_indices(values, args.selected_vartheta)

    rows: list[SweepRow] = []
    examples: list[dict[str, object]] = []

    for index, vartheta in enumerate(values):
        dt = vartheta / b_rad_s
        steps = max(2, int(round(max_time_s / dt)))
        actual_time_s = steps * dt
        model = dict(base_model)
        model["dt"] = dt

        times_s, compiled_values = compiled_gate_x_trace_fast(prep, model, steps, args.cutoff)
        exact_values = exact.x_trace(times_s)
        exact_spectrum = fft_spectrum(
            times_s,
            exact_values,
            pad_factor=args.fft_pad_factor,
            min_frequency_hz=args.min_frequency_hz,
        )
        compiled_global_spectrum = fft_spectrum(
            times_s,
            compiled_values,
            pad_factor=args.fft_pad_factor,
            min_frequency_hz=args.min_frequency_hz,
        )
        compiled_min_hz = max(args.min_frequency_hz, exact_spectrum.peak_hz - args.peak_window_hz)
        compiled_max_hz = exact_spectrum.peak_hz + args.peak_window_hz
        compiled_spectrum = fft_spectrum(
            times_s,
            compiled_values,
            pad_factor=args.fft_pad_factor,
            min_frequency_hz=compiled_min_hz,
            max_frequency_hz=compiled_max_hz,
        )

        difference = compiled_values - exact_values
        frequency_shift_hz = compiled_spectrum.peak_hz - exact_spectrum.peak_hz
        row = SweepRow(
            vartheta=vartheta,
            dt_us=dt * 1e6,
            steps=steps,
            actual_time_ms=actual_time_s * 1e3,
            xcd_gates=4 * steps,
            rz_gates=2 * steps,
            total_evolution_gates=6 * steps,
            exact_peak_hz=exact_spectrum.peak_hz,
            compiled_peak_hz=compiled_spectrum.peak_hz,
            compiled_global_peak_hz=compiled_global_spectrum.peak_hz,
            frequency_shift_hz=frequency_shift_hz,
            frequency_shift_rad_s=2.0 * math.pi * frequency_shift_hz,
            relative_frequency_shift=frequency_shift_hz / exact_spectrum.peak_hz,
            rms_x_error=float(np.sqrt(np.mean(difference**2))),
            max_abs_x_error=float(np.max(np.abs(difference))),
        )
        rows.append(row)

        if index in selected:
            examples.append(
                {
                    "row": row,
                    "times_s": times_s,
                    "exact_values": exact_values,
                    "compiled_values": compiled_values,
                    "exact_spectrum": exact_spectrum,
                    "compiled_spectrum": compiled_spectrum,
                }
            )

    return rows, examples


def write_csv(rows: list[SweepRow], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(SweepRow.__dataclass_fields__.keys())
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fieldnames})


def normalized_amplitudes(spectrum: Spectrum) -> np.ndarray:
    maximum = float(np.max(spectrum.amplitudes))
    if maximum <= 0.0:
        return spectrum.amplitudes
    return spectrum.amplitudes / maximum


def plot_spectra(
    examples: list[dict[str, object]],
    output: Path,
    max_frequency_hz: float,
) -> None:
    if not examples:
        return

    figure, axes = plt.subplots(
        len(examples),
        1,
        figsize=(7.2, 2.3 * len(examples)),
        sharex=True,
        constrained_layout=True,
    )
    if len(examples) == 1:
        axes = [axes]

    for axis, example in zip(axes, examples):
        row = example["row"]
        exact_spectrum = example["exact_spectrum"]
        compiled_spectrum = example["compiled_spectrum"]
        assert isinstance(row, SweepRow)
        assert isinstance(exact_spectrum, Spectrum)
        assert isinstance(compiled_spectrum, Spectrum)

        exact_mask = exact_spectrum.frequencies_hz <= max_frequency_hz
        compiled_mask = compiled_spectrum.frequencies_hz <= max_frequency_hz
        axis.plot(
            exact_spectrum.frequencies_hz[exact_mask],
            normalized_amplitudes(exact_spectrum)[exact_mask],
            color="#2f5aa6",
            linewidth=1.8,
            label=r"exact $H_{\mathrm{sim}}$",
        )
        axis.plot(
            compiled_spectrum.frequencies_hz[compiled_mask],
            normalized_amplitudes(compiled_spectrum)[compiled_mask],
            color="#b3453c",
            linewidth=1.4,
            alpha=0.9,
            label="compiled xCD/Rz",
        )
        axis.axvline(exact_spectrum.peak_hz, color="#2f5aa6", linewidth=0.9, alpha=0.5)
        axis.axvline(compiled_spectrum.peak_hz, color="#b3453c", linewidth=0.9, alpha=0.5)
        axis.set_ylabel("normalized FFT")
        axis.set_title(
            rf"$\vartheta={row.vartheta:.3g}$, "
            f"dt={row.dt_us:.2f} us, shift={row.frequency_shift_hz:+.3f} Hz"
        )
        axis.grid(alpha=0.15)

    axes[-1].set_xlabel("frequency (Hz)")
    axes[0].legend(frameon=False, loc="upper right")

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_traces(
    examples: list[dict[str, object]],
    output: Path,
    max_time_ms: float | None,
) -> None:
    if not examples:
        return

    figure, axes = plt.subplots(
        len(examples),
        1,
        figsize=(7.2, 2.2 * len(examples)),
        sharex=True,
        constrained_layout=True,
    )
    if len(examples) == 1:
        axes = [axes]

    for axis, example in zip(axes, examples):
        row = example["row"]
        times_ms = np.asarray(example["times_s"]) * 1e3
        exact_values = np.asarray(example["exact_values"])
        compiled_values = np.asarray(example["compiled_values"])
        assert isinstance(row, SweepRow)

        mask = np.ones_like(times_ms, dtype=bool)
        if max_time_ms is not None:
            mask = times_ms <= max_time_ms

        axis.plot(
            times_ms[mask],
            exact_values[mask],
            color="#2f5aa6",
            linewidth=1.7,
            label=r"exact $H_{\mathrm{sim}}$",
        )
        axis.plot(
            times_ms[mask],
            compiled_values[mask],
            color="#b3453c",
            linewidth=1.2,
            alpha=0.9,
            label="compiled xCD/Rz",
        )
        axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.25)
        axis.set_ylabel(r"$\langle x \rangle$")
        axis.set_title(rf"$\vartheta={row.vartheta:.3g}$, RMS error={row.rms_x_error:.4g}")
        axis.grid(alpha=0.15)

    axes[-1].set_xlabel("time (ms)")
    axes[0].legend(frameon=False, loc="upper right")

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


FREQ_TOML_MAP = {
    "output.jaqal": "jaqal",
    "frequency_sweep.csv": "output_csv",
    "frequency_sweep.spectra_png": "spectra_output",
    "frequency_sweep.trace_png": "trace_output",
    "frequency_sweep.max_time_ms": "max_time_ms",
    "frequency_sweep.vartheta_min": "vartheta_min",
    "frequency_sweep.vartheta_max": "vartheta_max",
    "frequency_sweep.vartheta_step": "vartheta_step",
    "frequency_sweep.selected_vartheta": "selected_vartheta",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep vartheta for the DMANH+ compiled-gate dynamics and estimate "
            "FFT peak shifts relative to exact H_sim propagation."
        )
    )
    add_experiment_arg(parser)
    parser.add_argument("--jaqal", type=Path, default=DEFAULT_JAQAL)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_CSV_OUTPUT)
    parser.add_argument("--spectra-output", type=Path, default=DEFAULT_SPECTRA_OUTPUT)
    parser.add_argument("--trace-output", type=Path, default=DEFAULT_TRACE_OUTPUT)
    parser.add_argument("--cutoff", type=int, default=FOCK_CUTOFF)
    parser.add_argument(
        "--dt-us",
        type=float,
        default=None,
        help=(
            "Base Jaqal timestep in microseconds. Defaults to the generated "
            f"Jaqal timestep comment, or {DEFAULT_DT_US:g} if absent."
        ),
    )
    parser.add_argument("--max-time-ms", type=float, default=160.0)
    parser.add_argument("--vartheta-min", type=float, default=0.1)
    parser.add_argument("--vartheta-max", type=float, default=3.0)
    parser.add_argument("--vartheta-step", type=float, default=0.1)
    parser.add_argument("--vartheta-values", type=float, nargs="+", default=None)
    parser.add_argument("--selected-vartheta", type=float, nargs="+", default=[0.1, 0.8, 1.6, 3.0])
    parser.add_argument("--fft-pad-factor", type=int, default=16)
    parser.add_argument("--min-frequency-hz", type=float, default=1.0)
    parser.add_argument(
        "--peak-window-hz",
        type=float,
        default=150.0,
        help="Track the compiled-gate peak inside this half-width around the exact peak.",
    )
    parser.add_argument("--spectra-max-frequency-hz", type=float, default=500.0)
    parser.add_argument(
        "--trace-max-time-ms",
        type=float,
        default=None,
        help="Optional right edge for the trace-example plot.",
    )
    apply_toml_defaults(
        parser, FREQ_TOML_MAP,
        path_dests={"jaqal", "output_csv", "spectra_output", "trace_output"},
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows, examples = run_sweep(args)
    write_csv(rows, args.output_csv)
    plot_spectra(examples, args.spectra_output, args.spectra_max_frequency_hz)
    plot_traces(examples, args.trace_output, args.trace_max_time_ms)

    print(args.output_csv)
    print(args.spectra_output)
    print(args.trace_output)
    if rows:
        first = rows[0]
        last = rows[-1]
        print(
            "frequency shift range: "
            f"{first.vartheta:g} -> {first.frequency_shift_hz:+.6g} Hz, "
            f"{last.vartheta:g} -> {last.frequency_shift_hz:+.6g} Hz"
        )


if __name__ == "__main__":
    main()
