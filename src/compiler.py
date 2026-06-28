from __future__ import annotations

from __future__ import annotations

import argparse
import cmath
import csv
import math
from pathlib import Path

import numpy as np

from experiment import add_experiment_arg, apply_toml_defaults
from gates import cosine_gate_lines, fmt, zcd_line


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "build" / "double_well.jaqal"
SANDIA_USEPULSES_MODULE = "Calibration_PulseDefinitions.QubitBosonPulses"
DEFAULT_MAX_TIME_MS = 4.0
DEFAULT_DT_US = 200.0
DEFAULT_DELTA_HZ = 500.0


def resolve_evolution_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.b_rad_s is not None:
        if math.isclose(args.b_rad_s, 0.0, rel_tol=0.0, abs_tol=1e-15):
            raise ValueError("--B-rad-s must be nonzero")

        dt_us_from_b = 1.0e6 * args.vartheta / args.b_rad_s
        if dt_us_from_b <= 0.0:
            raise ValueError("--B-rad-s and --vartheta must imply a positive timestep")
        if args.dt_us is not None and not math.isclose(
            args.dt_us,
            dt_us_from_b,
            rel_tol=1e-6,
            abs_tol=1e-3,
        ):
            raise ValueError("--dt-us disagrees with --B-rad-s and --vartheta")
        args.dt_us = dt_us_from_b
    elif args.dt_us is None:
        args.dt_us = DEFAULT_DT_US

    if args.dt_us <= 0.0:
        raise ValueError("--dt-us must be positive")

    if args.delta_rad_s is not None:
        if args.delta_hz is not None:
            delta_rad_s_from_hz = 2.0 * math.pi * args.delta_hz
            if not math.isclose(
                args.delta_rad_s,
                delta_rad_s_from_hz,
                rel_tol=1e-6,
                abs_tol=1e-6,
            ):
                raise ValueError("--delta-rad-s disagrees with --delta-hz")
        args.delta_hz = args.delta_rad_s / (2.0 * math.pi)
    else:
        if args.delta_hz is None:
            args.delta_hz = DEFAULT_DELTA_HZ
        args.delta_rad_s = 2.0 * math.pi * args.delta_hz

    if args.steps is not None:
        if args.steps <= 0:
            raise ValueError("--steps must be positive")
        max_time_ms_from_steps = args.steps * args.dt_us * 1e-3
        if args.max_time_ms is not None and not math.isclose(
            args.max_time_ms,
            max_time_ms_from_steps,
            rel_tol=1e-9,
            abs_tol=1e-6,
        ):
            raise ValueError("--max-time-ms disagrees with --steps and the resolved --dt-us")
        args.max_time_ms = max_time_ms_from_steps
    elif args.max_time_ms is None:
        args.max_time_ms = DEFAULT_MAX_TIME_MS

    if args.max_time_ms <= 0.0:
        raise ValueError("--max-time-ms must be positive")

    return args


def build_evolution_lines(args: argparse.Namespace) -> list[str]:
    max_time_s = args.max_time_ms / 1000.0
    dt_s = args.dt_us * 1e-6
    if getattr(args, "steps", None) is None:
        steps = round(max_time_s / dt_s)
        if not math.isclose(steps * dt_s, max_time_s, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("--max-time-ms must be an integer multiple of --dt-us")
    else:
        steps = args.steps
        max_time_s = steps * dt_s

    delta = getattr(args, "delta_rad_s", 2.0 * math.pi * args.delta_hz)
    amplitude = args.vartheta / dt_s
    prep_beta = args.x_min / math.sqrt(2.0)

    lines: list[str] = [
        "// Symmetric McGarry double-well using Sandia/QSCOUT gate vocabulary.",
        "// Sandia xCD arguments carry direct McGarry SDD displacement coordinates.",
        "// McGarry SQR rotations are decomposed into Sandia Rz gates.",
        "",
        "// Wavepacket preparation: |down>|0> -> |down>|x=-x_min> using zCD.",
        zcd_line(
            args.qubit_index,
            args.sideband_manifold,
            args.sideband_index,
            prep_beta,
            args.nf_start,
            args.nf_end,
        ),
        "",
        "// Symmetric double-well time evolution.",
        "// One timestep is one McGarry cosine gate Gc, expanded into two Q-blocks.",
        f"// alpha phase offset = {fmt(args.alpha_phase_offset)} rad.",
        f"// harmonic delta = {fmt(delta)} rad/s.",
        f"// cosine amplitude B = {fmt(amplitude)} rad/s.",
        f"// steps = {steps}.",
        f"// total evolution time = {fmt(max_time_s * 1000.0)} ms.",
        f"// timestep = {fmt(args.dt_us)} us.",
    ]

    for step in range(steps):
        zeta = args.alpha_phase_offset + delta * step * dt_s
        alpha = args.alpha0 * cmath.exp(1j * zeta)
        lines.append("")
        lines.append(
            f"// evolution step {step + 1}: k = {step}, "
            f"interval {fmt(step * dt_s * 1000.0)} -> {fmt((step + 1) * dt_s * 1000.0)} ms"
        )
        lines.extend(
            cosine_gate_lines(
                args.qubit_index,
                args.sideband_manifold,
                args.sideband_index,
                args.nf_start,
                args.nf_end,
                alpha,
                args.vartheta,
            )
        )

    return lines


def build_angles(args: argparse.Namespace) -> np.ndarray:
    """Return a (2, steps) array of gate parameters per evolution step.

    Row 0: -alpha.real  (first xCD real argument)
    Row 1: -alpha.imag  (first xCD imaginary argument)
    """
    max_time_s = args.max_time_ms / 1000.0
    dt_s = args.dt_us * 1e-6
    if getattr(args, "steps", None) is None:
        steps = round(max_time_s / dt_s)
    else:
        steps = args.steps

    delta = getattr(args, "delta_rad_s", 2.0 * math.pi * args.delta_hz)
    angles = np.empty((2, steps))
    for step in range(steps):
        zeta = args.alpha_phase_offset + delta * step * dt_s
        alpha = args.alpha0 * cmath.exp(1j * zeta)
        angles[0, step] = -alpha.real
        angles[1, step] = -alpha.imag
    return angles


def build_program(args: argparse.Namespace) -> str:
    if args.probe_qubit_index < 0:
        raise ValueError("--probe-qubit-index must be nonnegative")

    if args.probe_qubit_index == args.qubit_index:
        readout_layout = "// Readout uses the same qubit as preparation/evolution, matching McGarry."
    else:
        readout_layout = "// Readout uses a separate probe qubit, matching the QSCOUT notebook pattern."

    register_size = max(args.qubit_index, args.probe_qubit_index) + 1
    lines = [
        f"from {SANDIA_USEPULSES_MODULE} usepulses *",
        "",
        "// Hardware-facing characteristic-function measurement program.",
        "// McGarry methods:meas measures chi(beta) with xCD(beta/2), then measure_all.",
        readout_layout,
        "// Set imMeas=1 to prepend Rx(pi/2), written as R q[probe] 0 pi/2.",
        "// Readout variables are direct xCD pulse coordinates. Initialize at",
        "// the origin; external runners/notebooks should sweep imBeta.",
        "let reBeta 0",
        "let imBeta 0",
        f"let imMeas {args.imaginary_readout_loop}",
        "",
        f"register q[{register_size}]",
        "",
        "prepare_all",
        "",
        *build_evolution_lines(args),
        "",
        "// Characteristic-function readout.",
        "loop imMeas {",
        f"    R q[{args.probe_qubit_index}] 0 {fmt(math.pi / 2.0)}",
        "}",
        xcd_readout_line(
            args.probe_qubit_index,
            args.sideband_manifold,
            args.sideband_index,
            args.nf_start,
            args.nf_end,
        ),
        "measure_all",
    ]
    return "\n".join(lines) + "\n"


def xcd_readout_line(
    qubit_index: int,
    sideband_manifold: int,
    sideband_index: int,
    nf_start: int,
    nf_end: int,
) -> str:
    args = f"{sideband_manifold} {sideband_index} reBeta imBeta"
    if nf_start != 0 or nf_end != 1:
        args = f"{args} {nf_start} {nf_end}"
    return f"xCD q[{qubit_index}] {args}"


# Compiler parameters and their McGarry-paper counterparts.
#
# Local paper source:
#   ../reference/arXiv-2603.04744v1/anharmonicity.tex
#
# Evolution model:
#   --max-time-ms
#       Total simulated time T = K Delta t in eqs. time_trott/state_evo.
#       The default 4 ms is used when --steps is not provided. If --steps is
#       provided, the compiler derives T from K Delta t instead.
#   --steps
#       Integer number of McGarry timesteps K. This is the preferred way to
#       avoid asking for a total time that is not on the timestep grid.
#   --dt-us
#       Simulated timestep Delta t in eqs. time_trott, trig_gate, and
#       Gc_params. When --B-rad-s is provided, the compiler derives
#       Delta t = vartheta / B.
#   --B-rad-s
#       Cosine amplitude B in H_sim/Gc_params. This is not emitted as a
#       separate Jaqal value; it is realized by the pair vartheta and Delta t.
#   --delta-hz
#       delta / (2 pi) in H_sim and Gc_params. The compiler converts this to
#       angular frequency delta and uses it in zeta_k = zeta_0 + k delta
#       Delta t. McGarry's symmetric run uses 500 Hz.
#   --delta-rad-s
#       Direct angular-frequency form of delta. This is useful when the fit
#       reports delta in rad/s, as in Phil's DMANH+ notes.
#   --alpha0
#       alpha_0 in the SDD amplitude alpha = alpha_0 exp(i zeta), eqs. R_phi,
#       Q_def/Q, and Gc_params. For the single-Fourier double well,
#       alpha_0 = pi / (sqrt(2) Lambda). McGarry's symmetric run uses pi / 6.
#   --alpha-phase-offset
#       Compiler phase zeta_0 added to McGarry's zeta = k delta Delta t in
#       Gc_params. With direct xCD displacement coordinates, the default
#       -pi/2 makes the first cosine gate act on x in the local
#       Schrodinger-frame plotter.
#   --vartheta
#       SQR angle vartheta and trigonometric-gate strength, eqs. R_phi,
#       cosine_evo, and Gc_params. It satisfies vartheta = B Delta t. McGarry's
#       symmetric hardware value and the compiler default are 0.8.
#   --x-min
#       x_min in the initialisation method. The compiler prepares the left-well
#       state by emitting D(-sigma_x x_min / sqrt(2)); McGarry uses x_min = 1.5.
#
# Sandia/QSCOUT hardware address parameters:
#   --qubit-index
#       Addressed DV qubit used for McGarry preparation and evolution. The
#       paper is a single-ion experiment and does not assign Jaqal indices.
#   --sideband-manifold, --sideband-index
#       Sandia xCD motional-mode address. These identify the oscillator used
#       for McGarry's SDD gates; they are not TGIFS physics parameters.
#   --nf-start, --nf-end
#       Optional Sandia pulse-calibration Fock-state arguments for xCD. These
#       have no direct McGarry-theory counterpart.
#   --probe-qubit-index
#       Hardware readout qubit for the characteristic-function measurement in
#       methods:meas/fig:circuit. The default uses the same qubit as
#       --qubit-index, matching McGarry's single-ion experiment. Set this to a
#       different index, e.g. 1, for the two-qubit QSCOUT-notebook layout.
#
# Characteristic-function readout:
#   --imaginary-readout-loop
#       Toggle for the optional prepended R(0, pi/2) in the generated readout
#       block. The generated reBeta/imBeta variables are initialized at zero
#       and are intended to be swept externally.
#
# Output-only parameter:
#   --output
#       Filesystem path for the generated Jaqal. This has no McGarry-paper
#       counterpart.
#   --export-angles
#       Filesystem path for the exported notebook angle array. The compiler also
#       writes a CSV sidecar with the same basename. If a directory is supplied
#       for compatibility with older configs, the compiler writes angles.npy and
#       angles.csv inside it.
COMPILER_TOML_MAP = {
    "evolution.steps": "steps",
    "evolution.B_rad_s": "b_rad_s",
    "evolution.delta_rad_s": "delta_rad_s",
    "evolution.delta_hz": "delta_hz",
    "evolution.alpha0": "alpha0",
    "evolution.vartheta": "vartheta",
    "evolution.x_min": "x_min",
    "evolution.alpha_phase_offset": "alpha_phase_offset",
    "evolution.max_time_ms": "max_time_ms",
    "evolution.dt_us": "dt_us",
    "output.jaqal": "output",
    "output.angles": "export_angles",
    "output.notebook": "export_angles",
    "readout.imaginary_loop": "imaginary_readout_loop",
}


def _apply_compiler_toml(parser: argparse.ArgumentParser) -> None:
    apply_toml_defaults(parser, COMPILER_TOML_MAP, path_dests={"output", "export_angles"})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    add_experiment_arg(parser)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--max-time-ms",
        type=float,
        default=None,
        help=(
            f"Total simulated time. Defaults to {DEFAULT_MAX_TIME_MS:g} ms "
            "unless --steps is provided."
        ),
    )
    parser.add_argument("--steps", type=int, default=None, help="Integer McGarry timestep count K.")
    parser.add_argument(
        "--dt-us",
        type=float,
        default=None,
        help=(
            f"Timestep in microseconds. Defaults to {DEFAULT_DT_US:g} unless "
            "--B-rad-s is provided."
        ),
    )
    parser.add_argument(
        "--B-rad-s",
        "--b-rad-s",
        "--cosine-amplitude-rad-s",
        dest="b_rad_s",
        type=float,
        default=None,
        help="Cosine amplitude B in rad/s; derives dt_us = 1e6 * vartheta / B.",
    )
    parser.add_argument(
        "--delta-hz",
        type=float,
        default=None,
        help=f"Harmonic frequency delta / (2 pi). Defaults to {DEFAULT_DELTA_HZ:g} Hz.",
    )
    parser.add_argument("--delta-rad-s", type=float, default=None, help="Angular harmonic frequency delta.")
    parser.add_argument("--alpha0", type=float, default=math.pi / 6.0)
    parser.add_argument("--alpha-phase-offset", type=float, default=-math.pi / 2.0)
    parser.add_argument("--vartheta", type=float, default=0.8)
    parser.add_argument("--x-min", type=float, default=1.5)
    parser.add_argument("--qubit-index", type=int, default=0)
    parser.add_argument("--sideband-manifold", type=int, default=1)
    parser.add_argument("--sideband-index", type=int, default=1)
    parser.add_argument("--nf-start", type=int, default=0)
    parser.add_argument("--nf-end", type=int, default=1)
    parser.add_argument(
        "--probe-qubit-index",
        type=int,
        default=1,
        help=(
            "Probe qubit for hardware readout. Defaults to 1 for the "
            "two-qubit QSCOUT-notebook layout; use 0 for same-qubit McGarry readout."
        ),
    )
    parser.add_argument(
        "--imaginary-readout-loop",
        type=int,
        choices=[0, 1],
        default=1,
        help="Default imMeas loop value: 1 emits R q[probe] 0 pi/2 before xCD readout.",
    )
    parser.add_argument(
        "--export-angles",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Export the notebook angle array to PATH and write a CSV sidecar. "
            "If PATH is a directory, write angles.npy and angles.csv inside it "
            "for compatibility with older configs."
        ),
    )
    parser.add_argument(
        "--export-numpy",
        dest="export_angles",
        type=Path,
        default=None,
        metavar="PATH",
        help=argparse.SUPPRESS,
    )
    _apply_compiler_toml(parser)
    return resolve_evolution_args(parser.parse_args())


def export_angles(args: argparse.Namespace) -> None:
    output = args.export_angles
    if output.suffix == ".npy":
        angles_path = output
        csv_path = output.with_suffix(".csv")
    else:
        angles_path = output / "angles.npy"
        csv_path = output / "angles.csv"
    angles_path.parent.mkdir(parents=True, exist_ok=True)

    angles = build_angles(args)
    np.save(angles_path, angles)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "step",
                "evolution_step",
                "first_xcd_re",
                "first_xcd_im",
            ],
        )
        writer.writeheader()
        for step in range(angles.shape[1]):
            writer.writerow(
                {
                    "step": step,
                    "evolution_step": step + 1,
                    "first_xcd_re": f"{angles[0, step]:.17g}",
                    "first_xcd_im": f"{angles[1, step]:.17g}",
                }
            )

    print(angles_path)
    print(csv_path)


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    text = build_program(args)
    args.output.write_text(text)
    print(args.output)

    if args.export_angles is not None:
        export_angles(args)


if __name__ == "__main__":
    main()
