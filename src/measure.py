from __future__ import annotations

"""McGarry-style ideal measurement readout for the Sandia xCD program.
"""

import argparse
import math
from pathlib import Path

import numpy as np

from plots import (
    DEFAULT_CHI_TITLE,
    DEFAULT_COMPARISON_TITLE,
    plot_characteristic_slice_panels,
    plot_density_panels,
    plot_measurement_comparison,
)
from simulator import (
    DEFAULT_DT_US,
    DEFAULT_JAQAL,
    DEFAULT_TIMES_MS,
    FOCK_CUTOFF,
    REPO_ROOT,
    X_MAX,
    X_MIN,
    X_SAMPLES,
    annihilation_matrix,
    displacement_matrix,
    ho_basis_functions,
    jaqal_timestep_us,
    load_program_model_and_snapshots,
    rz_matrix,
)


DEFAULT_OUTPUT = REPO_ROOT / "build" / "sandia_double_well_measurement_panels.png"
DEFAULT_CHI_OUTPUT = REPO_ROOT / "build" / "sandia_double_well_chi_slice_panels.png"
DEFAULT_MEASUREMENT_TITLE = "McGarry Eq. 33 readout from chi(beta)"


def trapz(values: np.ndarray, samples: np.ndarray, axis: int = -1) -> np.ndarray:
    if hasattr(np, "trapezoid"):
        return np.trapezoid(values, samples, axis=axis)
    return np.trapz(values, samples, axis=axis)


def motional_density_matrix(
    state: np.ndarray,
    cutoff: int,
    readout_state: str,
) -> tuple[np.ndarray, float]:
    amplitudes = state.reshape((2, cutoff))
    if readout_state == "trace":
        rho = amplitudes.T @ amplitudes.conj()
        return rho / np.trace(rho), 1.0

    branch = {"up": 0, "down": 1}[readout_state]
    vector = amplitudes[branch]
    probability = float(np.vdot(vector, vector).real)
    if probability <= 0.0:
        raise ValueError(f"postselection probability vanished for qubit state {readout_state!r}")

    # McGarry et al. Sec. II.B.1-II.B.2 and Appendix C post-select on the
    # desired qubit outcome before oscillator readout. This projection is the
    # ideal-state-vector analogue of discarding failed shots.
    rho = np.outer(vector, vector.conj()) / probability
    return rho, probability


def position_density_from_rho(rho: np.ndarray, basis: np.ndarray) -> np.ndarray:
    density = np.einsum("nx,nm,mx->x", basis, rho, basis, optimize=True)
    return np.maximum(density.real, 0.0)


def characteristic_function(beta: complex, rho: np.ndarray, cutoff: int) -> complex:
    # McGarry et al. Sec. II.B.3 define the measured motional characteristic
    # function as chi(beta) = <D(beta)> with the expectation over the motional
    # state. The hardware estimates Re/Im chi through qubit Z readout; here we
    # evaluate the same expectation value exactly, without projection noise.
    return complex(np.trace(rho @ displacement_matrix(beta, cutoff)))


def probe_initial_vector(state: str) -> np.ndarray:
    if state == "up":
        return np.array([1.0, 0.0], dtype=np.complex128)
    if state == "down":
        return np.array([0.0, 1.0], dtype=np.complex128)
    raise ValueError(f"unsupported probe initial state {state!r}")


def conditional_displacement_unitary(beta: complex, cutoff: int) -> np.ndarray:
    # McGarry methods:meas reads chi(beta)=<D(beta)> by applying
    # D(sigma_x beta/2). Sandia xCD's direct displacement argument is therefore
    # beta/2 for this readout block.
    sandia_xcd_beta = 0.5 * beta
    plus = np.array([1.0, 1.0], dtype=np.complex128) / math.sqrt(2.0)
    minus = np.array([1.0, -1.0], dtype=np.complex128) / math.sqrt(2.0)
    plus_projector = np.outer(plus, plus.conj())
    minus_projector = np.outer(minus, minus.conj())
    return np.kron(plus_projector, displacement_matrix(sandia_xcd_beta, cutoff)) + np.kron(
        minus_projector,
        displacement_matrix(-sandia_xcd_beta, cutoff),
    )


def probe_z_expectation(
    rho: np.ndarray,
    cutoff: int,
    beta: complex,
    use_rz_rotation: bool,
    rotation_angle: float,
    probe_initial_state: str,
) -> float:
    probe_vector = probe_initial_vector(probe_initial_state)
    if use_rz_rotation and not math.isclose(
        rotation_angle,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-14,
    ):
        probe_vector = rz_matrix(rotation_angle) @ probe_vector

    probe_rho = np.outer(probe_vector, probe_vector.conj())
    joint_rho = np.kron(probe_rho, rho)
    readout_unitary = conditional_displacement_unitary(beta, cutoff)
    measured_rho = readout_unitary @ joint_rho @ readout_unitary.conj().T

    sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    observable = np.kron(sigma_z, np.eye(cutoff, dtype=np.complex128))
    return float(np.trace(measured_rho @ observable).real)


def probe_readout_check(
    rho: np.ndarray,
    cutoff: int,
    beta: complex,
    probe_initial_state: str,
) -> dict[str, float | complex]:
    # This explicitly emulates the ideal version of the generated readout block
    #   optional Rz(pi/2); xCD(beta/2); measure Z
    # and compares it to the direct McGarry characteristic function. The two
    # sign-error columns keep the Rz/qubit-basis convention visible.
    chi = characteristic_function(beta, rho, cutoff)
    z_no_rotation = probe_z_expectation(
        rho,
        cutoff,
        beta,
        use_rz_rotation=False,
        rotation_angle=0.0,
        probe_initial_state=probe_initial_state,
    )
    z_rz_pi_2 = probe_z_expectation(
        rho,
        cutoff,
        beta,
        use_rz_rotation=True,
        rotation_angle=math.pi / 2.0,
        probe_initial_state=probe_initial_state,
    )
    return {
        "probe_beta": beta,
        "probe_chi": chi,
        "probe_z_no_rotation": z_no_rotation,
        "probe_z_rz_pi_2": z_rz_pi_2,
        "probe_re_abs_error": abs(z_no_rotation - chi.real),
        "probe_rz_pi_2_minus_im_abs_error": abs(z_rz_pi_2 - chi.imag),
        "probe_rz_pi_2_plus_im_abs_error": abs(z_rz_pi_2 + chi.imag),
    }


def position_expectation_from_rho(rho: np.ndarray, cutoff: int) -> float:
    annihilation = annihilation_matrix(cutoff)
    position = (annihilation + annihilation.conj().T) / math.sqrt(2.0)
    return float(np.trace(rho @ position).real)


def position_expectation_2pfd(
    rho: np.ndarray,
    cutoff: int,
    h: float,
) -> tuple[float, complex]:
    if h <= 0.0:
        raise ValueError("--pfd-h must be positive")

    # McGarry et al. Eq. xexpect gives
    #   <x> = (1/sqrt(2)) d Im[chi(beta)] / d Im[beta] at beta=0.
    # Eq. xapprox estimates the slope using the origin and chi(i h).
    chi_ih = characteristic_function(1.0j * h, rho, cutoff)
    return float(chi_ih.imag / (math.sqrt(2.0) * h)), chi_ih


def positive_imaginary_beta_grid(beta_imag_max: float, positive_points: int) -> np.ndarray:
    if beta_imag_max <= 0.0:
        raise ValueError("--beta-imag-max must be positive")
    if positive_points < 2:
        raise ValueError("--positive-beta-points must be at least 2")
    return np.linspace(0.0, beta_imag_max, positive_points)


def characteristic_slice_imaginary_beta(
    rho: np.ndarray,
    cutoff: int,
    beta_imag_max: float,
    positive_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    positive_beta_imag = positive_imaginary_beta_grid(beta_imag_max, positive_points)
    positive_chi = np.array(
        [characteristic_function(1.0j * beta_imag, rho, cutoff) for beta_imag in positive_beta_imag],
        dtype=np.complex128,
    )

    # Half-plane scans are filled by Hermitian symmetry.
    # For a displacement characteristic function, D(beta)^dagger = D(-beta),
    # so a Hermitian rho gives chi(-beta)=chi(beta)^*. This is the relation
    # used here to fill the unmeasured negative imaginary beta axis.
    beta_imag = np.concatenate((-positive_beta_imag[:0:-1], positive_beta_imag))
    chi = np.concatenate((positive_chi[:0:-1].conj(), positive_chi))
    return beta_imag, chi


def position_density_from_characteristic_slice(
    xs: np.ndarray,
    beta_imag: np.ndarray,
    chi: np.ndarray,
) -> np.ndarray:
    # McGarry et al. Eq. 33:
    #   P(x) = (1 / 2 pi) integral du exp(-i u x) chi(i u / sqrt(2)).
    # Since beta = i u / sqrt(2), the imaginary beta coordinate is u/sqrt(2).
    u = math.sqrt(2.0) * beta_imag
    kernel = np.exp(-1.0j * np.outer(xs, u))
    density = trapz(kernel * chi[np.newaxis, :], u, axis=1) / (2.0 * math.pi)
    return np.maximum(density.real, 0.0)


def normalize_density(xs: np.ndarray, density: np.ndarray) -> np.ndarray:
    area = float(trapz(density, xs))
    if area <= 0.0:
        return density
    return density / area


def l2_relative_error(reference: np.ndarray, candidate: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    if denominator == 0.0:
        return float(np.linalg.norm(candidate))
    return float(np.linalg.norm(candidate - reference) / denominator)


def reconstruct_snapshot_density(
    state: np.ndarray,
    cutoff: int,
    basis: np.ndarray,
    xs: np.ndarray,
    readout_state: str,
    beta_imag_max: float,
    positive_beta_points: int,
    pfd_h: float,
    probe_initial_state: str,
    renormalize: bool,
) -> dict[str, object]:
    rho, postselection_probability = motional_density_matrix(state, cutoff, readout_state)
    direct_density = position_density_from_rho(rho, basis)
    x_direct = position_expectation_from_rho(rho, cutoff)
    x_2pfd, chi_ih = position_expectation_2pfd(rho, cutoff, pfd_h)
    readout_check = probe_readout_check(
        rho,
        cutoff,
        1.0j * pfd_h,
        probe_initial_state,
    )
    beta_imag, chi = characteristic_slice_imaginary_beta(
        rho,
        cutoff,
        beta_imag_max,
        positive_beta_points,
    )
    measured_density = position_density_from_characteristic_slice(xs, beta_imag, chi)

    if renormalize:
        direct_density = normalize_density(xs, direct_density)
        measured_density = normalize_density(xs, measured_density)

    chi0_error = abs(chi[positive_beta_points - 1] - 1.0)
    return {
        "rho": rho,
        "postselection_probability": postselection_probability,
        "beta_imag": beta_imag,
        "chi": chi,
        "direct_density": direct_density,
        "measured_density": measured_density,
        "x_direct": x_direct,
        "x_2pfd": x_2pfd,
        "x_2pfd_abs_error": abs(x_2pfd - x_direct),
        "chi_ih": chi_ih,
        **readout_check,
        "chi0_error": float(chi0_error),
        "max_abs_error": float(np.max(np.abs(measured_density - direct_density))),
        "relative_l2_error": l2_relative_error(direct_density, measured_density),
    }


def reconstruct_snapshots(
    snapshots: dict[float, np.ndarray],
    cutoff: int,
    times_ms: list[float],
    readout_state: str,
    beta_imag_max: float,
    positive_beta_points: int,
    pfd_h: float,
    probe_initial_state: str,
    renormalize: bool,
) -> dict[float, dict[str, object]]:
    xs = np.linspace(X_MIN, X_MAX, X_SAMPLES)
    basis = ho_basis_functions(cutoff, xs)
    return {
        time_ms: reconstruct_snapshot_density(
            snapshots[time_ms],
            cutoff,
            basis,
            xs,
            readout_state,
            beta_imag_max,
            positive_beta_points,
            pfd_h,
            probe_initial_state,
            renormalize=renormalize,
        )
        for time_ms in times_ms
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct McGarry-style position densities from the motional "
            "characteristic function, using the ideal xCD truncated-Fock simulator."
        )
    )
    parser.add_argument("--jaqal", type=Path, default=DEFAULT_JAQAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--title", default=DEFAULT_MEASUREMENT_TITLE)
    parser.add_argument(
        "--chi-output",
        type=Path,
        default=DEFAULT_CHI_OUTPUT,
        help=(
            "McGarry Fig. 7-style Re/Im chi(i y) slice plot. "
            "Use --no-chi-output to skip it."
        ),
    )
    parser.add_argument(
        "--no-chi-output",
        action="store_true",
        help="Do not emit the characteristic-function slice plot.",
    )
    parser.add_argument("--chi-title", default=DEFAULT_CHI_TITLE)
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=None,
        help="Optional direct-vs-readout comparison plot. Omitted by default.",
    )
    parser.add_argument("--comparison-title", default=DEFAULT_COMPARISON_TITLE)
    parser.add_argument("--cutoff", type=int, default=FOCK_CUTOFF)
    parser.add_argument("--times-ms", type=float, nargs="+", default=DEFAULT_TIMES_MS)
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
        "--readout-state",
        choices=["down", "up", "trace"],
        default="down",
        help=(
            "Motional state used for readout. 'down' follows the McGarry "
            "postselected qubit branch; 'trace' reproduces the older oracle "
            "plotter's unconditional reduced oscillator state."
        ),
    )
    parser.add_argument(
        "--beta-imag-max",
        type=float,
        default=5.0,
        help="Maximum Im[beta] for the one-dimensional chi scan. McGarry Fig. 7 uses 5.",
    )
    parser.add_argument(
        "--positive-beta-points",
        type=int,
        default=50,
        help=(
            "Number of nonnegative Im[beta] samples, including zero. McGarry "
            "Fig. 7 uses 50; increase for a denser ideal-readout comparison."
        ),
    )
    parser.add_argument(
        "--pfd-h",
        type=float,
        default=0.4,
        help="Two-point finite-difference offset h for McGarry Eq. xapprox.",
    )
    parser.add_argument(
        "--probe-initial-state",
        choices=["up", "down"],
        default="up",
        help=(
            "Initial probe qubit state for the ideal readout-block check. "
            "'up' makes no-rotation xCD readout match Re[chi] with the local "
            "simulator's Z convention."
        ),
    )
    parser.add_argument(
        "--no-renormalize-density",
        action="store_true",
        help="Do not renormalize direct and reconstructed P(x) over the plotted x window.",
    )
    parser.add_argument(
        "--interaction-frame",
        action="store_true",
        help="Use raw interaction-frame gate state instead of applying the inferred harmonic frame.",
    )
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
    results = reconstruct_snapshots(
        snapshots,
        args.cutoff,
        args.times_ms,
        args.readout_state,
        args.beta_imag_max,
        args.positive_beta_points,
        args.pfd_h,
        args.probe_initial_state,
        renormalize=not args.no_renormalize_density,
    )

    measurement_densities = {
        time_ms: results[time_ms]["measured_density"]
        for time_ms in args.times_ms
    }
    plot_density_panels(
        measurement_densities,
        model,
        args.times_ms,
        args.output,
        title=args.title,
    )
    if not args.no_chi_output:
        plot_characteristic_slice_panels(
            results,
            args.times_ms,
            args.chi_output,
            args.chi_title,
        )
    if args.comparison_output is not None:
        plot_measurement_comparison(
            results,
            model,
            args.times_ms,
            args.comparison_output,
            args.comparison_title,
        )

    print(args.output)
    if not args.no_chi_output:
        print(args.chi_output)
    if args.comparison_output is not None:
        print(args.comparison_output)
    print(
        "readout_state,time_ms,postselection_probability,chi0_error,"
        "max_abs_error,relative_l2_error,x_direct,x_2pfd,x_2pfd_abs_error,chi_ih_imag,"
        "probe_z_no_rotation,probe_re_abs_error,probe_z_rz_pi_2,"
        "probe_rz_pi_2_minus_im_abs_error,probe_rz_pi_2_plus_im_abs_error"
    )
    for time_ms in args.times_ms:
        result = results[time_ms]
        chi_ih = result["chi_ih"]
        print(
            f"{args.readout_state},"
            f"{time_ms:.12g},"
            f"{result['postselection_probability']:.12g},"
            f"{result['chi0_error']:.12g},"
            f"{result['max_abs_error']:.12g},"
            f"{result['relative_l2_error']:.12g},"
            f"{result['x_direct']:.12g},"
            f"{result['x_2pfd']:.12g},"
            f"{result['x_2pfd_abs_error']:.12g},"
            f"{chi_ih.imag:.12g},"
            f"{result['probe_z_no_rotation']:.12g},"
            f"{result['probe_re_abs_error']:.12g},"
            f"{result['probe_z_rz_pi_2']:.12g},"
            f"{result['probe_rz_pi_2_minus_im_abs_error']:.12g},"
            f"{result['probe_rz_pi_2_plus_im_abs_error']:.12g}"
        )


if __name__ == "__main__":
    main()
