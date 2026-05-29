from __future__ import annotations

import math


def fmt(value: float) -> str:
    if math.isclose(value, 0.0, rel_tol=0.0, abs_tol=1e-14):
        return "0"
    return f"{value:.17g}"


def gate_line(name: str, qubit_index: int, *args: float | int) -> str:
    rendered_args = " ".join(fmt(float(arg)) for arg in args)
    if rendered_args:
        return f"{name} q[{qubit_index}] {rendered_args}"
    return f"{name} q[{qubit_index}]"


def compact_rotation_angle(angle_rad: float) -> float:
    angle = (float(angle_rad) + math.pi) % (2.0 * math.pi) - math.pi
    if math.isclose(angle + math.pi, 0.0, rel_tol=0.0, abs_tol=1e-14):
        return math.pi
    return angle


def rz_lines(qubit_index: int, angle_rad: float) -> list[str]:
    angle_rad = compact_rotation_angle(angle_rad)
    if math.isclose(angle_rad, 0.0, rel_tol=0.0, abs_tol=1e-14):
        return []
    return [gate_line("Rz", qubit_index, angle_rad)]


def sqr_lines(qubit_index: int, vartheta: float, varphi: float) -> list[str]:
    if math.isclose(vartheta, 0.0, rel_tol=0.0, abs_tol=1e-14):
        return []
    return [
        *rz_lines(qubit_index, vartheta),
    ]


def xcd_line(
    qubit_index: int,
    sideband_manifold: int,
    sideband_index: int,
    alpha: complex,
    nf_start: int,
    nf_end: int,
) -> str:
    return xcd_beta_line(
        qubit_index,
        sideband_manifold,
        sideband_index,
        alpha,
        nf_start,
        nf_end,
    )


def zcd_line(
    qubit_index: int,
    sideband_manifold: int,
    sideband_index: int,
    beta: complex,
    nf_start: int,
    nf_end: int,
) -> str:
    args: list[float | int] = [sideband_manifold, sideband_index, beta.real, beta.imag]
    if nf_start != 0 or nf_end != 1:
        args.extend([nf_start, nf_end])
    return gate_line("zCD", qubit_index, *args)


def xcd_beta_line(
    qubit_index: int,
    sideband_manifold: int,
    sideband_index: int,
    beta: complex,
    nf_start: int,
    nf_end: int,
) -> str:
    """Emit Sandia xCD using its direct beta pulse argument."""
    args: list[float | int] = [sideband_manifold, sideband_index, beta.real, beta.imag]
    if nf_start != 0 or nf_end != 1:
        args.extend([nf_start, nf_end])
    return gate_line("xCD", qubit_index, *args)


def q_block_lines(
    qubit_index: int,
    sideband_manifold: int,
    sideband_index: int,
    nf_start: int,
    nf_end: int,
    alpha: complex,
    vartheta: float,
    varphi: float,
) -> list[str]:
    # McGarry Eq. 12 defines Q from two SDDs around one SQR.
    return [
        xcd_line(qubit_index, sideband_manifold, sideband_index, alpha, nf_start, nf_end),
        *sqr_lines(qubit_index, -vartheta, varphi),
        xcd_line(qubit_index, sideband_manifold, sideband_index, -alpha, nf_start, nf_end),
    ]


def cosine_gate_lines(
    qubit_index: int,
    sideband_manifold: int,
    sideband_index: int,
    nf_start: int,
    nf_end: int,
    alpha: complex,
    vartheta: float,
    varphi: float,
) -> list[str]:
    # McGarry Eq. 15 writes Gc = Q(alpha, theta, phi) Q(-alpha, theta, -phi).
    # Jaqal source is execution order, so the rightmost Q(-alpha, theta, -phi)
    # is listed first.
    return [
        *q_block_lines(
            qubit_index,
            sideband_manifold,
            sideband_index,
            nf_start,
            nf_end,
            -alpha,
            vartheta,
            -varphi,
        ),
        *q_block_lines(
            qubit_index,
            sideband_manifold,
            sideband_index,
            nf_start,
            nf_end,
            alpha,
            vartheta,
            varphi,
        ),
    ]
