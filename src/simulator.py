from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from jaqalpaq.core.block import BlockStatement, LoopStatement
from jaqalpaq.core.gate import GateStatement
from jaqalpaq.parser import parse_jaqal_file
from scipy.linalg import expm


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JAQAL = (
    REPO_ROOT
    / "build"
    / "double_well.jaqal"
)
DEFAULT_TIMES_MS = [0.0, 2.0, 4.0]
DEFAULT_DT_US = 20.0
FOCK_CUTOFF = 32
X_MIN = -4.0
X_MAX = 4.0
X_SAMPLES = 161
TIMESTEP_COMMENT_RE = re.compile(r"^\s*//\s*timestep\s*=\s*([0-9.eE+-]+)\s*us\.?\s*$")


@dataclass(frozen=True)
class StepParameters:
    alpha: complex
    vartheta: float
    varphi: float
    dt: float


def flatten_statements(statement) -> list[GateStatement]:
    if isinstance(statement, GateStatement):
        return [statement]
    if isinstance(statement, BlockStatement):
        if statement.parallel:
            raise ValueError("parallel Jaqal blocks are not supported by this CV-DV interpreter")
        gates: list[GateStatement] = []
        for child in statement.statements:
            gates.extend(flatten_statements(child))
        return gates
    if isinstance(statement, LoopStatement):
        gates = flatten_statements(statement.statements)
        return gates * int(statement.iterations)
    raise ValueError(f"unsupported Jaqal statement {statement!r}")


def gate_program(jaqal_file: Path) -> list[GateStatement]:
    circuit = parse_jaqal_file(str(jaqal_file), autoload_pulses=False, expand_macro=True)
    return flatten_statements(circuit.body)


def jaqal_timestep_us(jaqal_file: Path, override_dt_us: float | None = None) -> float:
    if override_dt_us is not None:
        return override_dt_us

    for line in jaqal_file.read_text().splitlines():
        match = TIMESTEP_COMMENT_RE.match(line)
        if match is not None:
            return float(match.group(1))

    return DEFAULT_DT_US


def gate_name(gate: GateStatement) -> str:
    return str(gate.name)


def gate_args(gate: GateStatement) -> list[object]:
    return gate.parameters_linear


def cd_beta(gate: GateStatement) -> complex:
    args = gate_args(gate)
    return complex(float(args[3]), float(args[4]))


def compact_rotation_angle(angle_rad: float) -> float:
    angle = (float(angle_rad) + math.pi) % (2.0 * math.pi) - math.pi
    if math.isclose(angle + math.pi, 0.0, rel_tol=0.0, abs_tol=1e-14):
        return math.pi
    return angle


def rz_angle(gate: GateStatement) -> float:
    args = gate_args(gate)
    return compact_rotation_angle(float(args[1]))


def rotation_axis_angle(gate: GateStatement) -> tuple[float, float]:
    args = gate_args(gate)
    return float(args[1]), compact_rotation_angle(float(args[2]))


def r_matrix(phase: float, angle: float) -> np.ndarray:
    sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    axis = math.cos(phase) * sigma_x + math.sin(phase) * sigma_y
    return expm(-0.5j * angle * axis)

def rz_matrix(angle: float) -> np.ndarray:
    return np.array(
        [
            [np.exp(-0.5j * angle), 0.0],
            [0.0, np.exp(0.5j * angle)],
        ],
        dtype=np.complex128,
    )


def infer_sqr_from_rz_gate(gate: GateStatement) -> tuple[float, float]:
    if gate_name(gate) != "Rz":
        raise ValueError(f"expected Rz for symmetric SQR, got {gate_name(gate)}")
    return rz_angle(gate), 0.0


PREP_PATTERNS = (
    ["Rz", "xCD", "Rz"],
    ["zCD"],
)
STEP_PATTERN = ["xCD", "Rz", "xCD", "xCD", "Rz", "xCD"]

def gate_names(gates: list[GateStatement]) -> list[str]:
    return [gate_name(gate) for gate in gates]


def step_blocks_from(gates: list[GateStatement], start: int) -> list[list[GateStatement]]:
    steps: list[list[GateStatement]] = []
    index = start
    while index + len(STEP_PATTERN) <= len(gates):
        candidate = gates[index : index + len(STEP_PATTERN)]
        if gate_names(candidate) != STEP_PATTERN:
            break
        steps.append(candidate)
        index += len(STEP_PATTERN)
    return steps


def split_program(gates: list[GateStatement]) -> tuple[list[GateStatement], list[list[GateStatement]]]:
    # The hardware-facing Jaqal includes prepare_all, readout xCD/Rz gates, and
    # measure_all around the McGarry evolution sequence. For ideal simulation,
    # locate the Rz/xCD/Rz prep followed by whole Gc blocks and ignore the
    # hardware-only wrapper statements.
    for start in range(len(gates)):
        for prep_pattern in PREP_PATTERNS:
            prep_end = start + len(prep_pattern)
            prep = gates[start:prep_end]
            if gate_names(prep) != prep_pattern:
                continue

            steps = step_blocks_from(gates, prep_end)
            if steps:
                return prep, steps

    raise ValueError(
        "could not find McGarry preparation followed by Sandia xCD/Rz Gc blocks "
        "inside the Jaqal program"
    )


def infer_step_parameters(steps: list[list[GateStatement]], dt_s: float) -> list[StepParameters]:
    parameters = []
    for step in steps:
        alpha = cd_beta(step[3])
        sqr_angle, varphi = infer_sqr_from_rz_gate(step[4])
        parameters.append(
            StepParameters(
                alpha=alpha,
                vartheta=-sqr_angle,
                varphi=varphi,
                dt=dt_s,
            )
        )
    return parameters


def infer_double_well(parameters: list[StepParameters]) -> dict[str, float]:
    first = parameters[0]
    alpha0 = abs(first.alpha)
    dt = first.dt
    vartheta = first.vartheta
    varphi = first.varphi
    lambda_period = math.pi / (math.sqrt(2.0) * alpha0)
    amplitude = vartheta / dt

    phases = np.unwrap(np.array([math.atan2(p.alpha.imag, p.alpha.real) for p in parameters]))
    alpha_phase_offset = float(phases[0])
    if len(phases) > 1:
        delta = float(np.mean(np.diff(phases) / np.array([p.dt for p in parameters[:-1]])))
    else:
        delta = 0.0

    return {
        "alpha0": alpha0,
        "dt": dt,
        "delta": delta,
        "amplitude": amplitude,
        "lambda_period": lambda_period,
        "varphi": varphi,
        "alpha_phase_offset": alpha_phase_offset,
    }


def annihilation_matrix(cutoff: int) -> np.ndarray:
    matrix = np.zeros((cutoff, cutoff), dtype=np.complex128)
    for level in range(1, cutoff):
        matrix[level - 1, level] = math.sqrt(level)
    return matrix


def position_matrix(cutoff: int) -> np.ndarray:
    annihilation = annihilation_matrix(cutoff)
    return (annihilation + annihilation.conj().T) / math.sqrt(2.0)


def ho_basis_functions(cutoff: int, xs: np.ndarray) -> np.ndarray:
    basis = np.zeros((cutoff, xs.size), dtype=np.float64)
    basis[0] = np.pi ** (-0.25) * np.exp(-(xs**2) / 2.0)
    if cutoff == 1:
        return basis

    basis[1] = math.sqrt(2.0) * xs * basis[0]
    for level in range(1, cutoff - 1):
        basis[level + 1] = (
            math.sqrt(2.0 / (level + 1.0)) * xs * basis[level]
            - math.sqrt(level / (level + 1.0)) * basis[level - 1]
        )
    return basis


def displacement_matrix(alpha: complex, cutoff: int) -> np.ndarray:
    annihilation = annihilation_matrix(cutoff)
    creation = annihilation.conj().T
    return expm(alpha * creation - alpha.conjugate() * annihilation)


def cosine_of_position(cutoff: int, wave_number: float, phase: float) -> np.ndarray:
    x_operator = position_matrix(cutoff)
    eigenvalues, eigenvectors = np.linalg.eigh(x_operator)
    cos_values = np.cos(wave_number * eigenvalues + phase)
    return (eigenvectors * cos_values[np.newaxis, :]) @ eigenvectors.conj().T


def hsim_hamiltonian(model: dict[str, float], cutoff: int) -> np.ndarray:
    levels = np.arange(cutoff, dtype=np.float64)
    harmonic = np.diag(model["delta"] * (levels + 0.5)).astype(np.complex128)
    wave_number = 2.0 * math.pi / model["lambda_period"]
    anharmonic = model["amplitude"] * cosine_of_position(cutoff, wave_number, model["varphi"])
    return harmonic + anharmonic


def initial_state(cutoff: int) -> np.ndarray:
    qubit_down = np.array([0.0, 1.0], dtype=np.complex128)
    oscillator_ground = np.zeros(cutoff, dtype=np.complex128)
    oscillator_ground[0] = 1.0
    return np.kron(qubit_down, oscillator_ground)


def postselected_motional_vector(
    state: np.ndarray,
    cutoff: int,
    readout_state: str = "down",
) -> np.ndarray:
    amplitudes = state.reshape((2, cutoff))
    branch = {"up": 0, "down": 1}[readout_state]
    vector = amplitudes[branch]
    norm = np.linalg.norm(vector)
    if norm <= 0.0:
        raise ValueError(f"postselected motional vector vanished for qubit state {readout_state!r}")
    return vector / norm


def motional_x_expectation(vector: np.ndarray, cutoff: int) -> float:
    return float(np.vdot(vector, position_matrix(cutoff) @ vector).real)


def postselected_x_expectation(
    state: np.ndarray,
    cutoff: int,
    readout_state: str = "down",
) -> float:
    return motional_x_expectation(
        postselected_motional_vector(state, cutoff, readout_state),
        cutoff,
    )


def apply_xcd(state: np.ndarray, sandia_alpha: complex, cutoff: int) -> np.ndarray:
    amplitudes = state.reshape((2, cutoff))
    z_to_x = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128) / math.sqrt(2.0)
    x_amplitudes = z_to_x.conj().T @ amplitudes
    # Current Sandia/DMANH convention: xCD(s) realizes mathematical D(-i s).
    # Previous local mathematical convention:
    # math_alpha = sandia_alpha
    math_alpha = -1.0j * sandia_alpha
    x_amplitudes[0] = displacement_matrix(math_alpha, cutoff) @ x_amplitudes[0]
    x_amplitudes[1] = displacement_matrix(-math_alpha, cutoff) @ x_amplitudes[1]
    return (z_to_x @ x_amplitudes).reshape(2 * cutoff)

def apply_zcd(state: np.ndarray, beta: complex, cutoff: int) -> np.ndarray:
    amplitudes = state.reshape((2, cutoff)).copy()
    amplitudes[0] = displacement_matrix(beta, cutoff) @ amplitudes[0]
    amplitudes[1] = displacement_matrix(-beta, cutoff) @ amplitudes[1]
    return amplitudes.reshape(2 * cutoff)


def apply_r(state: np.ndarray, phase: float, angle: float, cutoff: int) -> np.ndarray:
    return (r_matrix(phase, angle) @ state.reshape((2, cutoff))).reshape(2 * cutoff)


def apply_rz(state: np.ndarray, angle: float, cutoff: int) -> np.ndarray:
    return (rz_matrix(angle) @ state.reshape((2, cutoff))).reshape(2 * cutoff)


def apply_sqr_phi(state: np.ndarray, vartheta: float, varphi: float, cutoff: int) -> np.ndarray:
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    axis = math.sin(varphi) * sigma_y + math.cos(varphi) * sigma_z
    qubit_unitary = expm(-0.5j * vartheta * axis)
    return (qubit_unitary @ state.reshape((2, cutoff))).reshape(2 * cutoff)


def apply_gate(state: np.ndarray, gate: GateStatement, cutoff: int) -> np.ndarray:
    if gate_name(gate) == "xCD":
        return apply_xcd(state, cd_beta(gate), cutoff)
    if gate_name(gate) == "zCD":
        return apply_zcd(state, cd_beta(gate), cutoff)
    if gate_name(gate) == "Rz":
        return apply_rz(state, rz_angle(gate), cutoff)
    if gate_name(gate) == "R":
        phase, angle = rotation_axis_angle(gate)
        return apply_r(state, phase, angle, cutoff)
    raise ValueError(f"unsupported gate for Sandia CV-DV semantic simulator: {gate_name(gate)}")


def apply_prep(state: np.ndarray, prep: list[GateStatement], cutoff: int) -> np.ndarray:
    if len(prep) in (1, 3):
        for gate in prep:
            state = apply_gate(state, gate, cutoff)
        return state
    raise ValueError(f"unsupported preparation length: {len(prep)}")


def apply_step(state: np.ndarray, step: list[GateStatement], cutoff: int) -> np.ndarray:
    state = apply_gate(state, step[0], cutoff)  # xCD
    state = apply_gate(state, step[1], cutoff)  # Rz
    state = apply_gate(state, step[2], cutoff)  # xCD

    state = apply_gate(state, step[3], cutoff)  # xCD
    state = apply_gate(state, step[4], cutoff)  # Rz
    state = apply_gate(state, step[5], cutoff)  # xCD
    return state


def apply_q_block_from_parameters(
    state: np.ndarray,
    alpha: complex,
    vartheta: float,
    varphi: float,
    cutoff: int,
) -> np.ndarray:
    state = apply_xcd(state, alpha, cutoff)
    state = apply_sqr_phi(state, -vartheta, varphi, cutoff)
    state = apply_xcd(state, -alpha, cutoff)
    return state


def apply_cosine_step_from_parameters(
    state: np.ndarray,
    alpha: complex,
    vartheta: float,
    varphi: float,
    cutoff: int,
) -> np.ndarray:
    state = apply_q_block_from_parameters(state, -alpha, vartheta, -varphi, cutoff)
    state = apply_q_block_from_parameters(state, alpha, vartheta, varphi, cutoff)
    return state


def apply_harmonic_frame(state: np.ndarray, delta: float, time_s: float, cutoff: int) -> np.ndarray:
    levels = np.arange(cutoff, dtype=np.float64)
    oscillator_phases = np.exp(-1j * delta * (levels + 0.5) * time_s)
    amplitudes = state.reshape((2, cutoff)).copy()
    amplitudes *= oscillator_phases[np.newaxis, :]
    return amplitudes.reshape(2 * cutoff)


def oscillator_density(state: np.ndarray, cutoff: int, basis: np.ndarray) -> np.ndarray:
    amplitudes = state.reshape((2, cutoff))
    rho = amplitudes.T @ amplitudes.conj()
    density = np.einsum("nx,nm,mx->x", basis, rho, basis, optimize=True)
    return np.maximum(density.real, 0.0)


def run_snapshots(
    prep: list[GateStatement],
    steps: list[list[GateStatement]],
    parameters: list[StepParameters],
    times_ms: list[float],
    cutoff: int,
    use_schrodinger_frame: bool,
) -> dict[float, np.ndarray]:
    requested = [time_ms / 1000.0 for time_ms in times_ms]
    state = initial_state(cutoff)
    state = apply_prep(state, prep, cutoff)

    model = infer_double_well(parameters)
    snapshots = {}
    elapsed = 0.0
    for time_s, time_ms in zip(requested, times_ms):
        if math.isclose(time_s, 0.0, abs_tol=1e-12):
            snapshots[time_ms] = state.copy()

    for step, step_parameters in zip(steps, parameters):
        state = apply_step(state, step, cutoff)
        elapsed += step_parameters.dt

        for time_s, time_ms in zip(requested, times_ms):
            if time_ms in snapshots:
                continue
            if math.isclose(elapsed, time_s, rel_tol=0.0, abs_tol=step_parameters.dt / 10.0):
                snapshot = state.copy()
                if use_schrodinger_frame:
                    snapshot = apply_harmonic_frame(snapshot, model["delta"], elapsed, cutoff)
                snapshots[time_ms] = snapshot

    missing = [time_ms for time_ms in times_ms if time_ms not in snapshots]
    if missing:
        raise ValueError(f"program did not reach requested snapshot times: {missing}")
    return snapshots


def exact_hsim_x_trace(
    prep: list[GateStatement],
    model: dict[str, float],
    times_ms: np.ndarray,
    cutoff: int,
) -> np.ndarray:
    # McGarry Fig. 6's first trace is exact Schrodinger-picture evolution under
    # H_sim = delta/2 (x^2+p^2) + B cos(2 pi x / Lambda + Phi). This bypasses
    # the trigonometric-gate/Trotter sequence and evolves the prepared motional
    # state directly under that Hamiltonian.
    prepared_state = apply_prep(initial_state(cutoff), prep, cutoff)
    motional_state = postselected_motional_vector(prepared_state, cutoff, "down")

    hamiltonian = hsim_hamiltonian(model, cutoff)
    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    coefficients = eigenvectors.conj().T @ motional_state
    x_in_energy_basis = eigenvectors.conj().T @ position_matrix(cutoff) @ eigenvectors

    times_s = times_ms * 1e-3
    values = []
    for time_s in times_s:
        evolved_coefficients = np.exp(-1.0j * energies * time_s) * coefficients
        values.append(
            float(
                np.vdot(
                    evolved_coefficients,
                    x_in_energy_basis @ evolved_coefficients,
                ).real
            )
        )
    return np.array(values, dtype=np.float64)


def compiled_gate_x_trace(
    prep: list[GateStatement],
    model: dict[str, float],
    max_time_ms: float,
    cutoff: int,
    readout_state: str = "down",
) -> tuple[np.ndarray, np.ndarray]:
    # Ideal semantic simulation of the generated Sandia xCD/Rz gate pattern.
    # For times beyond the parsed file, extend the same McGarry step rule
    # inferred from the Jaqal: alpha_k = alpha0 exp(i(phi0 + delta k dt)).
    max_time_s = max_time_ms * 1e-3
    dt = model["dt"]
    steps = round(max_time_s / dt)
    if not math.isclose(steps * dt, max_time_s, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("--hsim-max-time-ms must be an integer multiple of the inferred timestep")

    state = apply_prep(initial_state(cutoff), prep, cutoff)
    times_ms = [0.0]
    values = [postselected_x_expectation(state, cutoff, readout_state)]

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
        times_ms.append(elapsed * 1000.0)
        values.append(postselected_x_expectation(snapshot, cutoff, readout_state))

    return np.array(times_ms, dtype=np.float64), np.array(values, dtype=np.float64)


def load_program_model_and_snapshots(
    jaqal: Path,
    times_ms: list[float],
    cutoff: int,
    dt_us: float,
    use_schrodinger_frame: bool,
) -> tuple[dict[str, float], dict[float, np.ndarray]]:
    gates = gate_program(jaqal)
    prep, steps = split_program(gates)
    parameters = infer_step_parameters(steps, dt_us * 1e-6)
    model = infer_double_well(parameters)
    snapshots = run_snapshots(
        prep,
        steps,
        parameters,
        times_ms,
        cutoff,
        use_schrodinger_frame=use_schrodinger_frame,
    )
    return model, snapshots


def load_program_for_hsim_trace(
    jaqal: Path,
    dt_us: float,
) -> tuple[dict[str, float], list[GateStatement]]:
    gates = gate_program(jaqal)
    prep, steps = split_program(gates)
    parameters = infer_step_parameters(steps, dt_us * 1e-6)
    return infer_double_well(parameters), prep


def potential_curve(xs: np.ndarray, model: dict[str, float]) -> np.ndarray:
    curve = 0.5 * model["delta"] * (xs**2)
    curve += model["amplitude"] * np.cos(
        (2.0 * math.pi * xs / model["lambda_period"]) + model["varphi"]
    )
    return curve
