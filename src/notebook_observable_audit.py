from __future__ import annotations

"""Audit the saved notebook's probability-bin observable.

This intentionally simulates the notebook path, not the intended McGarry
observable: q0 does preparation/evolution, q1 does readout, both qubits are
measured, and the saved notebook computes (prob[0] - prob[2])/(prob[0]+prob[2]).
The goal is to identify which notebook/hardware convention changes can reproduce
the saved repeat traces in an ideal no-noise model.
"""

import argparse
import csv
import json
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

from simulator import REPO_ROOT, displacement_matrix, r_matrix, rz_matrix


DEFAULT_REPEAT_DIR = REPO_ROOT / "results" / "20260625_Data" / "003_1000 repeats"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "build" / "notebook_observable_audit"
DEFAULT_IM_BETAS = np.linspace(-0.4, 0.4, 5)
DEFAULT_RE_BETA = -0.2
DEFAULT_PHI = -0.8
DEFAULT_PREP_BETA = 0.89021208217480396
DEFAULT_CUTOFF = 32
DEFAULT_SELECTED_STEPS = [0, 13, 26, 49]


@dataclass(frozen=True)
class Candidate:
    label: str
    prep_scale: float
    evolution_scale: float
    readout_scale: float
    re_beta: float
    rotation_sign: int
    gate_order: str
    bit_order: str
    down_bit: int
    bins: tuple[int, int]
    predicted: np.ndarray
    expz_rmse: float
    slope_rmse: float
    t0_rmse: float


class NotebookObservableSimulator:
    def __init__(self, cutoff: int, phi: float):
        self.cutoff = cutoff
        self.phi = phi
        self._displacements: dict[tuple[float, float], np.ndarray] = {}
        self._hadamard = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128) / math.sqrt(2.0)

    def displacement(self, beta: complex) -> np.ndarray:
        key = (round(float(beta.real), 15), round(float(beta.imag), 15))
        if key not in self._displacements:
            self._displacements[key] = displacement_matrix(beta, self.cutoff)
        return self._displacements[key]

    def initial_state(self) -> np.ndarray:
        state = np.zeros((2, 2, self.cutoff), dtype=np.complex128)
        state[1, 1, 0] = 1.0
        return state

    def apply_zsdf(self, state: np.ndarray, qubit: int, beta: complex) -> np.ndarray:
        output = np.empty_like(state)
        displacement_plus = self.displacement(beta)
        displacement_minus = self.displacement(-beta)
        if qubit == 0:
            output[0] = np.einsum("nm,bm->bn", displacement_plus, state[0], optimize=True)
            output[1] = np.einsum("nm,bm->bn", displacement_minus, state[1], optimize=True)
            return output

        output[:, 0] = np.einsum("nm,am->an", displacement_plus, state[:, 0], optimize=True)
        output[:, 1] = np.einsum("nm,am->an", displacement_minus, state[:, 1], optimize=True)
        return output

    def apply_xsdf(self, state: np.ndarray, qubit: int, beta: complex) -> np.ndarray:
        displacement_plus = self.displacement(beta)
        displacement_minus = self.displacement(-beta)
        if qubit == 0:
            x_state = np.einsum("xa,abn->xbn", self._hadamard.conj().T, state, optimize=True)
            x_state[0] = np.einsum("nm,bm->bn", displacement_plus, x_state[0], optimize=True)
            x_state[1] = np.einsum("nm,bm->bn", displacement_minus, x_state[1], optimize=True)
            return np.einsum("ax,xbn->abn", self._hadamard, x_state, optimize=True)

        x_state = np.einsum("xb,abn->axn", self._hadamard.conj().T, state, optimize=True)
        x_state[:, 0] = np.einsum("nm,am->an", displacement_plus, x_state[:, 0], optimize=True)
        x_state[:, 1] = np.einsum("nm,am->an", displacement_minus, x_state[:, 1], optimize=True)
        return np.einsum("bx,axn->abn", self._hadamard, x_state, optimize=True)

    def apply_rz_q0(self, state: np.ndarray) -> np.ndarray:
        rotation = rz_matrix(self.phi)
        return np.einsum("ax,xbn->abn", rotation, state, optimize=True)

    def apply_readout_rotation_q1(self, state: np.ndarray, rotation_sign: int) -> np.ndarray:
        rotation = r_matrix(0.0, rotation_sign * math.pi / 2.0)
        return np.einsum("bx,axn->abn", rotation, state, optimize=True)

    def probability_by_int(self, state: np.ndarray, bit_order: str, down_bit: int) -> np.ndarray:
        probabilities_by_qubit = np.sum(np.abs(state) ** 2, axis=2)
        probabilities = np.zeros(4, dtype=np.float64)
        for q0_internal in (0, 1):
            for q1_internal in (0, 1):
                q0_bit = down_bit if q0_internal == 1 else 1 - down_bit
                q1_bit = down_bit if q1_internal == 1 else 1 - down_bit
                if bit_order == "q0_lsb":
                    index = q0_bit + 2 * q1_bit
                elif bit_order == "q0_msb":
                    index = 2 * q0_bit + q1_bit
                else:
                    raise ValueError(f"unsupported bit order {bit_order!r}")
                probabilities[index] = probabilities_by_qubit[q0_internal, q1_internal]
        return probabilities

    @staticmethod
    def notebook_exp_z(probabilities: np.ndarray, bins: tuple[int, int]) -> float:
        denominator = probabilities[bins[0]] + probabilities[bins[1]]
        if denominator <= 1e-15:
            return 0.0
        return float((probabilities[bins[0]] - probabilities[bins[1]]) / denominator)

    def prepared_state(self, prep_scale: float) -> np.ndarray:
        state = self.initial_state()
        return self.apply_zsdf(state, 0, prep_scale * DEFAULT_PREP_BETA)

    def evolution_states(
        self,
        angles: np.ndarray,
        prep_scale: float,
        evolution_scale: float,
        gate_order: str,
    ) -> list[np.ndarray]:
        state = self.prepared_state(prep_scale)
        states = [state]
        for re_angle, im_angle in angles:
            alpha = evolution_scale * complex(float(re_angle), float(im_angle))
            if gate_order == "saved":
                sequence: list[complex | str] = [alpha, "rz", -alpha, -alpha, "rz", alpha]
            elif gate_order == "corrected":
                sequence = [-alpha, "rz", alpha, alpha, "rz", -alpha]
            else:
                raise ValueError(f"unsupported gate order {gate_order!r}")
            for item in sequence:
                if item == "rz":
                    state = self.apply_rz_q0(state)
                else:
                    state = self.apply_xsdf(state, 0, item)
            states.append(state)
        return states

    def predict(
        self,
        angles: np.ndarray,
        im_betas: np.ndarray,
        prep_scale: float,
        evolution_scale: float,
        readout_scale: float,
        re_beta: float,
        rotation_sign: int,
        gate_order: str,
        bit_order: str,
        down_bit: int,
        bins: tuple[int, int],
    ) -> np.ndarray:
        states = self.evolution_states(angles, prep_scale, evolution_scale, gate_order)
        predicted = np.empty((im_betas.size, len(states)), dtype=np.float64)
        for beta_index, im_beta in enumerate(im_betas):
            pulse_beta = readout_scale * complex(re_beta, float(im_beta))
            for step, state in enumerate(states):
                readout_state = self.apply_readout_rotation_q1(state, rotation_sign)
                readout_state = self.apply_xsdf(readout_state, 1, pulse_beta)
                probabilities = self.probability_by_int(readout_state, bit_order, down_bit)
                predicted[beta_index, step] = self.notebook_exp_z(probabilities, bins)
        return predicted


def parse_saved_angles(output_notebook: Path) -> np.ndarray:
    notebook = json.loads(output_notebook.read_text())
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if "print(angle)" not in source:
            continue
        for output in cell.get("outputs", []):
            if output.get("output_type") != "stream":
                continue
            text = "".join(output.get("text", []))
            match = re.search(r"\[\[.*?\]\]", text, flags=re.DOTALL)
            if match is None:
                continue
            values = np.fromstring(match.group(0).replace("[", " ").replace("]", " "), sep=" ")
            if values.size % 2 != 0:
                raise ValueError(f"could not parse angle pairs from {output_notebook}")
            return values.reshape((-1, 2))
    raise ValueError(f"could not find printed angle array in {output_notebook}")


def fit_slope_trace(values: np.ndarray, im_betas: np.ndarray) -> np.ndarray:
    design = np.column_stack([im_betas, np.ones_like(im_betas)])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return -coefficients[0] / math.sqrt(2.0)


def metrics(predicted: np.ndarray, data: np.ndarray, im_betas: np.ndarray) -> tuple[float, float, float]:
    expz_rmse = float(np.sqrt(np.mean((predicted - data) ** 2)))
    slope_rmse = float(np.sqrt(np.mean((fit_slope_trace(predicted, im_betas) - fit_slope_trace(data, im_betas)) ** 2)))
    t0_rmse = float(np.sqrt(np.mean((predicted[:, 0] - data[:, 0]) ** 2)))
    return expz_rmse, slope_rmse, t0_rmse


def candidate_label(candidate: Candidate) -> str:
    return (
        f"prep={candidate.prep_scale:g}, evol={candidate.evolution_scale:g}, "
        f"readout={candidate.readout_scale:g}, order={candidate.gate_order}"
    )


def scan_candidates(
    simulator: NotebookObservableSimulator,
    angles: np.ndarray,
    data: np.ndarray,
    im_betas: np.ndarray,
) -> list[Candidate]:
    candidates = []
    for prep_scale in [1.0, 0.75, 0.5, 0.45, 0.4, 0.25]:
        for evolution_scale in [1.0, 0.5, 0.25, 0.0]:
            for readout_scale in [1.0, 0.5, 0.25, -1.0, -0.5, -0.25]:
                for gate_order in ["saved", "corrected"]:
                    predicted = simulator.predict(
                        angles,
                        im_betas,
                        prep_scale=prep_scale,
                        evolution_scale=evolution_scale,
                        readout_scale=readout_scale,
                        re_beta=DEFAULT_RE_BETA,
                        rotation_sign=1,
                        gate_order=gate_order,
                        bit_order="q0_lsb",
                        down_bit=0,
                        bins=(0, 2),
                    )
                    expz_rmse, slope_rmse, t0_rmse = metrics(predicted, data, im_betas)
                    candidates.append(
                        Candidate(
                            label="grid",
                            prep_scale=prep_scale,
                            evolution_scale=evolution_scale,
                            readout_scale=readout_scale,
                            re_beta=DEFAULT_RE_BETA,
                            rotation_sign=1,
                            gate_order=gate_order,
                            bit_order="q0_lsb",
                            down_bit=0,
                            bins=(0, 2),
                            predicted=predicted,
                            expz_rmse=expz_rmse,
                            slope_rmse=slope_rmse,
                            t0_rmse=t0_rmse,
                        )
                    )
    return sorted(candidates, key=lambda candidate: candidate.expz_rmse)


def important_candidates(
    simulator: NotebookObservableSimulator,
    angles: np.ndarray,
    data: np.ndarray,
    im_betas: np.ndarray,
) -> list[Candidate]:
    configs = [
        ("as written direct xSDF arg", 1.0, 1.0, 1.0, "saved"),
        ("McGarry half-readout pulse", 1.0, 1.0, 0.5, "saved"),
        ("sign flipped direct xSDF arg", 1.0, 1.0, -1.0, "saved"),
        ("sign flipped McGarry half-readout", 1.0, 1.0, -0.5, "saved"),
        ("sign flipped + half prep", 0.5, 1.0, -0.5, "saved"),
        ("sign flipped + extra half readout", 1.0, 1.0, -0.25, "saved"),
        ("weak evolution best-family", 0.25, 0.5, -0.5, "saved"),
        ("no evolution best-family", 0.25, 0.0, -0.5, "saved"),
    ]
    candidates = []
    for label, prep_scale, evolution_scale, readout_scale, gate_order in configs:
        predicted = simulator.predict(
            angles,
            im_betas,
            prep_scale=prep_scale,
            evolution_scale=evolution_scale,
            readout_scale=readout_scale,
            re_beta=DEFAULT_RE_BETA,
            rotation_sign=1,
            gate_order=gate_order,
            bit_order="q0_lsb",
            down_bit=0,
            bins=(0, 2),
        )
        expz_rmse, slope_rmse, t0_rmse = metrics(predicted, data, im_betas)
        candidates.append(
            Candidate(
                label=label,
                prep_scale=prep_scale,
                evolution_scale=evolution_scale,
                readout_scale=readout_scale,
                re_beta=DEFAULT_RE_BETA,
                rotation_sign=1,
                gate_order=gate_order,
                bit_order="q0_lsb",
                down_bit=0,
                bins=(0, 2),
                predicted=predicted,
                expz_rmse=expz_rmse,
                slope_rmse=slope_rmse,
                t0_rmse=t0_rmse,
            )
        )
    return candidates


def write_summary(candidates: list[Candidate], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "label",
                "expz_rmse",
                "slope_rmse",
                "t0_rmse",
                "prep_scale",
                "evolution_scale",
                "readout_scale",
                "re_beta",
                "rotation_sign",
                "gate_order",
                "bit_order",
                "down_bit",
                "bins",
            ],
        )
        writer.writeheader()
        for rank, candidate in enumerate(candidates, 1):
            writer.writerow(
                {
                    "rank": rank,
                    "label": candidate.label,
                    "expz_rmse": f"{candidate.expz_rmse:.17g}",
                    "slope_rmse": f"{candidate.slope_rmse:.17g}",
                    "t0_rmse": f"{candidate.t0_rmse:.17g}",
                    "prep_scale": f"{candidate.prep_scale:.17g}",
                    "evolution_scale": f"{candidate.evolution_scale:.17g}",
                    "readout_scale": f"{candidate.readout_scale:.17g}",
                    "re_beta": f"{candidate.re_beta:.17g}",
                    "rotation_sign": candidate.rotation_sign,
                    "gate_order": candidate.gate_order,
                    "bit_order": candidate.bit_order,
                    "down_bit": candidate.down_bit,
                    "bins": f"{candidate.bins[0]} {candidate.bins[1]}",
                }
            )


def plot_slope_overlay(
    data: np.ndarray,
    im_betas: np.ndarray,
    candidates: list[Candidate],
    selected_steps: list[int],
    output: Path,
) -> None:
    times_ms = np.arange(data.shape[1], dtype=np.float64) * (0.8 / 5.09628e3 * 1.0e3)
    data_slope = fit_slope_trace(data, im_betas)

    figure, axis = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    axis.plot(times_ms, data_slope, color="#2d7f5e", linewidth=1.7, marker="o", markersize=3.0, label="saved notebook expZ slope")

    colors = ["#1f5aa6", "#b3433f", "#7b4a9e", "#d17827", "#6a6a6a"]
    linestyles = ["-", "--", ":", "-.", (0, (5, 2, 1, 2))]
    for index, candidate in enumerate(candidates):
        label = candidate.label if candidate.label != "grid" else candidate_label(candidate)
        axis.plot(
            times_ms,
            fit_slope_trace(candidate.predicted, im_betas),
            color=colors[index % len(colors)],
            linestyle=linestyles[index % len(linestyles)],
            linewidth=1.4,
            label=label,
        )

    for step in selected_steps:
        axis.axvline(times_ms[step], color="black", linewidth=0.6, alpha=0.14)

    axis.axhline(0.0, color="black", linewidth=0.8, alpha=0.28)
    axis.set_xlabel("time (ms)")
    axis.set_ylabel(r"notebook-extracted $\langle x\rangle$ proxy")
    axis.set_title("Saved notebook observable versus ideal notebook-observable hypotheses")
    axis.grid(alpha=0.16)
    axis.legend(loc="best", fontsize=8, frameon=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Simulate and audit the saved notebook probability-bin observable.")
    parser.add_argument("--repeat-dir", type=Path, default=DEFAULT_REPEAT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cutoff", type=int, default=DEFAULT_CUTOFF)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--plot-top", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_notebook = args.repeat_dir / "output.ipynb"
    data_path = args.repeat_dir / "expZ_imMeas.npy"
    data = np.load(data_path)
    angles = parse_saved_angles(output_notebook)
    im_betas = DEFAULT_IM_BETAS

    simulator = NotebookObservableSimulator(args.cutoff, DEFAULT_PHI)
    scanned = scan_candidates(simulator, angles, data, im_betas)
    important = important_candidates(simulator, angles, data, im_betas)

    summary_output = args.output_dir / f"{args.repeat_dir.name.replace(' ', '_')}_candidate_scan.csv"
    important_output = args.output_dir / f"{args.repeat_dir.name.replace(' ', '_')}_important_candidates.csv"
    plot_output = args.output_dir / f"{args.repeat_dir.name.replace(' ', '_')}_slope_overlay.png"

    write_summary(scanned, summary_output)
    write_summary(important, important_output)
    plot_candidates = important + scanned[: min(args.plot_top, len(scanned))]
    plot_slope_overlay(data, im_betas, plot_candidates, DEFAULT_SELECTED_STEPS, plot_output)

    print(summary_output)
    print(important_output)
    print(plot_output)
    print("top candidates:")
    for candidate in scanned[: args.top]:
        print(
            f"{candidate_label(candidate)} "
            f"expZ_rmse={candidate.expz_rmse:.4f} "
            f"slope_rmse={candidate.slope_rmse:.4f} "
            f"t0_rmse={candidate.t0_rmse:.4f}"
        )


if __name__ == "__main__":
    main()
