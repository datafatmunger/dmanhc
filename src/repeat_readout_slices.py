from __future__ import annotations

"""Plot saved repeat-count readout arrays against characteristic-function theory.

This script deliberately reads saved notebook outputs only. It does not execute
the notebooks in notebooks/001_50 repeats, notebooks/002_500 repeats, or
notebooks/003_1000 repeats.
"""

import argparse
import csv
import math
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cvdv-matplotlib-cache"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from measure import characteristic_function, motional_density_matrix, position_expectation_from_rho
from simulator import FOCK_CUTOFF, apply_rz, apply_xcd, apply_zcd, initial_state


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "build" / "repeat_readout_slices"
DEFAULT_ANGLES = REPO_ROOT / "build" / "notebook" / "angles.npy"
DEFAULT_REPEAT_DIRS = [
    REPO_ROOT / "notebooks" / "001_50 repeats",
    REPO_ROOT / "notebooks" / "002_500 repeats",
    REPO_ROOT / "notebooks" / "003_1000 repeats",
]
DEFAULT_IM_BETAS = np.linspace(-0.4, 0.4, 5)
DEFAULT_SANDIA_RE_BETA = -0.2
DEFAULT_PREP_X_MIN = 1.25895
DEFAULT_DT_MS = 0.15697724614817082
DEFAULT_SELECTED_STEPS = [0, 13, 26, 49]


def repeat_count(path: Path) -> int:
    head = path.name.split("_", maxsplit=2)[1]
    return int(head.split()[0])


def load_angles(path: Path) -> np.ndarray:
    angles = np.load(path)
    if angles.shape[0] == 3:
        angles = angles.T
    if angles.ndim != 2 or angles.shape[1] != 3:
        raise ValueError(f"expected angles as (3, steps) or (steps, 3), got {angles.shape}")
    return angles


def notebook_sandia_to_mcgarry_beta(sandia_beta: complex) -> complex:
    # Current local convention: xCD(s) realizes D(-i s), while the McGarry
    # readout applies D(sigma_x beta / 2). Therefore s = i beta / 2.
    return -2.0j * sandia_beta


def sandia_re_for_fig3f_y(y: float) -> float:
    # beta = i*y and s = i*beta/2 imply s = -y/2, a real Sandia xCD argument.
    return -0.5 * y


def apply_notebook_block(state: np.ndarray, angle: np.ndarray, cutoff: int) -> np.ndarray:
    alpha = complex(float(angle[0]), float(angle[1]))
    theta = float(angle[2])
    state = apply_xcd(state, alpha, cutoff)
    state = apply_rz(state, theta, cutoff)
    state = apply_xcd(state, -alpha, cutoff)
    state = apply_xcd(state, -alpha, cutoff)
    state = apply_rz(state, theta, cutoff)
    return apply_xcd(state, alpha, cutoff)


def prefix_rhos(
    angles: np.ndarray,
    cutoff: int,
    prep_x_min: float,
    readout_state: str,
) -> tuple[list[np.ndarray], np.ndarray]:
    state = initial_state(cutoff)
    if not math.isclose(prep_x_min, 0.0, rel_tol=0.0, abs_tol=1e-15):
        state = apply_zcd(state, prep_x_min / math.sqrt(2.0), cutoff)

    rhos = []
    x_values = []
    rho, _ = motional_density_matrix(state, cutoff, readout_state)
    rhos.append(rho)
    x_values.append(position_expectation_from_rho(rho, cutoff))

    for angle in angles:
        state = apply_notebook_block(state, angle, cutoff)
        rho, _ = motional_density_matrix(state, cutoff, readout_state)
        rhos.append(rho)
        x_values.append(position_expectation_from_rho(rho, cutoff))

    return rhos, np.array(x_values)


def theory_for_notebook_sweep(
    rhos: list[np.ndarray],
    sandia_re_beta: float,
    sandia_im_betas: np.ndarray,
    cutoff: int,
) -> tuple[np.ndarray, np.ndarray]:
    beta_values = np.array(
        [
            notebook_sandia_to_mcgarry_beta(complex(sandia_re_beta, sandia_im_beta))
            for sandia_im_beta in sandia_im_betas
        ],
        dtype=np.complex128,
    )
    values = np.empty((sandia_im_betas.size, len(rhos)), dtype=np.float64)
    for row, beta in enumerate(beta_values):
        for col, rho in enumerate(rhos):
            values[row, col] = characteristic_function(beta, rho, cutoff).imag
    return beta_values, values


def theory_for_fig3f(
    rhos: list[np.ndarray],
    y_values: np.ndarray,
    cutoff: int,
) -> np.ndarray:
    values = np.empty((y_values.size, len(rhos)), dtype=np.float64)
    for row, y in enumerate(y_values):
        beta = 1.0j * float(y)
        for col, rho in enumerate(rhos):
            values[row, col] = characteristic_function(beta, rho, cutoff).imag
    return values


def shot_stderr(values: np.ndarray, repeats: int) -> np.ndarray:
    # Approximate because the saved expZ array does not preserve raw counts or
    # the postselection denominator used by exp_z in the notebooks.
    clipped = np.clip(values, -1.0, 1.0)
    return np.sqrt(np.maximum(0.0, 1.0 - clipped**2) / repeats)


def plot_repeat(
    repeat_dir: Path,
    expz: np.ndarray,
    repeats: int,
    selected_steps: list[int],
    sandia_im_betas: np.ndarray,
    notebook_beta_values: np.ndarray,
    notebook_theory: np.ndarray,
    fig3f_y: np.ndarray,
    fig3f_theory: np.ndarray,
    x_values: np.ndarray,
    output: Path,
) -> None:
    fig, axes = plt.subplots(
        len(selected_steps),
        2,
        figsize=(11.5, 2.7 * len(selected_steps)),
        sharex="col",
        sharey=True,
        constrained_layout=True,
    )
    if len(selected_steps) == 1:
        axes = np.array([axes])

    for row, step in enumerate(selected_steps):
        left = axes[row, 0]
        right = axes[row, 1]
        time_ms = step * DEFAULT_DT_MS
        stderr = shot_stderr(expz[:, step], repeats)

        left.axhline(0.0, color="black", linewidth=0.8, alpha=0.22)
        left.errorbar(
            sandia_im_betas,
            expz[:, step],
            yerr=stderr,
            fmt="o",
            color="#1f5aa6",
            ecolor="#1f5aa6",
            elinewidth=0.9,
            capsize=2.5,
            label="saved expZ",
        )
        left.plot(
            sandia_im_betas,
            notebook_theory[:, step],
            color="#b3433f",
            linewidth=1.7,
            label="ideal, same Sandia sweep",
        )
        left.set_title(f"step {step}, t={time_ms:.3f} ms")
        left.set_ylabel("Im chi / expZ")
        left.grid(alpha=0.16)

        right.axhline(0.0, color="black", linewidth=0.8, alpha=0.22)
        right.plot(
            fig3f_y,
            fig3f_theory[:, step],
            color="#2d7f5e",
            linewidth=1.9,
            label="ideal Fig. 3(f) scan",
        )
        right.scatter(
            [0.4],
            [np.interp(0.4, fig3f_y, fig3f_theory[:, step])],
            color="#2d7f5e",
            s=26,
            zorder=3,
            label="h=0.4",
        )
        slope = math.sqrt(2.0) * x_values[step]
        local_line = slope * fig3f_y
        right.plot(
            fig3f_y,
            local_line,
            color="#626262",
            linewidth=1.0,
            linestyle=":",
            label="origin slope",
        )
        right.set_title(f"proper beta=i*y; <x>={x_values[step]:+.3f}")
        right.grid(alpha=0.16)

    axes[-1, 0].set_xlabel("notebook Sandia imBeta override")
    axes[-1, 1].set_xlabel("McGarry y = Im[beta]")
    axes[0, 0].legend(loc="best", fontsize=8, frameon=True)
    axes[0, 1].legend(loc="best", fontsize=8, frameon=True)
    mapped = ", ".join(
        f"{beta.real:+.1f}{beta.imag:+.1f}i" for beta in notebook_beta_values
    )
    fig.suptitle(
        (
            f"{repeat_dir.name}: saved sweep versus Fig. 3(f)-style theory\n"
            f"left maps s=-0.2+i*imBeta to McGarry beta values [{mapped}]"
        ),
        fontsize=11,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_step26_overlay(
    repeat_data: list[tuple[Path, int, np.ndarray]],
    sandia_im_betas: np.ndarray,
    notebook_theory: np.ndarray,
    fig3f_y: np.ndarray,
    fig3f_theory: np.ndarray,
    output: Path,
    step: int = 26,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True, sharey=True)

    axes[0].axhline(0.0, color="black", linewidth=0.8, alpha=0.22)
    for repeat_dir, repeats, expz in repeat_data:
        axes[0].errorbar(
            sandia_im_betas,
            expz[:, step],
            yerr=shot_stderr(expz[:, step], repeats),
            fmt="o-",
            linewidth=1.2,
            capsize=2.2,
            label=f"{repeats} repeats",
        )
    axes[0].plot(
        sandia_im_betas,
        notebook_theory[:, step],
        color="black",
        linewidth=1.8,
        linestyle="--",
        label="ideal same Sandia sweep",
    )
    axes[0].set_title(f"saved arrays at step {step}")
    axes[0].set_xlabel("notebook Sandia imBeta override")
    axes[0].set_ylabel("Im chi / expZ")
    axes[0].grid(alpha=0.16)
    axes[0].legend(loc="best", fontsize=8, frameon=True)

    axes[1].axhline(0.0, color="black", linewidth=0.8, alpha=0.22)
    axes[1].plot(
        fig3f_y,
        fig3f_theory[:, step],
        color="#2d7f5e",
        linewidth=2.0,
        label="ideal Fig. 3(f) scan",
    )
    axes[1].scatter(
        [0.4],
        [np.interp(0.4, fig3f_y, fig3f_theory[:, step])],
        color="#2d7f5e",
        s=28,
        label="h=0.4",
        zorder=3,
    )
    axes[1].set_title("proper McGarry coordinate")
    axes[1].set_xlabel("McGarry y = Im[beta]")
    axes[1].grid(alpha=0.16)
    axes[1].legend(loc="best", fontsize=8, frameon=True)

    fig.suptitle(
        f"Step {step} (t={step * DEFAULT_DT_MS:.3f} ms): saved sweep is not Fig. 3(f)'s beta=i*y scan",
        fontsize=11,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_summary(
    output: Path,
    repeat_data: list[tuple[Path, int, np.ndarray]],
    selected_steps: list[int],
    sandia_im_betas: np.ndarray,
    notebook_beta_values: np.ndarray,
    notebook_theory: np.ndarray,
    fig3f_y: np.ndarray,
    fig3f_theory: np.ndarray,
    x_values: np.ndarray,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        fieldnames = [
            "repeat_dir",
            "repeats",
            "step",
            "time_ms",
            "sandia_reBeta",
            "sandia_imBeta",
            "mapped_mcgarry_beta_re",
            "mapped_mcgarry_beta_im",
            "saved_expZ",
            "theory_same_sandia_sweep",
            "theory_fig3f_at_same_y",
            "x_direct",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for repeat_dir, repeats, expz in repeat_data:
            for step in selected_steps:
                for row, sandia_im_beta in enumerate(sandia_im_betas):
                    beta = notebook_beta_values[row]
                    fig3f_same_y = np.interp(
                        sandia_im_beta,
                        fig3f_y,
                        fig3f_theory[:, step],
                    )
                    writer.writerow(
                        {
                            "repeat_dir": repeat_dir.name,
                            "repeats": repeats,
                            "step": step,
                            "time_ms": f"{step * DEFAULT_DT_MS:.12g}",
                            "sandia_reBeta": f"{DEFAULT_SANDIA_RE_BETA:.17g}",
                            "sandia_imBeta": f"{sandia_im_beta:.17g}",
                            "mapped_mcgarry_beta_re": f"{beta.real:.17g}",
                            "mapped_mcgarry_beta_im": f"{beta.imag:.17g}",
                            "saved_expZ": f"{expz[row, step]:.17g}",
                            "theory_same_sandia_sweep": f"{notebook_theory[row, step]:.17g}",
                            "theory_fig3f_at_same_y": f"{fig3f_same_y:.17g}",
                            "x_direct": f"{x_values[step]:.17g}",
                        }
                    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot saved repeat-count expZ_imMeas arrays without running notebooks."
    )
    parser.add_argument("--angles", type=Path, default=DEFAULT_ANGLES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cutoff", type=int, default=FOCK_CUTOFF)
    parser.add_argument("--prep-x-min", type=float, default=DEFAULT_PREP_X_MIN)
    parser.add_argument("--readout-state", choices=["up", "down", "trace"], default="down")
    parser.add_argument("--fig3f-max", type=float, default=0.8)
    parser.add_argument("--fig3f-points", type=int, default=161)
    parser.add_argument("--selected-steps", type=int, nargs="+", default=DEFAULT_SELECTED_STEPS)
    parser.add_argument("--repeat-dir", type=Path, action="append", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    angles = load_angles(args.angles)
    sandia_im_betas = DEFAULT_IM_BETAS.copy()
    fig3f_y = np.linspace(-args.fig3f_max, args.fig3f_max, args.fig3f_points)

    rhos, x_values = prefix_rhos(
        angles,
        args.cutoff,
        args.prep_x_min,
        args.readout_state,
    )
    notebook_beta_values, notebook_theory = theory_for_notebook_sweep(
        rhos,
        DEFAULT_SANDIA_RE_BETA,
        sandia_im_betas,
        args.cutoff,
    )
    fig3f_theory = theory_for_fig3f(rhos, fig3f_y, args.cutoff)

    repeat_dirs = args.repeat_dir if args.repeat_dir else DEFAULT_REPEAT_DIRS
    repeat_data = []
    for repeat_dir in repeat_dirs:
        expz = np.load(repeat_dir / "expZ_imMeas.npy")
        if expz.shape != notebook_theory.shape:
            raise ValueError(
                f"{repeat_dir / 'expZ_imMeas.npy'} has shape {expz.shape}; "
                f"expected {notebook_theory.shape}"
            )
        repeats = repeat_count(repeat_dir)
        repeat_data.append((repeat_dir, repeats, expz))
        output = args.output_dir / f"{repeat_dir.name.split()[0]}_readout_slices.png"
        plot_repeat(
            repeat_dir,
            expz,
            repeats,
            args.selected_steps,
            sandia_im_betas,
            notebook_beta_values,
            notebook_theory,
            fig3f_y,
            fig3f_theory,
            x_values,
            output,
        )
        print(output)

    overlay_output = args.output_dir / "step26_overlay.png"
    plot_step26_overlay(
        repeat_data,
        sandia_im_betas,
        notebook_theory,
        fig3f_y,
        fig3f_theory,
        overlay_output,
    )
    print(overlay_output)

    csv_output = args.output_dir / "summary.csv"
    write_summary(
        csv_output,
        repeat_data,
        args.selected_steps,
        sandia_im_betas,
        notebook_beta_values,
        notebook_theory,
        fig3f_y,
        fig3f_theory,
        x_values,
    )
    print(csv_output)


if __name__ == "__main__":
    main()
