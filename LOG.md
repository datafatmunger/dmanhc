# DMANH Compiler Log

## 2026-06-10

- Checked McGarry/QSCOUT characteristic-function readout conventions against the local reference material.
  - McGarry readout measures `chi(beta)=<D(beta)>` by applying `D(sigma_x beta/2)` followed by Z-basis qubit measurement.
  - No prepended pulse gives `Re[chi(beta)]`; prepending `Rx(pi/2)` gives `Im[chi(beta)]`.
  - QSCOUT `R q 0 pi/2` is an X-axis equatorial rotation, not `Rz`.

- Fixed the ideal readout check in `src/measure.py`.
  - Previous check used an `Rz(pi/2)` matrix for the imaginary selector.
  - Current check uses the same generic X-axis `R(0, pi/2)` pulse emitted by the compiler.
  - CSV diagnostic column names were changed from `rz_pi_2` to `rx_pi_2`.

- Clarified generated evolution-step comments in `src/compiler.py`.
  - Previous comments looked like a measured `t=0` state: `// step 0: t = 0 ms`.
  - Current comments label the first block as an evolution interval: `evolution step 1: k = 0, interval 0 -> Delta t`.
  - No separate zero-evolution measurement circuit has been added yet.

- Updated DMANH+ initialization parameters.
  - Phil's Section 3 computes the well minimum from `V(x)=delta x^2/2 + B cos(2 pi x / Lambda)`.
  - With `delta=1.29817e3`, `B=5.09628e3`, and `alpha0=pi/(12 sqrt(2))`, the minima are near `x=+-1.25895`.
  - The corresponding preparation value is `x_min/sqrt(2)=0.8902`.
  - `Makefile` now uses `DMANH_X_MIN ?= 1.25895`.

- Switched the hardware-facing alpha phase convention.
  - Previous local convention used `ALPHA_PHASE_OFFSET=-pi/2` so the abstract mathematical displacement simulator plotted position-space dynamics correctly.
  - Phil's Sandia/DMANH convention treats the Jaqal `xCD` argument as the requested McGarry `alpha0 exp(i zeta)`, with `zeta=k delta Delta t`.
  - `ALPHA_PHASE_OFFSET` is now `0`, with the old `-pi/2` value left in a comment.

- Inlined the Sandia xCD convention at the relevant call sites.
  - Previous local mathematical convention is kept commented where the mapping is applied.
  - `src/simulator.py` maps a Sandia xCD argument `s` to the mathematical displacement parameter `-i s`.
  - `src/compiler.py` maps McGarry readout `beta` to emitted Sandia argument `s=i beta/2`.
  - `src/measure.py` builds the ideal readout unitary directly as mathematical `D(sigma_x beta/2)`.

- Updated readout mapping for the new Sandia/DMANH convention.
  - Since `xCD(s)` now emulates mathematical `D(-i s)`, McGarry readout `D(sigma_x beta/2)` requires Sandia argument `s=i beta/2`.
  - For the default 2PFD coordinate `beta=i 0.4`, generated Jaqal should use `reBeta=-0.2`, `imBeta=0`.
  - Added `MEASUREMENT_DERIVATION.md` with the hand calculation for `Re[chi]`,
    `Im[chi]`, and the `s=i beta/2` compiler mapping.

- Regenerated `build/dmanh.jaqal`.
  - Preparation now emits `zCD ... 0.89021208217480396 0`.
  - Evolution uses `alpha phase offset = 0`.
  - First evolution `xCD` is real, matching Phil's expectation for `k=0`.

- Expanded the `make dmanh` plotting outputs.
  - `src/plots.py` now emits both `build/dmanh.png` and `build/dmanh_hsim.png` for the DMANH+ target.
  - The DMANH+ `H_sim` trace uses `--hsim-max-time-ms 7.6918850612603702`, the exact `49 Delta t` endpoint required by the timestep-grid check.
  - `src/measure.py` now emits `build/dmanh_measurement_panels.png` and `build/dmanh_chi_slice_panels.png` for the same DMANH+ times.

- Replaced the attempted `dmanh-4ms`/`dmanh-compact` Phil-facing variants with `make dmanh-vartheta-1p6`.
  - This keeps the fitted DMANH+ well and dynamics fixed: `B=5.09628e3 rad/s`, `delta=1.29817e3 rad/s`, `alpha0=0.18512`, and `x_min=1.25895`.
  - It increases only `vartheta` from `0.8` to `1.6`, so `Delta t = vartheta / B = 313.954492 us`.
  - Because the original `49`-step endpoint at `7.691885 ms` is halfway between coarse steps, the target uses `25` steps and ends at `7.848862 ms`; the midpoint `4.081408 ms` is exactly `13` coarse steps.
  - It writes `build/dmanh_vartheta_1p6.jaqal`, `build/dmanh_vartheta_1p6.png`, `build/dmanh_vartheta_1p6_hsim.png`, `build/dmanh_vartheta_1p6_measurement_panels.png`, and `build/dmanh_vartheta_1p6_chi_slice_panels.png`.

- Added Phil's FFT frequency-error sweep as `make dmanh-frequency-sweep`.
  - `src/frequency_sweep.py` compares exact `H_sim` propagation with compiled xCD/Rz propagation while sweeping `vartheta=0.1..3.0`.
  - The sweep keeps `B`, `delta`, `alpha0`, and `x_min` fixed and changes only `Delta t = vartheta / B`.
  - It writes `build/dmanh_frequency_sweep.csv`, `build/dmanh_frequency_shift_vs_vartheta.png`, `build/dmanh_frequency_spectra.png`, and `build/dmanh_frequency_trace_examples.png`.

## 2026-06-19

- Agent/model: Codex, GPT-5.

- Added a read-only Sandia-notebook comparison tool as `src/notebook_readout_compare.py`.
  - The script parses `notebooks/output.ipynb` for the hard-coded angle list, `im_beta_list`, `reBeta`, and the printed `expZ_imMeas` matrix without modifying the notebooks.
  - It simulates the same prefix circuits with the existing `src/measure.py` and `src/simulator.py` conventions.
  - It also generates a separate prepped left-well Eq. `xapprox` comparison using `x_min=1.25895`.

- Added `make dmanh-readout-compare`.
  - This writes `build/notebook_readout_compare_output_vs_theory.png`, `build/notebook_readout_compare_eq35.png`, `build/notebook_readout_compare.csv`, and `build/notebook_readout_compare_eq35_prefix.jaqal`.
  - Verified with `.venv/bin/python3 -m py_compile src/notebook_readout_compare.py` and `make dmanh-readout-compare`.

- Convention finding from the Sandia output notebook:
  - The notebook sweeps the Jaqal variable `imBeta` while keeping `reBeta=-0.2`.
  - Under the current local convention, these are Sandia `xCD` arguments `s`, with McGarry `beta=-2 i s`.
  - Therefore the notebook's sweep maps to `beta = 2*imBeta + 0.4 i`: it varies `Re[beta]` at fixed `Im[beta]=0.4`.
  - A Fig. 3(f) / Eq. `xapprox` scan along `beta=i y` should instead keep `imBeta=0` and sweep `reBeta=-y/2`; the 2PFD point with `h=0.4` is `reBeta=-0.2`, `imBeta=0`.

- Current comparison result:
  - For the notebook's no-prep prefix circuit, ideal `measure.py` theory gives `<x>=0` and `Im[chi(i y)]=0` up to numerical roundoff.
  - The printed Sandia output matrix has the center `imBeta=0` trace near `1` for all three subcircuits, so it is not reproduced by the ideal no-prep symmetric circuit.
  - With a left-well `zCD` preparation, the generated Eq. `xapprox` comparator gives direct `<x>` around `-1.257`, `-1.201`, and `-1.096` after one, two, and three prefix blocks; the `h=0.4` 2PFD estimates are about `-1.064`, `-1.028`, and `-0.959`.

Open questions:

- Whether to generate separate hardware programs for `0`, `26`, and `49` evolution steps before readout.
- Whether the `zCD` simulator path should also be routed through a named convention helper, or remain direct because it is being used as a preparation shortcut.
- Whether to clean the appended scratch notes currently present at the end of `README.md`.
