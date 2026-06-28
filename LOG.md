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

- Removed the obsolete standalone `dmanh-vartheta-1p6` experiment target and
  TOML. Coarse-timestep checks now live in the frequency sweep rather than as
  a separate named experiment.

- Added Phil's FFT frequency-error sweep as `make dmanh-frequency-sweep`.
  - `src/frequency_sweep.py` compares exact `H_sim` propagation with compiled xCD/Rz propagation while sweeping `vartheta=0.1..3.0`.
  - The sweep keeps `B`, `delta`, `alpha0`, and `x_min` fixed and changes only `Delta t = vartheta / B`.
  - It writes `build/dmanh_frequency_sweep.csv`, `build/dmanh_frequency_spectra.png`, and `build/dmanh_frequency_trace_examples.png`.

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

## 2026-06-27

- Agent/model: Codex, GPT-5.

- Added a read-only slope-extraction analysis as `src/readout_x_trace.py` and
  `make readout-x-trace`.
  - The script reads the saved `expZ_imMeas.npy` arrays without modifying
    `results/20260625_Data`.
  - It fits `Im[chi] = slope * Im[beta] + intercept` at each subcircuit step
    and reports the McGarry-style estimate
    `<x> = slope / (measurement_im_chi_sign * sqrt(2))`.
  - For the saved notebook arrays, `measurement_im_chi_sign = -1` because the
    local ideal readout check gives `R q[probe] 0 pi/2` with a down-initialized
    probe as `-Im[chi]`, not `Im[chi]`.
  - It writes `build/readout_x_trace/x_trace_from_slopes.png`,
    `build/readout_x_trace/x_trace_from_slopes.csv`, and
    `build/readout_x_trace/selected_step_slope_fits.png`.

- Added Phil's requested expected-versus-actual overlay for the 1000-repeat data.
  - The overlay compares the slope-extracted data against both exact
    `H_sim` propagation and the ideal compiled-gate simulation inferred from
    `build/dmanh.jaqal`.
  - It writes `build/readout_x_trace/actual_vs_expected_1000_repeats.png`
    and `build/readout_x_trace/actual_vs_expected_1000_repeats.csv`.
  - Under the local simulator/Jaqal convention, the t=0 preparation should be
    `<x> = -1.25895`.
  - The raw 1000-repeat slope extraction gives `slope/sqrt(2) = +0.51790`.
    After the checked readout-sign correction, the inferred t=0 value is
    `<x> = -0.51790`.
  - The corrected t=0 sign agrees with the left-well preparation, but the
    magnitude is still much smaller than the ideal expected value.

- Checked the readout sign and beta-axis mappings explicitly.
  - For the t=0 left-well state, the local ideal characteristic function gives
    `chi(i 0.4) = 0.698749 - 0.603236 i`.
  - With probe initial state `up`, `R(0, +pi/2)` readout returns
    `-0.603236`, matching `Im[chi]`.
  - With probe initial state `down`, `R(0, +pi/2)` readout returns
    `+0.603236`, matching `-Im[chi]`.
  - The saved notebook t=0 data have raw positive slope, so they match the
    sign of the down-probe / `-Im[chi]` convention, not the McGarry Eq. 35
    `Im[chi]` convention.
  - The current local Sandia mapping `beta = -2 i s` with notebook
    `s = -0.2 + i imBeta` predicts an even t=0 curve with zero slope, which is
    not what the saved data show. The saved t=0 curve behaves like a Fig. 3(f)
    `Im[beta]` scan up to the readout-sign correction and a contrast scale.
  - The minus sign in `readout_x_trace.py` is therefore a conversion from the
    saved notebook observable to the McGarry Eq. 35 convention, not a fit fudge:
    if the saved observable is `Z_notebook = -Im[chi]`, then
    `<x> = - slope(Z_notebook) / sqrt(2)`.
  - Files alone do not fully identify whether `Z_notebook = -Im[chi]` comes
    from the probe being effectively in the opposite Z eigenstate, from the
    notebook's `state0 - state1` bit-order convention, or from an opposite
    beta-axis convention. The t=0 left-well calibration diagnoses the net sign,
    but the hardware/logical-bit convention should be confirmed before
    presenting this as a settled hardware readout sign.
  - Teddy confirmed that `xSDF` is the spin-dependent force and is equivalent
    to the conditional-displacement `xCD` gate, with the same argument
    convention. Therefore the saved notebooks' use of `xSDF` rather than
    `xCD` should not be treated as a separate explanation for a sign or phase
    change.
  - `build/notebook/betas.npy` contains the intended McGarry beta scan
    `(0,-0.4)..(0,0.4)`, but the saved result notebooks do not load it; they
    hard-code `reBeta=-0.2` and override `imBeta`. This is another reason to
    keep the readout/beta-axis convention explicit in any report to Phil.

- Audited git history for the suspected negative preparation displacement.
  - Every committed `src/compiler.py` version that contains the explicit
    preparation has `prep_beta = args.x_min / sqrt(2.0)`.
  - No committed `zCD`/`zSDF` preparation line with a literal negative
    displacement argument was found.
  - The older sign-related convention was `alpha_phase_offset = -pi/2` through
    commit `dbfe289`; commit `797ea3f` changed that default to `0` and also
    changed the readout mapping from direct `beta/2` to Sandia argument
    `s = i beta / 2`.
  - The committed Sandia notebooks in `c9953ae` use
    `zSDF q[0] 0.89021208217480396 0` for preparation.

- Checked the notebook evolution-angle sign convention against the generated
  compiler output.
  - `build/dmanh.jaqal` starts evolution step 1 with
    `xCD q[0] 1 1 -0.18512012242326523 0`.
  - `build/notebook/angles.npy` also stores the first exported xSDF argument
    as `[-0.18512012, -0.0]`, and the notebook line
    `angle = (-compiled_angles[:2].T).tolist()` changes that in-memory value to
    `[+0.18512012, +0.0]`.
  - In the current `notebooks/Indiana_DMANH_bosonic_or_v4.ipynb`, that
    in-memory sign flip is canceled by the template starting the evolution
    block with `{m_reBeta}`/`{m_imBeta}`. The printed first evolution pulse is
    therefore negative and matches `build/dmanh.jaqal`.
  - In the saved `results/20260625_Data/.../output.ipynb` notebooks, the older
    template starts the evolution block with `{reBeta}`/`{imBeta}`. With the
    same `angle = -compiled_angles[:2].T` line, those printed circuits start
    evolution step 1 with positive `xSDF q[0] 0.18512012242326523 0.0`, the
    opposite sign from `build/dmanh.jaqal`.
  - This evolution-pulse sign discrepancy can affect evolved subcircuits
    (`t > 0`), but it cannot explain the sign of the zero-evolution
    measurement because the t=0 subcircuit has no evolution pulses.
  - McGarry's theory defines
    `Q(alpha, vartheta, phi)=D(-sigma_x alpha) R_phi(-vartheta) D(sigma_x alpha)`
    and `G_c tilde(alpha)=Q(alpha) Q(-alpha)`. Since Jaqal is listed in
    execution order, the rightmost theoretical factor `Q(-alpha)` must be
    emitted first, giving the pulse signs `-alpha, +alpha, +alpha, -alpha`.
    Thus `build/dmanh.jaqal` and `notebooks/Indiana_DMANH_bosonic_or_v4.ipynb`
    are the McGarry-order version. The saved result notebooks emit
    `+alpha, -alpha, -alpha, +alpha`, i.e. the reversed product
    `Q(-alpha) Q(alpha)`. That reverse order has the same first-order
    continuum target but different finite-`vartheta` Trotter error, so it is
    not the sequence stated in the paper/compiler theory.

- Recovered context after VS Code was accidentally closed and reran the current
  targets.
  - `make dmanh` now regenerates `build/dmanh.jaqal` with the current readout
    mapping: McGarry `beta=i 0.4` emits Sandia `reBeta=-0.2`, `imBeta=0`.
    A stale ignored build artifact had still shown the older `reBeta=0`,
    `imBeta=0.2` mapping.
  - `make readout-x-trace` runs cleanly after regeneration.
  - The regenerated ideal checks report direct compiled-gate
    `<x>=-1.25895` at `t=0`, `<x>=+0.0235949` at step `26`
    (`4.081408 ms`), and `<x>=+1.241030` at step `49`
    (`7.691885 ms`).
  - `make readout-x-trace` now emits expected-versus-actual overlays and CSVs
    for every discovered repeat dataset, not only the 1000-repeat dataset:
    `actual_vs_expected_50_repeats.*`,
    `actual_vs_expected_500_repeats.*`, and
    `actual_vs_expected_1000_repeats.*`.
  - The sign-corrected `t=0` slope estimates are `-0.38386 +/- 0.03819`
    for 50 repeats, `-0.39966 +/- 0.07874` for 500 repeats, and
    `-0.51790 +/- 0.03273` for 1000 repeats, versus the expected
    `-1.25895`.

- Current answer to Phil's follow-up questions:
  1. From the generated Jaqal and local simulator convention, the `t=0`
     preparation is not supposed to be positive. The prep line is
     `zSDF/zCD q[0] ... 0.89021208217480396 0`; for the down branch in the
     local convention this is `<x>=-sqrt(2)*0.890212...=-1.25895`.
     The saved 1000-repeat data have raw `slope/sqrt(2)=+0.51790`, but the
     local ideal readout check says the saved notebook observable has the
     opposite sign from McGarry `Im[chi]`, so the McGarry-converted value is
     `<x>=-0.51790`. Thus the sign can be made consistent with left-well
     preparation only after a net readout-sign conversion; the remaining
     magnitude discrepancy is real. The files alone still do not identify
     whether the net sign comes from probe initial state, `state0-state1`
     bit ordering, or beta-axis convention.
  2. The 1000-repeat expected-versus-actual overlay is
     `build/readout_x_trace/actual_vs_expected_1000_repeats.png`; the numeric
     table is `build/readout_x_trace/actual_vs_expected_1000_repeats.csv`.
     The overlay includes both exact `H_sim` and ideal compiled-gate traces.
  3. The comparison is off immediately at `t=0`: expected is `-1.25895`,
     while the sign-corrected extracted value is `-0.51790 +/- 0.03273`.
     This is not a late-time decoherence-only discrepancy. Evolution-gate
     ordering in the saved result notebooks is also reversed relative to the
     current McGarry-order compiler, but that cannot explain `t=0` because the
     zero-evolution subcircuit has no evolution pulses.

## 2026-06-28

- Agent/model: Codex, GPT-5.

- Read Phil's 2026-06-27 response asking for the expected result under the
  actual saved-measurement condition `reBeta=-0.2` for all measurements.

- Follow-up correction: the requested artifact is the `readout_x_trace`
  expected overlay.
  - Updated `src/readout_x_trace.py` so the `actual_vs_expected_*` outputs
    simulate the same fixed-`reBeta` readout sweep before fitting the expected
    traces. The expected compiled-gate and exact-`H_sim` curves are now:
    simulate `Im chi(reBeta+i imBeta)` at the five saved `imBeta` points,
    fit slope versus the notebook `imBeta` axis, then compare that slope trace
    to the saved data's slope extraction.
  - Defaults for the expected overlay now follow the corrected convention:
    `reBeta=0` and `beta_mapping=direct`. The old saved-notebook condition can
    still be requested explicitly with `--expected-readout-re-beta -0.2`, and
    the legacy local mapping remains available through
    `--expected-beta-mapping sandia-local`.
  - Regenerated `build/readout_x_trace/actual_vs_expected_50_repeats.*`,
    `actual_vs_expected_500_repeats.*`, and
    `actual_vs_expected_1000_repeats.*`.
  - For the 1000-repeat corrected slope trace at steps `0`, `13`, `26`, and
    `49`, the saved data are `[-0.5179,-0.2500,-0.1142,-0.0868]`; the
    compiled-gate fixed-readout expected values are
    `[-1.0730,-0.8116,+0.0092,+1.0183]`; the exact-`H_sim` fixed-readout
    expected values are `[-1.0730,-0.6662,+0.1503,+0.9662]`.
  - The CSV keeps the direct state expectations as reference columns
    `compiled_gate_direct_x` and `exact_hsim_direct_x`, but these are no
    longer the plotted expected curves.
  - After Phil confirmed the readout sign issue, treat the sign conversion as
    a bug fix rather than an alternate "before correction" view. The
    `actual_vs_expected_*` plots now show only the sign-fixed readout slope
    trace. Raw slope columns remain in the CSV for audit/debugging, but are no
    longer plotted as a separate top panel.

- Simplified the compiler's notebook export.
  - `make dmanh` now writes only the angle array needed by the notebook:
    `build/dmanh_angles.npy`.
  - It also writes `build/dmanh_angles.csv` as a human-readable sidecar. CSV
    row `step=k` corresponds to column `k` of `build/dmanh_angles.npy`.
  - The exported angle artifacts now contain only the two per-step SDF/xCD
    displacement coordinates. The constant notebook Rz angle/vartheta is not
    exported.
  - The old `build/notebook/betas.npy` export was removed because the five
    readout beta/imBeta values are fixed in the notebooks and should not be a
    compiler artifact.
  - `experiments/dmanh.toml` now uses `output.angles =
    "build/dmanh_angles.npy"` instead of `output.notebook = "build/notebook"`.
  - The compiler still accepts old directory-style `--export-numpy` /
    `output.notebook` paths as compatibility aliases, but those aliases write
    only `angles.npy` and `angles.csv`, not `betas.npy`.

- Corrected the xCD/readout coordinate convention after checking the paper and
  the Sandia pulse definitions.
  - The McGarry paper requires the readout operation
    `D(sigma_x beta/2)` and the position scan has `Re[beta]=0`.
  - The previous local compiler/simulator convention treated `xCD(s)` as
    `D(-i s)` and therefore emitted `s=i beta/2`; this was not justified by
    the paper or by `Calibration_PulseDefinitions_curated.py`.
  - The compiler now emits direct readout pulse coordinates and the default
    generated template starts at the origin: `let reBeta 0`, `let imBeta 0`.
    Compiler-level readout-coordinate flags were removed so the coordinate
    sweep lives only in the notebook/runner.
  - The ideal simulator and frequency-sweep simulator now interpret xCD
    arguments as direct displacement coordinates.
  - Restored `alpha_phase_offset=-pi/2` in the experiment TOML files so the
    direct-xCD compiled-gate plots continue to match the local
    Schrodinger-frame `H_sim` comparisons.
  - Updated the active Indiana DMANH notebook header to `let reBeta 0`.
  - The saved-data slope overlay keeps the old rotated mapping as a legacy
    option for the already-run `reBeta=-0.2` notebooks.

- Updated the measurement simulation so experiment-facing `chi` and Eq. 33
  density plots execute the parsed readout block from the generated Jaqal.
  - `src/simulator.py` now preserves the final readout gates instead of only
    returning preparation/evolution blocks. It also canonicalizes `xSDF` as
    `xCD` for the semantic simulator.
  - `src/measure.py` now computes plotted `Re[chi]`/`Im[chi]` from ideal
    probe-Z readout experiments on each prefix state. Direct
    `chi(beta)=<D(beta)>` evaluation remains only as a roundoff audit column.
  - `src/readout_x_trace.py` now uses the parsed readout block for expected
    saved-data overlays before converting the ideal result back to McGarry
    `Im[chi]` for slope fitting.

- Added an xSDF-spelled overlay to the `H_sim` comparison plots.
  - `reference/Calibration_PulseDefinitions_curated.py` has separate
    `gate_xSDF` and `gate_xCD` code paths: `gate_xSDF` uses the older direct
    SDF helper, while `gate_xCD` routes through `gate_CD` and the newer
    sideband-manifold/mode/Fock-state calibration plumbing.
  - Both phase-0 paths call the same `find_t_phi_pi2` displacement-angle
    helper after their calibration-specific scaling/sign choices, so the local
    truncated-Fock simulator now plots the xSDF/Rz semantic sequence against
    the xCD/Rz semantic sequence.
  - This is an ideal semantic comparison only. It does not simulate the actual
    red/blue sideband waveforms, Rabi rates, calibration amplitudes, or
    `disp_scale_factor` choices in the curated pulse-definition file.
  - Fixed the semantic parser so notebook-style `xSDF q[0] reBeta imBeta`
    uses the shorter two-coordinate argument layout rather than the
    `xCD q[0] manifold mode reBeta imBeta` layout.

- Added temporary diagnostic hypothesis overlays to `readout_x_trace`.
  - These do not replace the correct exact-`H_sim` or compiled-gate expected
    traces. They are explicitly labeled "hypothesis" curves meant to answer
    Phil's question in reverse: what kind of simulation-side change would be
    needed to resemble the saved repeat data?
  - The generated `actual_vs_expected_*` plots/CSVs now include:
    half-strength preparation only, half-strength preparation plus half
    evolution displacement, and half-strength preparation plus the old saved
    notebook gate order.
  - For the 1000-repeat data at steps `0`, `13`, `26`, and `49`, the saved
    slope trace is `[-0.5179,-0.2500,-0.1142,-0.0868]`. The half-prep-only
    exact-`H_sim` hypothesis gives `[-0.5779,-0.3815,-0.0175,+0.5866]`; the
    half-prep plus half-alpha hypothesis gives
    `[-0.5779,+0.3725,+0.0942,+0.3786]`; and the half-prep plus saved-order
    compiled-gate hypothesis gives `[-0.5779,-0.4878,-0.1234,+0.6399]`.
  - Conclusion so far: a factor-of-two preparation scale can explain the
    reduced `t=0` magnitude, but it does not explain the late-time flatness.
    The old saved-notebook gate order also does not explain the final trace by
    itself. The mismatch therefore still points to an additional readout
    contrast/axis issue, a missing experimental operation, or a more severe
    hardware/runtime discrepancy rather than just the simple SDF scale
    convention.

- Added `src/notebook_observable_audit.py` and `make notebook-observable-audit`.
  - This simulates the saved notebook's observable directly, rather than the
    intended McGarry observable: q0 preparation/evolution, q1 readout,
    `measure_all`, and the notebook postprocessing
    `(prob[0]-prob[2])/(prob[0]+prob[2])`.
  - The audit uses the angle array printed in
    `results/20260625_Data/003_1000 repeats/output.ipynb`, so it follows the
    saved notebook's old `+alpha,-alpha,-alpha,+alpha` order and
    `reBeta=-0.2` readout condition.
  - Under the local ideal xSDF semantics, the notebook as written has the
    wrong raw readout sign relative to the saved `expZ_imMeas` slope. Flipping
    the readout sign gives the correct sign, but the full-strength prep/evolve
    trace is still much too large.
  - For the 1000-repeat data, the saved notebook-extracted slope trace at
    steps `0`, `13`, `26`, and `49` is
    `[-0.5179,-0.2500,-0.1142,-0.0868]`. The sign-flipped half-prep
    notebook-observable hypothesis gives
    `[-0.5664,-0.3495,+0.8261,-0.1328]`; the sign-flipped full-prep with an
    extra half readout scale gives `[-0.6048,-0.5325,+0.5825,+0.0825]`.
    Both can explain the initial scale neighborhood, but both predict a large
    mid-run excursion that is absent from the saved data.
  - The best coarse no-noise fits to the full saved notebook observable require
    evolution to be nearly absent or strongly suppressed, often with an even
    smaller effective readout scale; e.g. a representative weak-evolution
    family gives `[-0.2870,-0.2857,-0.2406,-0.0134]`. This matches the
    flatness better but no longer matches the `t=0` scale.
  - Current interpretation: the first-displacement problem and the late-time
    flatness are not explained by one sign or one factor-of-two convention.
    The saved notebook likely combines at least a readout sign/coordinate
    problem with a separate suppression/misexecution of preparation/evolution
    strength or a probability-bin/postselection observable that is not the
    intended McGarry readout.

- Agent/model: Claude, Opus 4.8.

- Answered Phil's request to re-simulate the expected result under the actual
  saved-measurement condition `reBeta=-0.2` for all measurements.
  - Ran `src/readout_x_trace.py --expected-readout-re-beta -0.2` into
    `build/readout_x_trace_reBeta_-0.2/` (separate dir so the corrected
    `reBeta=0` build is not clobbered). No data was re-taken.
  - Result: `reBeta=-0.2` is essentially identical to `reBeta=0`. The
    1000-repeat exact-`H_sim` expected trace at steps `0`, `13`, `26`, `49`
    is `[-1.0730,-0.6662,+0.1503,+0.9662]` at `reBeta=-0.2` versus
    `[-1.0947,-0.6685,+0.1530,+0.9757]` at `reBeta=0` (~2% difference).
    Physical reason: `reBeta` only shifts the real part of the readout
    displacement, but the slope is extracted along the `Im[beta]` axis, so the
    `Im[chi]` slope is insensitive to `reBeta`. The `reBeta=-0.2` bug is
    therefore NOT the source of the scale/shape mismatch; the data should be
    compared to ~the same curve we already had.

- Clarified what the 50/500/1000 "repeats" datasets are: identical circuit and
  hardware program, differing only in shot count (statistical averaging). Data
  arrays are all shape `(5, 50)` (5 `imBeta` points x 50 timesteps). Implied
  per-point shot noise ~1/sqrt(N): 0.14 (50), 0.045 (500), 0.032 (1000).
  - Consequence: the late-time "flattening getting worse with more repeats" is
    not a physical change. Decoherence acts within a single shot's evolution
    and is independent of repeat count. The flattening is the true signal,
    only resolved once the noise floor drops; the 50-repeat apparent match to
    the swinging hypothesis is noise (late-time error bars +/-0.15..0.27).

- Settled on a realistic two-part hypothesis for the data and added it to the
  `readout_x_trace` expected overlays.
  - Knob 1, "0.5x prep": the preparation conditional displacement delivers half
    the intended amplitude, so the wavepacket starts at `|alpha0|/2`. Since
    `<x> = sqrt(2) Re<a>` and `<a>=alpha`, this halves `<x>` at all times and
    explains the factor-of-2 at `t=0` (data `-0.518` vs ideal `-1.073`, ratio
    0.48). It only rescales amplitude; the shape still swings negative to
    positive, so alone it gives RMS ~0.33 to the 1000-repeat data.
  - Knob 2, `H_sim` damping: added a phenomenological Lindblad dissipator
    `gamma*D[a]` on top of the unitary `H_sim` in a new `damped_hsim_rhos`
    that evolves the postselected motional density matrix via the column-
    stacked Liouvillian. It makes `<a>(t) = alpha0 e^{-i w t} e^{-gamma t/2}`,
    so `<x>(t)` keeps oscillating inside an `exp(-gamma t/2)` envelope that
    decays the swing toward zero -> the late-time flattening. This decay is a
    function of evolution TIME, not repeat count.
  - Honest naming caveat: the implemented term is amplitude damping (energy
    relaxation toward the motional ground state), not literally pure dephasing
    (`D[a^dagger a]`). For this `<x>` observable the two are nearly
    indistinguishable (both give the `exp(-gamma t/2)` envelope); they differ
    only in `<n>`/energy, which this data does not measure. Switching to
    `D[a^dagger a]` would barely change the fitted curve.
  - Default `gamma = 500/s` (time constant `2/gamma ~ 4 ms`, comparable to the
    7.7 ms run), tunable via `--hsim-damping-rate`. Chosen coarsely; "we do
    not need to be so precise." A scan showed full-trace RMS still inching down
    past `gamma ~ 650/s`, but `500/s` is a reasonable round value.
  - Combined "0.5x prep + H_sim damping" full-trace RMS to data: 0.31 (50
    repeats, noise-dominated), 0.10 (500), 0.09 (1000). At steps `0,13,26,49`
    the 1000-repeat data `[-0.518,-0.250,-0.114,-0.087]` versus the combined
    hypothesis `[-0.566,-0.261,-0.040,+0.051]`.

- `src/readout_x_trace.py` changes for this round (kept minimal):
  - New `DEFAULT_HSIM_DAMPING_RATE = 500.0` and `--hsim-damping-rate` CLI flag.
  - New `damped_hsim_rhos(...)` Lindblad propagator; imports
    `annihilation_matrix` from `simulator`.
  - `HypothesisTrace` gained a `damping_rate` field; CSV gained a
    `hypothesis_<i>_damping_rate` column.
  - Diagnostic hypotheses reduced to two: `hypothesis: 0.5x prep` and
    `hypothesis: 0.5x prep + H_sim damping`. The earlier `0.5x prep + 0.5x
    alpha` and `0.5x prep + saved gate order` curves were removed at Phil's
    request, and the damping legend label no longer prints the `gamma` value.

- Still unexplained: the early-time blip/dip, which lives in a different
  quantity (the notebook's directly-plotted observable, not the slope-extracted
  `<x>`) and is untouched by either knob.

Open questions:

- Whether to generate separate hardware programs for `0`, `26`, and `49` evolution steps before readout.
- Whether the `zCD` simulator path should also be routed through a named convention helper, or remain direct because it is being used as a preparation shortcut.
- Whether to clean the appended scratch notes currently present at the end of `README.md`.
- Confirm the hardware sign conventions for `zSDF` preparation, the notebook
  `state0 - state1` postselected `exp_z` convention, and whether notebook
  `reBeta`/`imBeta` should be interpreted as Sandia pulse arguments or
  McGarry characteristic-function coordinates in the slope extraction.
