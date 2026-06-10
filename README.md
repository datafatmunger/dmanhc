# dmanhc - CV-DV DMANH Compiler

This repository contains a small hardware-facing compiler and ideal simulator for Phil's DMANH+ double-well experiment in Sandia/QSCOUT Jaqal vocabulary.

The main target is to generate inspectable Jaqal for a displaced motional wavepacket evolving under the DMANH+ fit parameters. A McGarry-style symmetric double-well benchmark remains available for comparison and uses the same parameter interface.

## Layout

```text
src/
  gates.py       # Jaqal gate-text helpers
  compiler.py    # emits hardware-facing Jaqal
  simulator.py   # ideal truncated-Fock gate semantics
  measure.py     # ideal characteristic-function readout utilities
  plots.py       # plotting helpers and diagnostics
build/
  *.jaqal        # generated Jaqal programs
  *.png          # generated plots
```

## Basic usage

Install dependencies:

```sh
python -m pip install -r requirements.txt
```

Generate Phil's DMANH+ Jaqal, density plot, `H_sim` trace, and readout/chi plots:

```sh
make dmanh
```

Generate the symmetric double-well benchmark and diagnostics:

```sh
make plots
```

Run the direct compiler/plot commands for Phil's DMANH+ target:

```sh
python src/compiler.py \
  --output build/dmanh.jaqal \
  --steps 49 \
  --B-rad-s 5.09628e3 \
  --delta-rad-s 1.29817e3 \
  --alpha0 0.18512 \
  --vartheta 0.8 \
  --x-min 1.25895 \
  --alpha-phase-offset 0

python src/plots.py \
  --jaqal build/dmanh.jaqal \
  --output build/dmanh.png \
  --title 'DMANH+' \
  --times-ms 0 4.081408 7.691885 \
  --hsim-output build/dmanh_hsim.png \
  --hsim-title 'DMANH+ exact $H_{\mathrm{sim}}$ versus compiled-gate dynamics' \
  --hsim-max-time-ms 7.6918850612603702

python src/measure.py \
  --jaqal build/dmanh.jaqal \
  --output build/dmanh_measurement_panels.png \
  --title 'DMANH+ McGarry Eq. 33 readout from chi(beta)' \
  --times-ms 0 4.081408 7.691885 \
  --chi-output build/dmanh_chi_slice_panels.png \
  --chi-title 'DMANH+ characteristic-function slice'
```

Both Makefile experiment targets use the same physical parameter form:

| Target | `K` steps | `B` rad/s | `delta` rad/s | `alpha0` | `x_min` | snapshots ms |
|---|---:|---:|---:|---:|---:|---|
| `make dmanh` | `49` | `5.09628e3` | `1.29817e3` | `0.18512` | `1.25895` | `0, 4.081408, 7.691885` |
| `make plots` | `20` | `4.0e3` | `3.141592653589793e3` | `0.5235987755982989` | `1.5` | `0, 2.00, 4.00` |

## Generated Jaqal structure

The compiler emits a program with four conceptual parts:

```text
prepare_all
zCD preparation
xCD/Rz McGarry evolution blocks
optional probe readout
measure_all
```

The current symmetric evolution block uses only `xCD` and `Rz`:

```text
xCD
Rz
xCD
xCD
Rz
xCD
```

Preparation uses `zCD` so the initially prepared qubit, assumed to be in a Z-basis state, gives an ordinary displacement of the motional wavepacket rather than an X-basis conditional cat state.

Readout, when included, follows the QSCOUT characteristic-function pattern: an optional probe-qubit X rotation followed by a Sandia `xCD` argument mapped from McGarry's mathematical `D(sigma_x beta/2)` and Z-basis measurement. In Jaqal the rotation may be emitted as the notebook-style generic rotation

```text
R q[probe] 0 theta
```

where `theta = pi/2` selects the imaginary part of the characteristic function. `Rz` is used for the symmetric McGarry evolution, not for the Re/Im readout selector.

## Core parameters

### Physics parameters

| Parameter | Meaning | Default |
|---|---:|---:|
| `--steps` | integer McGarry timestep count `K`; derives `T = K Delta t` | unset |
| `--B-rad-s` | cosine amplitude `B`, with `Delta t = vartheta / B` | unset |
| `--delta-rad-s` | angular harmonic frequency `delta` | unset |
| `--alpha0` | displacement amplitude scale `alpha_0` | `pi / 6` |
| `--vartheta` | cosine-gate angle, with `B = vartheta / Delta t` | `0.8` |
| `--x-min` | initial left-well displacement target | `1.5` |
| `--alpha-phase-offset` | optional extra phase added to McGarry `zeta = k Delta t`; default follows Sandia/DMANH direct xCD convention | `0` |

The Makefile targets use angular units directly. The compiler still accepts `--max-time-ms`, `--dt-us`, and `--delta-hz` for manual compatibility runs, but the experiment targets avoid those aliases so the fit-to-compiler map stays explicit.

The compiler can either compute the number of McGarry timesteps from a requested total time,

```text
K = T / Delta t
```

or derive the total time from an explicit integer `K`:

```text
T = K Delta t
```

For Phil's DMANH+ fit, `B = 5.09628e3 rad/s` and `vartheta = 0.8` imply
`Delta t = 156.977246 us`. With `--steps 49`, the total evolution time is
`7.691885 ms`, with a split-like intermediate snapshot at step 26
(`4.081408 ms`).

At timestep `k`, the compiler emits

```text
alpha_k = alpha0 * exp(i * (alpha_phase_offset + delta * k * Delta t))
```

This phase advance is part of the McGarry rotating-frame construction. It is not the sideband index and not a Jaqal loop counter.

### Hardware address parameters

| Parameter | Meaning | Default |
|---|---:|---:|
| `--qubit-index` | qubit used for preparation and evolution | `0` |
| `--probe-qubit-index` | qubit used for characteristic-function readout | `1` |
<!-- | `--sideband-manifold` | motional manifold address | `1` |
| `--sideband-index` | motional mode address within the manifold | `1` |
| `--nf-start`, `--nf-end` | optional pulse-calibration Fock-state arguments | `0`, `1` | -->

These are hardware addressing/calibration values. They do not change the abstract McGarry timestep count `K`.

### Readout parameters

| Parameter | Meaning | Default |
|---|---:|---:|
| `--readout-re-beta` | real part of requested characteristic-function coordinate | `0.0` |
| `--readout-im-beta` | imaginary part of requested characteristic-function coordinate | `0.4` |
| `--imaginary-readout-loop` | when `1`, prepend `R q[probe] 0 pi/2` before readout | `1` |

The generated readout `xCD` is mapped from half of the requested characteristic-function coordinate because the probe protocol applies `D(sigma_x beta / 2)`. With the current Sandia/DMANH convention, the Sandia argument `s` is converted to the mathematical displacement `D(-i s)`.

## Gate conventions

The compiler distinguishes three roles:

```text
preparation: zCD
symmetric evolution: xCD + Rz
readout/probe: R q[probe] 0 theta + xCD
```

The ideal simulator treats `xCD`, `zCD`, and `Rz` as abstract gates on a truncated oscillator Hilbert space. It does not simulate red/blue sideband tones, Rabi rates, Lamb-Dicke factors, pulse calibration, mode crosstalk, heating, or measurement noise.

The sideband manifold and mode index are passed through to the generated Jaqal as hardware addresses. They are not used by the ideal simulator except to parse the displacement arguments.







My read, item by item:

  1. Initial zCD displacement
     Current build/dmanh.jaqal:23 uses 0.85418499167334927, which is 1.208 / sqrt(2). That comes directly from src/compiler.py:91.

     But if I recompute the well minimum from the DMANH potential parameters, I get:
     x_min ≈ 1.2589101853, hence x_min / sqrt(2) ≈ 0.8901839289.

     So if Overleaf section 3.1 derives the minimum from the fitted potential, then yes, our magnitude is likely stale/wrong. There is also a separate
     sign-convention question for zCD.

  2. Elapsed time t=0
     Feedback is correct if t=0 means “prepare, then immediately measure.” That circuit should have no evolution displacements.

     The generated comment // step 0: t = 0 ms currently means “first Trotter interval, evaluating the potential at time 0,” not “measurement at elapsed
     time 0.” That first interval does contain xCD gates. We probably need separate generated circuits for measurement after 0, 26, and 49 steps, rather
     than one final-time program.

  3. First xCD is imaginary
     This is due to our explicit alpha_phase_offset = -pi/2. At step 0:

     alpha = alpha0 * exp(i * (-pi/2)) = -i * 0.18512

     and the first emitted gate uses -alpha, so it becomes +i * 0.18512.

     If we follow the paper’s literal zeta = k delta Delta t with no local offset, the first emitted xCD would be real, roughly -0.18512 0. Also, B Delta t
     = 0.8 is vartheta, not zeta; the phase advance is delta Delta t ≈ 0.203783.

  So the likely real issues are: x_min magnitude, separate zero-evolution measurement circuits, and whether our -pi/2 phase offset is the convention they
  want.
