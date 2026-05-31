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

Generate Phil's DMANH+ Jaqal and density plot:

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
  --x-min 1.208 \
  --alpha-phase-offset -1.5707963267948966

python src/plots.py \
  --jaqal build/dmanh.jaqal \
  --output build/dmanh.png \
  --title 'Phil DMANH' \
  --times-ms 0 4.081408 7.691885 \
  --no-hsim-output
```

Both Makefile experiment targets use the same physical parameter form:

| Target | `K` steps | `B` rad/s | `delta` rad/s | `alpha0` | `x_min` | snapshots ms |
|---|---:|---:|---:|---:|---:|---|
| `make dmanh` | `49` | `5.09628e3` | `1.29817e3` | `0.18512` | `1.208` | `0, 4.081408, 7.691885` |
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

Readout, when included, follows the QSCOUT characteristic-function pattern: an optional probe-qubit X rotation followed by `xCD(beta/2)` and Z-basis measurement. In Jaqal this may be emitted as the notebook-style generic rotation

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
| `--alpha-phase-offset` | local gate-phase correction aligning the first cosine kick with plotted `x` | `-pi/2` |

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

The generated readout `xCD` uses half of the requested characteristic-function coordinate because the probe protocol applies `D(sigma_x beta / 2)`.

## Gate conventions

The compiler distinguishes three roles:

```text
preparation: zCD
symmetric evolution: xCD + Rz
readout/probe: R q[probe] 0 theta + xCD
```

The ideal simulator treats `xCD`, `zCD`, and `Rz` as abstract gates on a truncated oscillator Hilbert space. It does not simulate red/blue sideband tones, Rabi rates, Lamb-Dicke factors, pulse calibration, mode crosstalk, heating, or measurement noise.

The sideband manifold and mode index are passed through to the generated Jaqal as hardware addresses. They are not used by the ideal simulator except to parse the displacement arguments.

