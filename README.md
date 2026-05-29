# dmanhc - CV-DV Double-Well Compiler

This repository contains a small hardware-facing compiler and ideal simulator for a McGarry-style CV-DV double-well experiment in Sandia/QSCOUT Jaqal vocabulary.

The main target is to generate inspectable Jaqal for a displaced motional wavepacket evolving under a symmetric anharmonic double-well sequence, with optional characteristic-function readout through a probe qubit.

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

Generate the default 4 ms symmetric double-well program:

```sh
python src/compiler.py --output build/symmetric_double_well_4ms_sandia.jaqal
```

Generate plots from the emitted Jaqal:

```sh
python src/plots.py --jaqal build/symmetric_double_well_4ms_sandia.jaqal
```

Run the ideal characteristic-function reconstruction:

```sh
python src/measure.py --jaqal build/symmetric_double_well_4ms_sandia.jaqal
```

Run with parameters:

```sh
python src/compiler.py \
		--output $(DMANH_JAQAL) \
		--max-time-ms 4 \
		--dt-us 200 \
		--delta-hz 754.95 \
		--alpha0 0.49365 \
		--vartheta 0.8 \
		--varphi 0 \
		--x-min 1.208 \
		--probe-qubit-index 0

	$(PYTHON) $(SRC)/plots.py \
		--jaqal $(DMANH_JAQAL) \
		--output $(DMANH_DIRECT_PNG) \
		--title '$(DMANH_DIRECT_TITLE)' \
		--no-hsim-output
```

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
| `--max-time-ms` | total simulated evolution time `T` | `4.0` |
| `--dt-us` | timestep `Delta t` | `200.0` |
| `--delta-hz` | harmonic frequency `delta / 2pi` | `500.0` |
| `--alpha0` | displacement amplitude scale `alpha_0` | `pi / 6` |
| `--vartheta` | cosine-gate angle, with `B = vartheta / Delta t` | `0.8` |
| `--varphi` | cosine phase; currently intended for symmetric `varphi = 0` | `0.0` |
| `--x-min` | initial left-well displacement target | `1.5` |
| `--alpha-phase-offset` | initial phase offset for rotating `alpha_k` | `-pi/2` |

The number of McGarry timesteps is computed internally:

```text
K = T / Delta t
```

For example, `T = 4 ms` and `Delta t = 200 us` gives `K = 20` evolution blocks.

At timestep `k`, the compiler emits

```text
alpha_k = alpha0 * exp(i * (alpha_phase_offset + delta * k * Delta t))
```

This phase advance is part of the McGarry rotating-frame construction. It is not the sideband index and not a Jaqal loop counter.

### Hardware address parameters

| Parameter | Meaning | Default |
|---|---:|---:|
| `--qubit-index` | qubit used for preparation and evolution | `0` |
| `--probe-qubit-index` | qubit used for characteristic-function readout | `0` or `1`, depending on compiler default |
| `--sideband-manifold` | motional manifold address | `1` |
| `--sideband-index` | motional mode address within the manifold | `1` |
| `--nf-start`, `--nf-end` | optional pulse-calibration Fock-state arguments | `0`, `1` |

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




