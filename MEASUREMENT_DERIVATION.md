# Characteristic-Function Readout Derivation

This note records the hand calculation for the McGarry characteristic-function
readout and how it maps to the compiler/simulator code.

## Definitions

The mathematical oscillator displacement used by the simulator is

```text
D(alpha) = exp(alpha a^\dagger - alpha^* a).
```

McGarry defines the motional characteristic function as

```text
chi(beta) = <D(beta)> = Tr[rho D(beta)].
```

The readout circuit measures `chi(beta)` by applying the conditional
displacement

```text
D(sigma_x beta / 2)
```

to a probe qubit and then measuring the probe in the Z basis.

## xCD Convention

The convention is written directly at the relevant call sites:
`src/simulator.py` converts evolution `xCD` arguments, `src/compiler.py`
maps readout `beta` to emitted Sandia variables, and `src/measure.py` builds
the ideal mathematical readout unitary.

Current Sandia/DMANH convention:

```text
math_alpha = -i * sandia_alpha
```

Previous local mathematical convention:

```text
math_alpha = sandia_alpha
```

Therefore, to implement a desired mathematical displacement `D(gamma)` with a
Sandia `xCD` argument `s`, use

```text
s = i * gamma.
```

For McGarry readout, `gamma = beta / 2`, so

```text
s = i * beta / 2.
```

For the default 2PFD x-readout point,

```text
beta = i h = i 0.4,
gamma = i 0.2,
s = i * i 0.2 = -0.2.
```

Thus the generated Jaqal should contain

```text
let reBeta -0.2
let imBeta 0
let imMeas 1
```

where `imMeas=1` prepends the X-axis `R q[probe] 0 pi/2` pulse that selects
the imaginary part.

## Probe Expectation Values

Let

```text
gamma = beta / 2
U = |+><+| D(gamma) + |-><-| D(-gamma),
```

where `|+>` and `|->` are the probe qubit's `sigma_x` eigenstates. In the
`sigma_x` basis,

```text
sigma_z = |+><-| + |-><+|.
```

After the conditional displacement,

```text
U^\dagger sigma_z U
  = |+><-| D(-beta) + |-><+| D(beta).
```

For an initial probe state `|up> = (|+> + |->) / sqrt(2)`,

```text
<sigma_z>
  = (chi(beta) + chi(-beta)) / 2
  = Re[chi(beta)],
```

using Hermitian symmetry `chi(-beta) = chi(beta)^*`.

For an initial probe state `R_x(theta)|up>`,

```text
<sigma_z>
  = (e^{-i theta} chi(beta) + e^{i theta} chi(-beta)) / 2.
```

At `theta = pi/2`,

```text
<sigma_z> = Im[chi(beta)].
```

This is the McGarry statement: no prepended pulse measures the real part, and
prepended `R_x(pi/2)` measures the imaginary part.

The sign above assumes the simulator's probe basis convention:

```text
|up> = [1, 0], |down> = [0, 1], sigma_z = |up><up| - |down><down|.
```

If the same circuit starts the probe in `|down>` instead of `|up>`, both signs
flip:

```text
no rotation:      <sigma_z> = -Re[chi(beta)]
R_x(pi/2):        <sigma_z> = -Im[chi(beta)].
```

The QSCOUT notebook computes a postselected `state0 - state1` probability
difference. Its comments identify `state1` as probe spin-up in the AJC
postselection helper, so the notebook's reported sign may correspond to
`P_down - P_up` rather than the simulator's `P_up - P_down`. That is a basis
labeling issue, not a displacement issue; when comparing with hardware data,
check whether the reported expectation is `<Z>` or `-<Z>` in the simulator's
basis.

## Position Expectation

For the one-dimensional position readout used by 2PFD,

```text
<x> = (1 / sqrt(2)) d Im[chi(beta)] / d Im[beta] |_{beta=0}.
```

The two-point finite-difference approximation with `beta = i h` is

```text
<x> ~= Im[chi(i h)] / (sqrt(2) h).
```

For the default `h=0.4`, the compiler emits the mapped Sandia readout argument
`s=-0.2` and `imMeas=1`, so the measured probe expectation estimates
`Im[chi(i 0.4)]`.

With the simulator's positive sign convention,

```text
<x>_2PFD = probe_z_rx_pi_2 / (sqrt(2) h).
```

If a hardware analysis reports `P_down - P_up`, use the negative of that
reported value in the same formula.

## DMANH+ Hand-Check Rows

For the generated `build/dmanh.jaqal` readout:

```text
requested McGarry coordinate: beta = i 0.4
desired conditional displacement: D(sigma_x i 0.2)
Sandia xCD argument: s = i beta / 2 = -0.2
generated variables: reBeta = -0.2, imBeta = 0, imMeas = 1
```

The following rows are the direct hand-check quantities from the same emitted
Jaqal, using the postselected `down` motional branch and the simulator's
`|up>` probe sign convention.

| time ms | chi(i0.4) | probe no rotation | probe Rx(pi/2) | direct <x> | 2PFD <x> |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000000 | 0.698749069671 - 0.603235879735 i | 0.698749069671 | -0.603235879735 | -1.258950000000 | -1.066380453040 |
| 4.081408 | 0.750814044745 + 0.003615150741 i | 0.750814044745 | 0.003615150741 | 0.023601780168 | 0.006390744009 |
| 7.691885 | 0.616784476273 + 0.564369115177 i | 0.616784476273 | 0.564369115177 | 1.241030047270 | 0.997673071085 |

These rows verify the measurement identity used by `src/measure.py`:

```text
probe_z_no_rotation = Re[chi(i0.4)]
probe_z_rx_pi_2     = Im[chi(i0.4)].
```

The gap between `direct <x>` and `2PFD <x>` is the finite-difference error at
`h=0.4`, not a readout sign or compiler mapping error.

## Code Map

- `src/compiler.py`
  - Computes the emitted readout `reBeta/imBeta` directly as `i beta / 2`.
  - Emits `R q[probe] 0 pi/2` when `imMeas=1`.

- `src/simulator.py`
  - Applies evolution `xCD` by converting each Sandia argument `s` to
    mathematical displacement parameter `-i s`.

- `src/measure.py`
  - Computes direct `chi(beta) = Tr[rho D(beta)]`.
  - Builds the ideal readout unitary as mathematical `D(sigma_x beta/2)`.
  - Prints diagnostic columns comparing no-rotation readout to `Re[chi]` and
    X-rotated readout to `Im[chi]`.
