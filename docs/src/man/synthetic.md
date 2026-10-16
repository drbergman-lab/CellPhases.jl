# Synthetic data

Generators of values with known answers, in the shape `cellPhases` takes, for checking the methods
and their settings on your own agents.

!!! tierbrief
    `planeWave`, `pulseTrain` and `oscillatorField` make values with a known answer on positions
    and IDs you supply, usually a real run's, to check the method and its settings on that
    geometry. All values are positive, so `transform = log` also applies.

On a real run, keep its agents and replace the signal:

```julia
real = cellPhases(sequence; value = :notch, cell_types = "rgc", band = (100, 300))
v = planeWave(real.positions; dt = 10, period = 180, speed = 25/60, direction = 90)
w = PhaseWaves(cellPhases(real.positions, v, real.ids; dt = 10, band = (100, 300)); r = 30)
# recovered speed and direction, and cleanShare(w), say whether r = 30 works on these cells
```

The examples below use toy agents from the test suite instead.

```@example synthetic
using CellPhases, Random
include(joinpath(pkgdir(CellPhases), "test", "agent_field.jl"))   # toy agents
field = agentField(; n = 300, L = 250.0, nframes = 150, step = 0.5, rng = Xoshiro(1))
```

## Plane waves: known speed and direction

!!! tierfull
    `exp(baseline + cos(ωt − k·x) + noise·ε)`, travelling at `speed` position units per time
    unit in `direction` degrees. The phase `cellPhases` should return is `ωt − k·x`, which is what the
    sign-convention tests check.

```@example synthetic
wave = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = -90,
                 rng = Xoshiro(2))
w = PhaseWaves(cellPhases(field.positions, wave, field.ids; dt = 10, band = (100, 300)); r = 30)
roseData(w).meandirection
```

## Pulse trains: a broadband signal

!!! tierfull
    Each agent fires a brief pulse once per cycle, with its own slightly different period and a
    random start. Most of a pulse train's power is in harmonics of the period, so it tests the
    band-pass filter in a way a sine cannot: a sine keeps its phase when attenuated, while a
    pulse train read through the wrong band gives a harmonic's period instead of the true one.

```@example synthetic
pulses = pulseTrain(field.ids; dt = 10, period = 170, rng = Xoshiro(3))
CycleLengths(cellPhases(field.positions, pulses, field.ids; dt = 10, band = (100, 300)))
```

## An oscillator field: known coherence

!!! tierfull
    One period for all agents, with phase offsets drawn from a Gaussian random field of standard
    deviation `phasesd` and correlation length `correlationlength`, so the expected coherence at
    distance `d` is `exp(−phasesd²(1 − exp(−d²/2ℓ²)))`. With `correlationlength = 0` the offsets
    are independent: the null case. `diffusion` adds slow phase noise.

```@example synthetic
osc = oscillatorField(field.positions, field.ids; dt = 10, period = 180, rng = Xoshiro(4))
c = PhaseCoherence(cellPhases(field.positions, osc, field.ids; dt = 10, band = (100, 300));
                   dr = 10, rmax = 40, frames = 1:4:150, rng = Xoshiro(5))
c.excess
```

!!! tierdev
    The measured coherence of an `oscillatorField` falls short of the formula at small distances:
    the value noise and the filter both add phase noise, and the offset is fixed at the agent's
    first position while it moves. Tests check ordering and the null, not the formula.
