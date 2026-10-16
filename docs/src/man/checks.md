# Checking the phase

Two checks to run before reading anything into the phases: does the phase follow the value, and
how long are the cycles.

```@setup checks
using CellPhases, Random
include(joinpath(pkgdir(CellPhases), "test", "agent_field.jl"))   # toy agents
field = agentField(; n = 200, L = 200.0, nframes = 150, rng = Xoshiro(1))
v = pulseTrain(field.ids; dt = 10, period = 170, rng = Xoshiro(2))
p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300))
```

## Does the phase follow the value?

!!! tierbrief
    `PhaseValueCheck` regresses each agent-frame's transformed value, z-scored within its track,
    on `cos φ` and `sin φ`.

!!! tierfull
    `r2` is the share of the value's variance explained by the phase alone; the rest is
    fluctuation outside the band, mostly bursts. `peak` is the phase at which the fitted value is
    highest. For a smooth oscillation it sits at 0, where the band-passed signal peaks. A pulse
    or a signal with a broad plateau and a sharp trough still has a clear relation to the phase,
    but `peak` can move off 0 and the binned medians (`centers`, `median`, `q25`, `q75`) show the
    shape.

    In the motivating cortex model, the phase explained 44% of notch and NICD variance and 63–71%
    of the pulse-like Delta and Ngn.

```@example checks
PhaseValueCheck(p)
```

## Cycle lengths

!!! tierbrief
    A cycle ends each time φ crosses 0 going upward; `CycleLengths` returns the time between
    successive crossings, per agent.

!!! tierfull
    Crossing times are interpolated linearly between frames, so cycle lengths are not quantised
    to the frame interval. Crossings are counted only within unbroken stretches of phase, so a
    masked frame never merges two cycles into one. A distribution centred inside the band with
    room on both sides is the sign that the band brackets the rhythm; one piled against a band
    edge means the band is cutting it.

```@example checks
c = CycleLengths(p)
```

!!! tierdev
    A step from φ = −3 to φ = +3 is a wrap, not an upward crossing of 0; crossings require
    `b − a < π`.
