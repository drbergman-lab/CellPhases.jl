# CellPhases.jl

Give every tracked agent in an agent-based simulation a phase for an oscillating value, and
measure how those phases relate in space.

!!! tierbrief
    Pass per-frame positions, values and IDs with the frame interval and the band of periods to
    keep; `cellPhases` returns a phase per agent per frame. The readouts take that result.

!!! tierfull
    The motivating case is a Notch–Delta model of the developing cortex: every radial glial cell
    cycles its Notch-pathway signals every ~170 minutes, and the questions are how far neighbours
    are in step and whether phase waves travel through the tissue. A pooled correlation answers
    neither well, because opposite local waves cancel when pooled. Following each cell's own phase
    keeps the local structure.

    The pipeline per agent: optionally transform the trace (a log, say), band-pass it to the rhythm of
    interest, and take the angle of its analytic signal (Hilbert transform). The readouts then
    compare phases across agents: coherence with distance, a local order parameter, and local
    phase-gradient waves with their direction and speed. Each is compared with a null from
    time-shift surrogates.

```@example index
using CellPhases, Random

# toy data: 400 slowly moving agents carrying a wave of period 180 that travels +y at
# 0.5 position units per time unit; frames 10 time units apart
include(joinpath(pkgdir(CellPhases), "test", "agent_field.jl"))   # toy agents
field = agentField(; n = 400, L = 300.0, nframes = 200, step = 0.5, rng = Xoshiro(1))
v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 90,
              rng = Xoshiro(2))

p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300))
```

## Reading the result

`p.phase[t][i]` is the phase of the agent in row `i` of frame `t`, in (−π, π], with 0 at the peaks
of the band-passed signal; `NaN` where the agent has none. The readouts summarise it:

```@example index
CycleLengths(p)
```

```@example index
PhaseCoherence(p; dr = 10, rmax = 60, frames = 1:4:200, rng = Xoshiro(3))
```

```@example index
w = PhaseWaves(p; r = 30, frames = 1:4:200)
rose = roseData(w)
(direction = rose.meandirection, speed = rose.speed[argmax(rose.share)],
 clean_percent = rose.clean)
```

## Where to go next

- [Phases](man/phases.md): traces, the band-pass filter and the Hilbert phase.
- [Checking the phase](man/checks.md): does the phase follow the value, and cycle lengths.
- [Synchrony in space](man/synchrony.md): coherence with distance and the local order parameter,
  with their null.
- [Phase waves](man/waves.md): local direction and speed, and the limits of the method.
- [PhysiCell and replicates](man/physicell.md): straight from an output folder or PCMM
  simulations, and comparing conditions.
- [Plotting](man/plotting.md): plot-ready data and recipes.
- [Synthetic data](man/synthetic.md): generators with known answers.
- [Pitfalls](man/pitfalls.md): mistakes already made once.
- [API reference](reference.md).
