# Synchrony in space

How far apart agents can be and still be in step, measured two ways, each against a null.

```@setup sync
using CellPhases, Random
field = agentField(; n = 400, L = 300.0, nframes = 200, step = 0.5, rng = Xoshiro(1))
v = oscillatorField(field.positions, field.ids; dt = 10, period = 180,
                    correlationlength = 25.0, phasesd = 1.0, rng = Xoshiro(2))
p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300))
```

## The null: time-shift surrogates

!!! tierbrief
    `timeShift(p)` rotates each track's phase series in time by its own random offset. Every
    agent keeps its rhythm; agents stop being aligned with each other.

!!! tierfull
    The null asks: how in step would these agents look by chance, given their own rhythms? If
    the whole tissue is partly synchronised at some moment, even distant pairs score above 0, and
    a single agent's slow drift in frequency makes it look partly aligned with anything. Shifting
    each track as a whole keeps all of that and destroys only the alignment between agents.

    **Why not shuffle phases among agents within each frame?** That keeps the mean right but
    makes the null band far too narrow. A fresh shuffle in each frame is independent of the
    shuffle in the next, so averaging over frames shrinks its spread as if every frame were an
    independent observation. The real phases are not: they are band-passed, so each agent's
    phase moves only about `2π·dt/period` per frame and consecutive frames are correlated. The
    observed value therefore varies much more than the shuffled one, and p-values come out much
    too small.

    `PhaseCoherence` and `OrderParameter` draw `nsurrogates` surrogates (19 by default, so the
    smallest p-value is 0.05) and report their mean as `null`, their 2.5–97.5% range as
    `nulllow`/`nullhigh`, and a one-sided `pvalue`. Raise `nsurrogates` for finer p-values; the
    cost grows linearly.

!!! tierjournal "2026-10-15 — Time-shift surrogates as the null"
    Decided: the null for coherence and the order parameter is the time-shift surrogate, drawn
    per track.
    Rejected: the prototype's per-frame shuffle of phases among agents. Its mean is a fine
    baseline to subtract, but its shuffles are independent from frame to frame while the real
    phases are correlated across frames, so its band is far too narrow for p-values.

## Coherence with distance

!!! tierbrief
    `PhaseCoherence` averages `cos(φᵢ − φⱼ)` over pairs of agents in the same frame, binned by
    the distance between them.

!!! tierfull
    `cos(φᵢ − φⱼ)` is 1 for two agents at the same point in their cycle, 0 a quarter cycle apart
    and −1 opposite. Bins are `dr` wide and centred on `0, dr, …, rmax`. The first bin holds only
    the closest pairs; in PhysiCell these are compressed or newly divided cells and are real.

    `excess = coherence − null` is the coherence due to proximity. `dr` and `rmax` are in the
    positions' units and have no defaults. `frames` selects which frames to pool. Because the
    phases are band-passed, consecutive frames are correlated when the band's shortest period
    spans many frames (with `dt = 10` and `band = (100, 300)`, a phase moves at most about 36°
    per frame), so `frames = 1:2:n` or sparser saves time for little loss. Compare two
    selections to check.

```@example sync
c = PhaseCoherence(p; dr = 10, rmax = 80, frames = 1:4:200, rng = Xoshiro(3))
```

## The local order parameter

!!! tierbrief
    `OrderParameter` measures, around each agent, how aligned the phases within radius `r` are:
    `R = |mean e^{iφ}|`, averaged over agents and frames.

!!! tierfull
    Each phase is a unit arrow; `R` is the length of the average arrow over the agents within
    `r`, the centre agent included. `R = 1` means all point the same way. Random phases do not
    give 0: the average of `n` random unit arrows has length about `1/√n`, so with few neighbours
    at small `r`, chance alone gives a large `R`. Read `R` against `null`. `R` below `null` means
    neighbours are less aligned than chance, as in lateral inhibition.

```@example sync
o = OrderParameter(p; radii = [10, 20, 40, 60], frames = 1:4:200, rng = Xoshiro(3))
```

!!! tierdev
    Both readouts draw all surrogates up front (`_phasesets`) and accumulate observed and
    surrogate sums in the same pass over pairs, so the pair search runs once per frame. Memory is
    `(1 + nsurrogates)` copies of the phase arrays.
