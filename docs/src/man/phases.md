# Phases

`cellPhases` turns per-frame positions, values and IDs into a phase for every tracked agent.

!!! tierbrief
    Frames must be equally spaced; `dt` is the interval in minutes. `band = (shortest, longest)`
    is the range of periods, in minutes, that holds the rhythm.

```@example phases
using CellPhases, Random
field = agentField(; n = 200, L = 200.0, nframes = 150, rng = Xoshiro(1))
v = pulseTrain(field.ids; dt = 10, period = 170, rng = Xoshiro(2))

p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300))
```

The input is three per-frame vectors: `positions[t]` (tuples, vectors or `SVector`s), `values[t]`
and `ids[t]`, one entry per agent present in frame `t`. IDs can be of any type.

!!! tierjournal "2026-10-15 — A separate package from SpaceTimeCorrelations, no shared base"
    Decided: CellPhases is its own package, with no common base package.
    Why: SpaceTimeCorrelations computes the pooled correlation from per-frame values and needs no
    IDs and no oscillation. CellPhases needs tracked agents with long tracks and one dominant
    rhythm, and pulls in DSP.jl. The two share input conventions only (per-frame positions,
    values, IDs; PhysiCell snapshots). Analyses that use both, such as checking that the pooled
    correlation peak does not drift, belong in analysis scripts.

## 1. Each agent's trace

!!! tierbrief
    Each agent's frames are collected in time order. A missing frame splits the trace; pieces
    shorter than `minlength` minutes are dropped.

!!! tierfull
    The filter in step 2 needs evenly spaced samples, so a gap cannot be bridged. In PhysiCell
    output a gap appears when a cell dies, or changes into a cell type that was not selected.

    The default `minlength` is `2trim + 2·band[2]`: after trimming both ends (step 3), two of the
    slowest cycles in the band remain. With the default band of the motivating case
    (100–300 min) and 10-minute frames that is 960 minutes, 96 frames.

    The trace is transformed (`transform = log` by default) and its mean subtracted. The log
    keeps bursts from dominating a heavy-tailed signal: a jump from 2 to 20 counts the same as
    one from 1 to 10. With `log`, values must be positive; pass `transform = identity` for a
    signal that can be zero or negative.

Each kept piece is a [`Track`](@ref):

```@example phases
p.tracks[1]
```

!!! tierdev
    Tracks are built in `src/phases.jl`: `_byagent` gathers `(frames, rows)` per ID, `_runs`
    splits at gaps. A track's `frames` and `rows` cover only the frames with a phase, after
    trimming; surrogates and the readouts walk tracks, never IDs, so an agent with two tracks is
    two independent series.

## 2. Band-pass filter

!!! tierbrief
    A zero-phase Butterworth band-pass keeps periods between `band[1]` and `band[2]` minutes.

!!! tierfull
    The phase in step 3 is only meaningful for one dominant rhythm. Fast bursts would make the
    angle jump, and slow drift would move the analytic signal off the origin. The filter removes
    both: with the default order 3 and band 100–300 minutes, a 60-minute burst component and a
    500-minute drift are both cut to about 1% of their amplitude, while 150–250-minute periods
    pass at more than 90%.

    It runs forward and then backward (`filtfilt`). The two passes cancel the filter's delay, so
    peaks stay where they were, and they square the gain: the cutoffs, the half-power points of
    one pass, become half-amplitude points.

    **What a Butterworth filter is.** Each filtered value is a weighted sum of recent raw values
    minus feedback from recent filtered values,
    `y[n] = b₀x[n] + … + b₆x[n−6] − a₁y[n−1] − … − a₆y[n−6]`. The family (Butterworth: the
    flattest possible pass band, no ripple), the order and the cutoffs fix all the coefficients;
    a band-pass doubles the order, to 6 for order 3. The `b`s block constants outright, which is
    why each agent's baseline disappears. The `a`s set three pole pairs; for 100–300 minutes at
    10-minute frames they resonate at periods of about 108, 211 and 295 minutes and together make
    the flat top.

    **Why order 3, applied once.** Steeper filters, or repeated passes, narrow the band, ring for
    longer, and impose a clean rhythm on irregular agents, hiding the cycle-to-cycle variability
    that is often the point.

```@example phases
using DSP
f = bandpassFilter((100, 300); dt = 10)
round.(abs.(freqresp(f, 2π * 10 .* [1/60, 1/100, 1/170, 1/300, 1/500])); digits = 3)
```

!!! tierdev
    `bandpassFilter` always passes `fs = 1/dt` to `DSP.digitalfilter`, so the cutoffs are in
    cycles per minute. Without `fs`, DSP reads them as fractions of the Nyquist frequency. The
    `filter` test pins the half-power points at the band edges, which is the check that catches
    a missing `fs`; a pure sine does not (see [Pitfalls](pitfalls.md)).

!!! tierjournal "2026-10-15 — Physical units only; fs always passed"
    Decided: the API takes the band, `trim` and `minlength` in minutes and `dt` in minutes per
    frame; frames appear only in `frames` selections and in results. `bandpassFilter` always
    passes `fs`.
    Why: the prototype once built the filter without `fs`; DSP read the cutoffs as fractions of
    Nyquist and passed 200–600 minutes instead of 100–300. The spurious "upward wave" it produced
    survived every pure-sine check. Taking minutes everywhere leaves no frame-based number for a
    user to mistranslate.

## 3. Hilbert phase

!!! tierbrief
    The Hilbert transform turns the band-passed trace `y` into a point `z = y + i·H[y]` circling
    the origin; the phase is its angle and the amplitude its distance from the origin.

!!! tierfull
    `y(t)` alone is ambiguous: the same value occurs on the way up and on the way down. `H[y]` is
    `y` shifted by a quarter cycle (each frequency component by a quarter of its own period), so
    for `y = A cos φ`, `H[y] = A sin φ`. The phase is `φ = arg z` in (−π, π], with φ = 0 at the
    peaks of `y`; no frequency is assumed, and the frequency comes out as the rate at which φ
    advances.

    The ends are trimmed: both the filter and `H` use neighbouring times, so the first and last
    `trim` minutes of every track are dropped. `trim` defaults to the filter's memory, the time
    for its slowest pole to decay to 1/e, rounded up to whole frames: 180 minutes for 100–300
    minutes at 10-minute frames.

    **Masking.** Where the amplitude is close to zero the angle is mostly noise. `mask = q` sets
    the phase to `NaN` in each track's frames whose amplitude is below that track's
    `q`-quantile; `mask = 0.2` drops the weakest fifth. It is off by default.

```@example phases
p.trim, p.minlength
```

!!! tierdev
    `hilbert` is `DSP.hilbert` (FFT, zero the negative frequencies, double the positive, inverse
    FFT). The trim default comes from `_filtermemory(filter)`, `−1/log(max |pole|)` frames.

!!! tierjournal "2026-10-15 — trim and minlength derived from the band"
    Decided: `trim` defaults to the filter's memory and `minlength` to `2trim + 2·band[2]`, both
    in minutes.
    Why: the prototype fixed 18 and 90 frames for 10-minute frames and a 100–300-minute band.
    Both are really properties of the filter: 18 frames was "about the filter's memory", and 90
    frames "five cycles plus trim". Deriving them keeps them right when the band or the frame
    interval changes. For the prototype's settings the derived values are 180 minutes (18
    frames, as before) and 960 minutes (96 frames rather than 90).
