# Phases

`cellPhases` turns per-frame positions, values and IDs into a phase for every tracked agent.

!!! tierbrief
    Frames must be equally spaced, `dt` apart. `band = (pmin, pmax)` is the range of periods
    that holds the rhythm, in the same time unit as `dt`.

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
    shorter than `minlength` are dropped.

!!! tierfull
    The filter in step 2 needs evenly spaced samples, so a gap cannot be bridged. In PhysiCell
    output a gap appears when a cell dies, or changes into a cell type that was not selected.

    The default `minlength` is `2trim + 2·band[2]`: after trimming both ends (step 3), two of the
    slowest cycles in the band remain. With `dt = 10` and `band = (100, 300)` that is 960,
    96 frames.

    The trace is passed through `transform` (`identity` by default) and its mean subtracted.
    For a heavy-tailed, bursty signal, `transform = log` keeps bursts from dominating: a jump
    from 2 to 20 then counts the same as one from 1 to 10. With `log`, values must be positive.
    The transform also sets what `PhaseValueCheck` compares the phase with.

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
    A zero-phase Butterworth band-pass keeps periods between `pmin` and `pmax`, given as
    `band = (pmin, pmax)`, for frames `dt` apart. Both are required.

!!! tierfull
    **Units.** No time unit is assumed. `dt`, `band`, `trim` and `minlength` only have to be in
    the same unit, whatever it is: minutes from PhysiCell, hours from a microscope, or plain
    frames with `dt = 1`. Everything derived inherits that unit: cycle lengths are in it, and
    speeds are position units per it.

!!! tierfull
    **The conversion, worked through.** With `dt = 10` and `band = (100, 300)`, writing `u` for
    the time unit:

    | Quantity | Value |
    |:--|:--|
    | `dt` | 10 u per frame |
    | `band = (pmin, pmax)` | (100, 300) u, periods |
    | sampling rate `fs = 1/dt` | 0.1 samples per u |
    | low cutoff `1/pmax` | 1/300 cycles per u |
    | high cutoff `1/pmin` | 1/100 cycles per u |

    `bandpassFilter` makes exactly this conversion and passes `fs` to `DSP.digitalfilter`; it
    never hands DSP normalised frequencies. `pmin` must exceed two frames (`2dt`), the shortest
    period the frame interval can resolve, and every track must last at least two of the
    longest periods (`minlength ≥ 2pmax`).

    **What the cutoffs mean as applied.** `pmin` and `pmax` are the half-power points of one pass
    of the filter, where its amplitude gain is `1/√2`. The filter runs forward and then backward,
    which squares the gain, so in the filtered signal a component with period `pmin` or `pmax`
    keeps half its amplitude.

!!! tierfull
    The phase in step 3 is only meaningful for one dominant rhythm. Fast bursts would make the
    angle jump, and slow drift would move the analytic signal off the origin. The filter removes
    both: with order 3 and `band = (100, 300)`, a component of period 60 and a drift of period
    500 are both cut to about 1% of their amplitude, while periods 150–250 pass at more than 90%.

    It runs forward and then backward (`filtfilt`). The two passes cancel the filter's delay, so
    peaks stay where they were, and they square the gain: the cutoffs, the half-power points of
    one pass, become half-amplitude points.

    **What a Butterworth filter is.** Each filtered value is a weighted sum of recent raw values
    minus feedback from recent filtered values,
    `y[n] = b₀x[n] + … + b₆x[n−6] − a₁y[n−1] − … − a₆y[n−6]`. The family (Butterworth: the
    flattest possible pass band, no ripple), the order and the cutoffs fix all the coefficients;
    a band-pass doubles the order, to 6 for order 3. The `b`s block constants outright, which is
    why each agent's baseline disappears. The `a`s set three pole pairs; for `band = (100, 300)`
    and `dt = 10` they resonate at periods of about 108, 211 and 295 and together make the flat
    top.

    **Why order 3, applied once.** Steeper filters, or repeated passes, narrow the band, ring for
    longer, and impose a clean rhythm on irregular agents, hiding the cycle-to-cycle variability
    that is often the point.

The amplitude gain as applied (squared, for the forward and backward passes) at a few periods:

```@example phases
using DSP
dt = 10
f = bandpassFilter((100, 300); dt)
periods = [60, 100, 170, 300, 500]                           # same unit as dt
gain = abs2.(freqresp(f, 2π .* dt ./ periods))               # radians per sample = 2π·dt/period
round.(gain; digits = 3)
```

!!! tierdev
    `bandpassFilter` always passes `fs = 1/dt` to `DSP.digitalfilter`, so the cutoffs are in
    cycles per time unit of `dt`. Without `fs`, DSP reads them as fractions of the Nyquist frequency. The
    `filter` test pins the half-power points at the band edges, which is the check that catches
    a missing `fs`; a pure sine does not (see [Pitfalls](pitfalls.md)).

!!! tierjournal "2026-10-15 — One time unit, the user's; fs always passed; band required"
    Decided: the band is given as periods `(pmin, pmax)`, together with `trim`, `minlength` and
    the frame interval `dt`, all in one time unit the user chooses; no unit is named or assumed.
    Frames appear only in `frames` selections and in results. `bandpassFilter` converts to cycles
    per time unit and always passes `fs = 1/dt`.
    `band` has no default (it is a time scale of the system studied); `order = 3` does. The
    filter test measures the realised response at `pmin`, `pmax` and mid-band for four frame
    intervals.
    Why: the prototype once built the filter without `fs`; DSP read the cutoffs as fractions of
    Nyquist and passed periods of about 200–600 minutes instead of 100–300. The spurious
    "upward wave" it produced survived every pure-sine check. The mistake was mixing units (normalised frequency against
    periods), not the choice of unit; requiring everything in the unit of `dt` leaves nothing to
    mistranslate.

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
    `trim` of every track are dropped. `trim` defaults to the filter's memory, the time for its
    slowest pole to decay to 1/e, rounded up to whole frames: 180 (18 frames) for
    `band = (100, 300)` and `dt = 10`.

    **Masking.** When the analytic signal `z = y + i·H[y]` spirals in close to the origin, the
    band-passed signal is weak compared with the noise, and its angle is ill-defined: a small
    fluctuation swings the phase by a large amount. That shows up as phase slips, phase running
    backward (`ω ≤ 0`), spuriously short or long cycle lengths, and noisy gradient fits.
    `mask = q` sets the phase to `NaN` in each track's frames whose amplitude is below that
    track's own `q`-quantile; `mask = 0.2` drops each agent's weakest fifth. It is off by
    default, because it has a cost: the results that remain are biased toward moments when
    agents oscillate strongly, which can overstate how regular and how synchronised the rhythm
    is. Compare with and without it before relying on it.

```@example phases
p.trim, p.minlength
```

!!! tierdev
    `hilbert` is `DSP.hilbert` (FFT, zero the negative frequencies, double the positive, inverse
    FFT). The trim default comes from `_filtermemory(filter)`, `−1/log(max |pole|)` frames.

!!! tierjournal "2026-10-15 — No transform by default"
    Decided: `transform = identity` by default; `log` is opt-in.
    Why: the prototype logged notch because notch is heavy-tailed, a property of that signal,
    not of the method. A log default would also make `cellPhases` fail on any value ≤ 0. Like
    the band, it is the user's modelling choice.

!!! tierjournal "2026-10-15 — Masking off by default"
    Decided: low-amplitude masking is opt-in (`mask = q`, per-track amplitude quantile).
    Why: it removes ill-defined phases but biases what remains toward strongly oscillating
    moments; that trade-off is the user's to make. Phase-defect (spiral core) counting was also
    considered and left out of the first version.

!!! tierjournal "2026-10-15 — trim and minlength derived from the band"
    Decided: `trim` defaults to the filter's memory and `minlength` to `2trim + 2·band[2]`, both
    in the unit of `dt`.
    Why: the prototype fixed 18 and 90 frames for 10-minute frames and a 100–300-minute band.
    Both are really properties of the filter: 18 frames was "about the filter's memory", and 90
    frames "five cycles plus trim". Deriving them keeps them right when the band or the frame
    interval changes. For the prototype's settings the derived values are 18 frames, as
    before, and 96 frames rather than 90.
