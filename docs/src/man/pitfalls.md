# Pitfalls

Mistakes that were made once while developing this analysis, and what guards against them.

!!! tierbrief
    Each item names the mistake, what it did, and the guard.

**The filter's sampling rate.**

!!! tierfull
    Without the sampling rate, DSP.jl reads cutoffs as fractions of the Nyquist frequency. A
    prototype filter passed 200–600 minutes instead of 100–300, suppressed the real rhythm, and
    produced a spurious "upward wave". Pure-sine tests missed it, because a single frequency
    keeps its phase when attenuated. `bandpassFilter` takes the band in minutes and always passes
    the sampling rate; the tests pin the half-power points and use a pulse train.

**Phase wrapping at large radius.**

!!! tierfull
    Neighbourhoods wider than about half a wavelength wrap phase differences past ±π, and the
    gradient fit fails. Slow waves have short wavelengths. Compare two radii, and watch the share
    of clean fits. See [Phase waves](waves.md#The-neighbourhood-radius).

**Pooling hides directions.**

!!! tierfull
    A pooled space-time correlation map averages opposite local waves to zero drift. Local phase
    gradients see them; a pooled map does not.

**Raw bursty signals.**

!!! tierfull
    A bursty signal decorrelates within tens of minutes through bursts that do not travel. Analyse
    the band-passed oscillation for waves, not the raw value.

**Selection by fit quality.**

!!! tierfull
    Directions and speeds come only from clean fits, which can be a small minority of agents (2–24%
    in the motivating model). Report `cleanShare` next to any direction or speed.

**Movie colour scales.**

!!! tierfull
    A colour scale that rescales every frame makes each frame's strongest burst look equally
    strong and exaggerates apparent turbulence. Fix the colour range across frames.
