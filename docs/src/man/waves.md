# Phase waves

Around each agent, fit a local plane wave to the neighbours' phases: its gradient gives the
direction of travel, and with the agent's own phase rate, the speed.

```@setup waves
using CellPhases, Random, Statistics
field = agentField(; n = 400, L = 300.0, nframes = 200, step = 0.5, rng = Xoshiro(1))
v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 90,
              rng = Xoshiro(2))
p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300))
```

## Direction

!!! tierbrief
    `phaseWaves` fits `wrap(φⱼ − φᵢ) ≈ ∇φ · (xⱼ − xᵢ)` by least squares over the neighbours
    within `r` (at least `minneighbours`), for every agent with a phase in the selected frames.
    The wave travels along −∇φ.

!!! tierfull
    Wrapping maps each difference into (−π, π], so 350° against 10° counts as 20° apart. A phase
    point, say the peak, moves toward agents that have not reached it yet, which have lower φ:
    hence −∇φ. `direction` is that angle in radians, 0 along +x and π/2 along +y.

    Each fit has an R². Where no single plane explains the neighbours (a wave front, a defect,
    or no wave at all) R² is low. Summaries use the fits with R² ≥ 0.5, and the share of such
    fits (`cleanShare`) must be reported with them: in the motivating model only 2–24% of cells
    had clean fits, depending on the signal.

```@example waves
w = phaseWaves(p; r = 30, frames = 1:4:200)
good = w.r2 .>= 0.5
rad2deg(angle(mean(cis.(w.direction[good])))), cleanShare(w)
```

## Speed

!!! tierbrief
    `speed = ω / |∇φ|`, with `ω` the agent's own phase advance to the next frame, in radians per
    minute. Speeds are in position units per minute, relative to the agents.

!!! tierfull
    On a moving front the phase is constant. A point moving `c·dt` along the direction of travel
    sees its phase change by `ω·dt − |∇φ|·c·dt = 0`, so `c = ω/|∇φ|`. `∇φ` comes from one frame;
    `ω` follows the same agent to the next frame by ID. Frames where the phase runs backward
    (`ω ≤ 0`) have no defined wave and get `speed = NaN`; the fit is kept for its direction.

    Because ω follows the agent, the speed is relative to the tissue. The speed in the lab adds
    the agents' own motion along the wave.

```@example waves
60 * median(filter(!isnan, w.speed[good]))     # per hour; the true speed is 30
```

!!! tierjournal "2026-10-15 — Speeds relative to the agents"
    Decided: `ω` is each agent's own phase advance, so speeds are relative to the agents.
    Why: the question is how a signal travels through the tissue, not how far the tissue itself
    moves. The prototype also computed a speed with one global period for all cells; per-agent ω
    is the honest version and is the one kept.

## The neighbourhood radius

!!! tierbrief
    Waves shorter than about `2r` cannot be read: neighbours more than half a wavelength away wrap
    past ±π and the fit fails.

!!! tierfull
    The plane fit works only while the true phase differences inside the neighbourhood stay
    within ±π, that is, while the wavelength exceeds about twice the radius. Slow waves have
    short wavelengths (`λ = c·T`), so they need a small radius; a small radius has fewer
    neighbours and noisier fits. Run two radii: if direction and clean share agree, the radius is
    not the limit. The share of clean fits collapses where the limit is crossed, so it warns of
    the problem.

    Synthetic plane waves (period 180 minutes, cortex-like density, slowly moving agents)
    through the whole pipeline, as in the test suite:

    | True speed, direction | Wavelength | r = 20 | r = 30 |
    |:--|:--|:--|:--|
    | 10 µm/h, down | 30 µm | 10.3 µm/h, −90°, 21% clean | fails (1% clean) |
    | 15 µm/h, down | 45 µm | 15.0 µm/h, −90° | fails (11% clean) |
    | 25 µm/h, up | 75 µm | 25.0 µm/h, 90° | 25.0 µm/h, 90° |
    | 40 µm/h, right | 120 µm | 39.8 µm/h, 0° | 39.9 µm/h, 0° |

!!! tierdev
    The fit is `_planefit` in `src/waves.jl`; it has no intercept, since the centre agent's own
    difference is 0 by construction. The R² is computed against zero, not against the mean
    difference, for the same reason.
