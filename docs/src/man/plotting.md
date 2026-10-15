# Plotting

Each readout has a function that returns plot-ready numbers, and a recipe that draws them.

!!! tierbrief
    `roseData`, `coherenceData`, `orderParameterData`, `cycleLengthData`, `phaseValueData`,
    `traceData` and `portraitData` return plain numbers and polygons, for any plotting package.
    With Plots.jl loaded, `plot(result)` draws the default view.

!!! tierfull
    The plot-data functions hold every decision about what to draw, so a hand-built figure and a
    recipe show the same thing. The rose diagram is the main one: wedges of `bin` degrees whose
    length is `√share`, so that wedge *area* is proportional to the number of fits; reference
    rings at 5% and 10%, labelled along the shortest wedge where they overlap least; and an arrow
    for the mean direction whose length is the mean resultant length `R` times the outer radius.
    Wedges are coloured by median speed.

```@setup plotting
using Plots
gr(fmt = :png)
```

```@example plotting
using CellPhases, Random, Plots
field = agentField(; n = 400, L = 300.0, nframes = 200, step = 0.5, rng = Xoshiro(1))
v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 60, noise = 0.6,
              rng = Xoshiro(2))
p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300))
w = phaseWaves(p; r = 30, frames = 1:4:200)
c = phaseCoherence(p; dr = 10, rmax = 80, frames = 1:4:200, rng = Xoshiro(3))
o = orderParameter(p; radii = [10, 20, 30, 40, 60], frames = 1:4:200, rng = Xoshiro(3))

plot(plot(w), plot(c), plot(o), plot(cycleLengths(p)), plot(phaseValueCheck(p)),
     plot(p, :portrait, 7); layout = (2, 3), size = (1200, 750))
```

One agent's trace and phase:

```@example plotting
plot(p, :trace, 7; size = (800, 400))
```

The data behind the rose:

```@example plotting
rose = roseData(w)
(; rose.meandirection, rose.resultant, rose.clean, rose.up, rose.down)
```

!!! tierdev
    Recipes are in `src/recipes.jl` and use only RecipesBase; each draws from a plot-data
    function and adds nothing of its own beyond styling. A second `Symbol` argument picks among
    views (`plot(p, :trace, id)`, `plot(p, :portrait, id)`, `plot(w, :rose)`). Tests run the
    recipes through `RecipesBase.apply_recipe` without a plotting package.

!!! tierjournal "2026-10-15 — Plot data in core, recipes from it, no packaged composite"
    Decided: the core has plot-data functions and RecipesBase recipes built on them; no
    multi-panel figure is packaged.
    Why: plot data is plain numbers, so it carries no plotting dependency and serves any
    package. A composite figure is an analysis choice (which signals, which panels, which cross-
    check with SpaceTimeCorrelations) and stays in analysis scripts.
