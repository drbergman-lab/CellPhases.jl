# Makie recipe functions. The methods live in ext/CellPhasesMakieExt.jl and load with Makie (any
# backend: `using CairoMakie`, `using GLMakie`, …). `plot(result)` in Makie also works and picks
# the matching recipe.

"""
    roseplot(w::PhaseWaves; bin = 20, minr2 = 0.5, rings = (5, 10), colormap = :Blues)
    roseplot!(ax, w; …)

Makie recipe: the directions of local phase waves as a rose diagram, drawn from
[`roseData`](@ref). Wedges are filled and coloured by median speed, on an equal-aspect axis with
hidden decorations; compass labels right/up/left/down, reference rings labelled along the shortest
wedge, the mean-direction arrow with a light halo, and the clean-fit share. Requires a Makie
backend.

# Examples
```julia
using CairoMakie, CellPhases
roseplot(PhaseWaves(p; r = 30))
```
"""
function roseplot end

"""
    roseplot!(ax, w::PhaseWaves; kwargs...)

In-place form of [`roseplot`](@ref).
"""
function roseplot! end

"""
    coherenceplot(c::PhaseCoherence; skipzero = true, color, bandcolor)
    coherenceplot!(ax, c; …)

Makie recipe: excess phase coherence against distance with the surrogate 95% band, from
[`coherenceData`](@ref). Requires a Makie backend.
"""
function coherenceplot end

"""
    coherenceplot!(ax, c::PhaseCoherence; kwargs...)

In-place form of [`coherenceplot`](@ref).
"""
function coherenceplot! end

"""
    orderparameterplot(o::OrderParameter; color, bandcolor)
    orderparameterplot!(ax, o; …)

Makie recipe: the local order parameter against radius with its surrogate null and 95% band,
from [`orderParameterData`](@ref). Requires a Makie backend.
"""
function orderparameterplot end

"""
    orderparameterplot!(ax, o::OrderParameter; kwargs...)

In-place form of [`orderparameterplot`](@ref).
"""
function orderparameterplot! end

"""
    cyclelengthplot(c::CycleLengths; edges = nothing, color)
    cyclelengthplot!(ax, c; …)

Makie recipe: a histogram of cycle lengths with the median marked, from
[`cycleLengthData`](@ref). Requires a Makie backend.
"""
function cyclelengthplot end

"""
    cyclelengthplot!(ax, c::CycleLengths; kwargs...)

In-place form of [`cyclelengthplot`](@ref).
"""
function cyclelengthplot! end

"""
    phasevalueplot(v::PhaseValueCheck; colormap = :Blues)
    phasevalueplot!(ax, v; …)

Makie recipe: the z-scored value against phase as a density, with the binned median and
quartiles, from [`phaseValueData`](@ref). Requires a Makie backend.
"""
function phasevalueplot end

"""
    phasevalueplot!(ax, v::PhaseValueCheck; kwargs...)

In-place form of [`phasevalueplot`](@ref).
"""
function phasevalueplot! end

"""
    portraitplot(p::PhaseResult, id; colormap = :Blues)
    portraitplot!(ax, p, id; …)

Makie recipe: one agent's analytic signal in the complex plane, coloured by time, from
[`portraitData`](@ref). Requires a Makie backend.
"""
function portraitplot end

"""
    portraitplot!(ax, p::PhaseResult, id; kwargs...)

In-place form of [`portraitplot`](@ref).
"""
function portraitplot! end
