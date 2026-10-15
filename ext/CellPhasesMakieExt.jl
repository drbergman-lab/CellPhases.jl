# Makie recipes, loaded with Makie. Each draws from a plot-data function in src/plotdata.jl and
# sets its preferred axis through `Makie.preferred_axis_attributes` (equal aspect and no
# decorations for the rose). `plot(result)` maps each result type to its recipe.

module CellPhasesMakieExt

using CellPhases
using CellPhases: PhaseWaves, PhaseCoherence, OrderParameter, CycleLengths, PhaseValueCheck,
                  PhaseResult
import CellPhases: roseplot, roseplot!, coherenceplot, coherenceplot!, orderparameterplot,
                   orderparameterplot!, cyclelengthplot, cyclelengthplot!, phasevalueplot,
                   phasevalueplot!, portraitplot, portraitplot!
using Makie: Makie, @recipe, Point2f, Vec2f, poly!, lines!, scatter!, text!, arrows2d!, band!, scatterlines!,
             barplot!, vlines!, hlines!, heatmap!, DataAspect, Axis

const INK = Makie.RGBf(0.12, 0.12, 0.12)
const MUTED = Makie.RGBf(0.42, 0.42, 0.40)
const SURFACE = Makie.RGBf(0.99, 0.99, 0.98)

_points(pts) = [Point2f(x, y) for (x, y) in pts]

#==================================================================================================#
# Rose

@recipe RosePlot (waves,) begin
    "Wedge width in degrees; must divide 360."
    bin = 20
    "Only fits with at least this R² are drawn."
    minr2 = 0.5
    "Shares (%) at which reference rings are drawn."
    rings = (5, 10)
    "Colormap for the median speed of each wedge."
    colormap = :Blues
    "Colour of the background the halo and wedge outlines blend into."
    surfacecolor = SURFACE
    "Font size of the compass labels; other labels are smaller."
    fontsize = 12
end

Makie.preferred_axis_attributes(::Type{<:Axis}, ::RosePlot) =
    (aspect = DataAspect(), xgridvisible = false, ygridvisible = false,
     leftspinevisible = false, rightspinevisible = false, topspinevisible = false,
     bottomspinevisible = false, xticksvisible = false, yticksvisible = false,
     xticklabelsvisible = false, yticklabelsvisible = false)

function Makie.plot!(p::RosePlot)
    map!(p.attributes, [:waves, :bin, :minr2, :rings], :rose) do w, bin, minr2, rings
        roseData(w; bin, minr2, rings)
    end
    map!(p.attributes, :rose, [:wedgepolys, :wedgespeed, :axes, :ringlines, :ringpos, :ringlabels,
                               :arrowtip, :compasspos, :compasslabels, :compassalign, :notepos,
                               :note, :frame]) do r
        θ = r.labelangle
        L = r.rlim * 1.42
        speed = [isnan(s) ? 0.0 : s for s in r.speed]
        axes = Point2f[]
        for t in (0, 90, 180, 270)
            push!(axes, Point2f(0, 0), Point2f(r.rlim * cosd(t), r.rlim * sind(t)), Point2f(NaN))
        end
        # NaN-separated; push! rather than `[pts; Point2f(NaN)]`, which would splice the point's
        # coordinates in as separate elements.
        rings = Point2f[]
        for (ρ, _) in r.rings
            append!(rings, _points(CellPhases._circle(ρ)))
            push!(rings, Point2f(NaN))
        end
        return ([_points(w) for w in r.wedges], speed, axes, rings,
                [Point2f(ρ * cosd(θ), ρ * sind(θ)) for (ρ, _) in r.rings], [l for (_, l) in r.rings],
                Vec2f(r.arrow...), [Point2f(x, y) for (x, y, _, _, _) in r.compass],
                [l for (_, _, l, _, _) in r.compass], [(h, v) for (_, _, _, h, v) in r.compass],
                Point2f(-L, -L), "$(round(Int, r.clean))% of fits with a clean gradient",
                [Point2f(-L, -L), Point2f(L, L)])
    end
    # Text does not count toward axis limits; transparent corners keep the labels inside.
    scatter!(p, p.frame; color = :transparent, markersize = 0)
    lines!(p, p.axes; color = (MUTED, 0.3), linewidth = 1)
    poly!(p, p.wedgepolys; color = p.wedgespeed, colormap = p.colormap, strokecolor = p.surfacecolor,
          strokewidth = 1.5)
    lines!(p, p.ringlines; color = (INK, 0.3), linewidth = 0.8)
    text!(p, p.ringpos; text = p.ringlabels, color = MUTED, fontsize = 10, align = (:left, :top),
          offset = (3, -2))
    for (c, sw, tl, tw) in ((p.surfacecolor, 6, 14, 15), (INK, 2.5, 10, 10))
        arrows2d!(p, [Point2f(0, 0)], map(v -> [v], p.arrowtip); color = c, shaftwidth = sw,
                  tiplength = tl, tipwidth = tw)
    end
    text!(p, p.compasspos; text = p.compasslabels, align = p.compassalign, color = INK,
          fontsize = p.fontsize)
    text!(p, p.notepos; text = p.note, color = MUTED, fontsize = 10, align = (:left, :bottom))
    return p
end

#==================================================================================================#
# Coherence and order parameter

@recipe CoherencePlot (coherence,) begin
    "Drop the first distance bin, which holds only the closest pairs."
    skipzero = true
    color = Makie.wong_colors()[1]
    bandcolor = (MUTED, 0.2)
    linewidth = 2
end

Makie.preferred_axis_attributes(::Type{<:Axis}, ::CoherencePlot) =
    (xlabel = "distance", ylabel = "excess phase coherence")

function Makie.plot!(p::CoherencePlot)
    map!(p.attributes, [:coherence, :skipzero], [:x, :y, :low, :high]) do c, skipzero
        d = coherenceData(c; skipzero)
        return d.x, d.y, d.low, d.high
    end
    band!(p, p.x, p.low, p.high; color = p.bandcolor, label = "null 95%")
    hlines!(p, [0]; color = MUTED, linewidth = 1)
    lines!(p, p.x, p.y; color = p.color, linewidth = p.linewidth, label = "⟨cos Δφ⟩ − null")
    return p
end

@recipe OrderParameterPlot (orderparameter,) begin
    color = Makie.wong_colors()[1]
    bandcolor = (MUTED, 0.2)
    linewidth = 2
end

Makie.preferred_axis_attributes(::Type{<:Axis}, ::OrderParameterPlot) =
    (xlabel = "radius", ylabel = "order parameter R", limits = (nothing, (0, 1)))

function Makie.plot!(p::OrderParameterPlot)
    map!(p.attributes, :orderparameter, [:radii, :R, :null, :low, :high]) do o
        d = orderParameterData(o)
        return d.radii, d.R, d.null, d.low, d.high
    end
    band!(p, p.radii, p.low, p.high; color = p.bandcolor, label = "null 95%")
    lines!(p, p.radii, p.null; color = MUTED, linestyle = :dash, linewidth = p.linewidth,
           label = "time-shifted")
    scatterlines!(p, p.radii, p.R; color = p.color, linewidth = p.linewidth, label = "actual")
    return p
end

#==================================================================================================#
# Cycle lengths, phase against value, portrait

@recipe CycleLengthPlot (cyclelengths,) begin
    "Histogram bin edges; 30 bins over the central 98% by default."
    edges = nothing
    color = Makie.wong_colors()[1]
end

Makie.preferred_axis_attributes(::Type{<:Axis}, ::CycleLengthPlot) =
    (xlabel = "cycle length", ylabel = "cycles")

function Makie.plot!(p::CycleLengthPlot)
    map!(p.attributes, [:cyclelengths, :edges], [:centres, :counts, :widths, :median]) do c, edges
        d = cycleLengthData(c; edges)
        return (d.edges[1:end-1] .+ d.edges[2:end]) ./ 2, d.counts, diff(d.edges), [d.median]
    end
    barplot!(p, p.centres, p.counts; width = p.widths, gap = 0, color = p.color,
             strokecolor = SURFACE, strokewidth = 1)
    vlines!(p, p.median; color = INK, linestyle = :dash, linewidth = 1.5)
    return p
end

@recipe PhaseValuePlot (check,) begin
    colormap = :Blues
end

Makie.preferred_axis_attributes(::Type{<:Axis}, ::PhaseValuePlot) =
    (xlabel = "phase φ", ylabel = "value, z-scored within each track",
     xticks = ([-π, -π / 2, 0, π / 2, π], ["−π", "−π/2", "0", "π/2", "π"]))

function Makie.plot!(p::PhaseValuePlot)
    map!(p.attributes, :check, [:φedges, :zscoreedges, :density, :centers, :median, :q25, :q75]) do v
        d = phaseValueData(v)
        return d.phaseedges, d.zscoreedges, d.density, d.centers, d.median, d.q25, d.q75
    end
    heatmap!(p, p.φedges, p.zscoreedges, p.density; colormap = p.colormap)
    band!(p, p.centers, p.q25, p.q75; color = (INK, 0.12))
    lines!(p, p.centers, p.median; color = INK, linewidth = 2.5)
    return p
end

@recipe PortraitPlot (phases, id) begin
    colormap = :Blues
    linewidth = 1.5
end

Makie.preferred_axis_attributes(::Type{<:Axis}, ::PortraitPlot) =
    (aspect = DataAspect(), xlabel = "band-passed signal", ylabel = "Hilbert quadrature")

function Makie.plot!(p::PortraitPlot)
    map!(p.attributes, [:phases, :id], [:trace, :time]) do ph, id
        d = portraitData(ph, id)
        return Point2f.(d.re, d.im), d.time
    end
    lines!(p, p.trace; color = p.time, colormap = p.colormap, linewidth = p.linewidth)
    return p
end

#==================================================================================================#
# plot(result) picks the recipe

Makie.plottype(::PhaseWaves) = RosePlot
Makie.plottype(::PhaseCoherence) = CoherencePlot
Makie.plottype(::OrderParameter) = OrderParameterPlot
Makie.plottype(::CycleLengths) = CycleLengthPlot
Makie.plottype(::PhaseValueCheck) = PhaseValuePlot

end
