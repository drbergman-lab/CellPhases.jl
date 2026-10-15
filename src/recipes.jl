# Plots.jl recipes (through RecipesBase), drawn from the plot-data functions in plotdata.jl.
# `plot(result)` gives the default view of each readout; a second Symbol argument picks another.

@recipe function f(w::PhaseWaves, kind::Symbol = :rose; bin = 20, minr2 = 0.5)
    kind === :rose || throw(ArgumentError("PhaseWaves plots: :rose; got $(repr(kind))"))
    rose = roseData(w; bin, minr2)
    L = rose.rlim * 1.42
    aspect_ratio := :equal
    framestyle := :none
    grid := false
    legend := false
    xlims := (-L, L)
    ylims := (-L, L)
    # Subplot attributes, set before the series so every series carries them.
    annotations := [[(r * cosd(rose.labelangle), r * sind(rose.labelangle), (label, 8, :left, :gray40))
                     for (r, label) in rose.rings];
                    [(x, y, (label, 10, h)) for (x, y, label, h, _) in rose.compass];
                    (-L, -L, ("$(round(Int, rose.clean))% of fits with R² ≥ $minr2", 8, :left, :bottom, :gray40))]
    for θ in (0, 90, 180, 270)
        @series begin
            seriestype := :path
            linecolor := :gray85
            [0, rose.rlim * cosd(θ)], [0, rose.rlim * sind(θ)]
        end
    end
    @series begin
        seriestype := :shape
        linecolor := :white
        fill_z := rose.speed
        seriescolor --> :blues
        colorbar_title --> "median speed"
        _nanjoin(rose.wedges)
    end
    for (r, _) in rose.rings
        @series begin
            seriestype := :path
            linecolor := :black
            linealpha := 0.3
            _xy(_circle(r))
        end
    end
    for (c, lw) in ((:white, 5), (:black, 2))
        @series begin
            seriestype := :path
            linecolor := c
            linewidth := lw
            arrow := true
            [0, rose.arrow[1]], [0, rose.arrow[2]]
        end
    end
    return nothing
end

@recipe function f(c::PhaseCoherence; skipzero = true)
    d = coherenceData(c; skipzero)
    xguide --> "distance"
    yguide --> "excess phase coherence"
    @series begin
        seriestype := :path
        fillrange := d.low
        fillalpha := 0.2
        linealpha := 0
        seriescolor := :gray
        label := "null 95%"
        d.x, d.high
    end
    @series begin
        seriestype := :path
        linewidth --> 2
        label --> "⟨cos Δφ⟩ − null"
        d.x, d.y
    end
end

@recipe function f(o::OrderParameter)
    d = orderParameterData(o)
    xguide --> "radius"
    yguide --> "order parameter R"
    ylims --> (0, 1)
    @series begin
        seriestype := :path
        fillrange := d.low
        fillalpha := 0.2
        linealpha := 0
        seriescolor := :gray
        label := "null 95%"
        d.radii, d.high
    end
    @series begin
        seriestype := :path
        linestyle := :dash
        seriescolor := :gray
        label := "time-shifted"
        d.radii, d.null
    end
    @series begin
        seriestype := :path
        markershape --> :circle
        linewidth --> 2
        label --> "actual"
        d.radii, d.R
    end
end

@recipe function f(c::CycleLengths; edges = nothing)
    d = cycleLengthData(c; edges)
    xguide --> "cycle length"
    yguide --> "cycles"
    legend --> false
    @series begin
        seriestype := :bar
        bar_width := diff(d.edges)
        (d.edges[1:end-1] .+ d.edges[2:end]) ./ 2, d.counts
    end
    @series begin
        seriestype := :vline
        linestyle := :dash
        linecolor := :black
        [d.median]
    end
end

@recipe function f(v::PhaseValueCheck)
    d = phaseValueData(v)
    xguide --> "phase φ"
    yguide --> "value, z-scored within each track"
    xlims --> (-π, π)
    legend --> false
    @series begin
        seriestype := :heatmap
        seriescolor --> :blues
        (d.phaseedges[1:end-1] .+ d.phaseedges[2:end]) ./ 2,
        (d.zedges[1:end-1] .+ d.zedges[2:end]) ./ 2, permutedims(d.density)
    end
    @series begin
        seriestype := :path
        fillrange := d.q25
        fillalpha := 0.15
        linealpha := 0
        seriescolor := :black
        d.centers, d.q75
    end
    @series begin
        seriestype := :path
        linewidth := 2
        linecolor := :black
        d.centers, d.median
    end
end

@recipe function f(p::PhaseResult, kind::Symbol, id)
    if kind === :trace
        d = traceData(p, id)
        layout := (2, 1)
        legend := false
        @series begin
            subplot := 1
            yguide := "value"
            d.time, d.value
        end
        @series begin
            subplot := 2
            seriestype := :scatter
            markersize := 2
            xguide := "time"
            yguide := "phase"
            ylims := (-π, π)
            d.time, d.phase
        end
    elseif kind === :portrait
        d = portraitData(p, id)
        aspect_ratio := :equal
        legend := false
        xguide --> "band-passed signal"
        yguide --> "Hilbert quadrature"
        @series begin
            line_z := d.time
            seriescolor --> :blues
            d.re, d.im
        end
    else
        throw(ArgumentError("PhaseResult plots: :trace, :portrait; got $(repr(kind))"))
    end
end

_xy(pts) = (first.(pts), last.(pts))

function _nanjoin(polys)
    x, y = Float64[], Float64[]
    for (k, poly) in enumerate(polys)
        k > 1 && (push!(x, NaN); push!(y, NaN))
        append!(x, first.(poly))
        append!(y, last.(poly))
    end
    return x, y
end
