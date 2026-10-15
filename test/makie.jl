# The Makie extension, with Makie alone (no backend): each recipe builds its plot and sets its
# preferred axis. Rendering is checked in the docs build with CairoMakie.

import Makie

@testset "Makie extension" begin
    v = planeWave(FIELD.positions; dt = DT, period = 180, speed = 0.5, direction = 90, rng = Xoshiro(31))
    p = phases(v)
    w = phaseWaves(p; r = 30, frames = 1:8:200)
    c = phaseCoherence(p; dr = 10, rmax = 40, frames = 1:8:200, nsurrogates = 3)
    o = orderParameter(p; radii = [10, 20], frames = 1:8:200, nsurrogates = 3)

    fig = Makie.Figure()
    ax, rose = roseplot(fig[1, 1], w)
    @test ax.aspect[] isa Makie.DataAspect
    @test !ax.xticklabelsvisible[] && !ax.bottomspinevisible[]
    @test length(rose.plots) == 9           # frame, axes, wedges, rings, ring labels, 2 arrows, compass, note
    @test rose.wedgespeed[] ≈ replace(roseData(w).speed, NaN => 0.0)
    @test eltype(rose.ringlines[]) == Makie.Point2f
    ax2, _ = coherenceplot(fig[1, 2], c)
    @test ax2.xlabel[] == "distance"
    orderparameterplot(fig[1, 3], o)
    cyclelengthplot(fig[2, 1], cycleLengths(p))
    phasevalueplot(fig[2, 2], phaseValueCheck(p))
    ax6, port = portraitplot(fig[2, 3], p, 7)
    @test ax6.aspect[] isa Makie.DataAspect
    @test length(port.trace[]) == length(portraitData(p, 7).re)
    ax7 = Makie.Axis(fig[3, 1])
    @test roseplot!(ax7, w) isa Makie.Plot
    # plot(result) picks the recipe
    @test Makie.plot(fig[3, 2], w).plot isa Makie.Plot{roseplot}
    @test Makie.plot(fig[3, 3], c).plot isa Makie.Plot{coherenceplot}
end
