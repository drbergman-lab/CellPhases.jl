using CellPhases
using Test
using Random, Statistics
using RecipesBase: RecipesBase
using DSP: freqresp

# Real-like geometry: a cortex-like density (~15 µm spacing), agents drifting slowly, 10-min frames.
const FIELD = agentField(; n = 400, L = 300.0, nframes = 200, step = 0.5, rng = Xoshiro(1))
const DT = 10.0
const BAND = (100, 300)
phases(v; kw...) = cellPhases(FIELD.positions, v, FIELD.ids; dt = DT, band = BAND, kw...)
circmean(θ) = angle(mean(cis.(θ)))
angdiff(a, b) = abs(angle(cis(a - b)))

@testset "CellPhases" begin

@testset "filter" begin
    f = bandpassFilter(BAND; dt = DT)
    # Half-power points of one pass sit at the band edges, given in minutes: this is what goes
    # wrong when the sampling rate is not passed to DSP.
    H = abs.(freqresp(f, 2π * DT .* [1 / 100, 1 / 300, 1 / 170, 1 / 50, 1 / 1000]))
    @test H[1] ≈ 1 / √2 atol = 1e-3
    @test H[2] ≈ 1 / √2 atol = 1e-3
    @test H[3] > 0.99
    @test H[4] < 0.1 && H[5] < 0.1
    @test_throws ArgumentError bandpassFilter((300, 100); dt = DT)
    @test_throws ArgumentError bandpassFilter((15, 300); dt = DT)      # shorter than 2 frames
    @test_throws ArgumentError bandpassFilter(BAND; dt = DT, order = 0)
    @test 17 < CellPhases._filtermemory(f) < 18
end

@testset "tracks, trim and minlength" begin
    v = planeWave(FIELD.positions; dt = DT, period = 180, speed = 0.5, direction = 0, rng = Xoshiro(2))
    p = phases(v)
    @test p.trim == 180.0                     # filter memory, rounded up to whole frames
    @test p.minlength == 2 * 180 + 2 * 300
    @test length(p.tracks) == 400
    @test all(tr -> tr.frames == 19:182, p.tracks)
    @test all(t -> all(isnan, p.phase[t]), [1:18; 183:200])
    @test !any(isnan, p.phase[100])
    @test all(x -> -π <= x <= π, p.phase[100])

    # an agent missing from frames 90–91 is split in two; both pieces are too short
    ids = deepcopy(FIELD.ids)
    pos = deepcopy(FIELD.positions)
    vv = deepcopy(v)
    for t in 90:91
        deleteat!(ids[t], 1); deleteat!(pos[t], 1); deleteat!(vv[t], 1)
    end
    q = cellPhases(pos, vv, ids; dt = DT, band = BAND)
    @test [tr.frames for tr in q.tracks if tr.id == 1] == [110:182]   # frames 1–89 are too short
    q = cellPhases(pos, vv, ids; dt = DT, band = BAND, minlength = 600, trim = 100)
    one = filter(tr -> tr.id == 1, q.tracks)
    @test [tr.frames for tr in one] == [11:79, 102:190]
    @test q.trim == 100.0

    # explicit settings, errors
    @test_throws ArgumentError cellPhases(pos, vv, [fill(1, length(i)) for i in ids]; dt = DT, band = BAND)
    @test_throws DimensionMismatch cellPhases(pos[1:3], vv, ids; dt = DT, band = BAND)
    neg = deepcopy(v); neg[5][3] = -1.0
    @test_throws ArgumentError phases(neg)
    @test phases(neg; transform = identity) isa PhaseResult
    @test_throws ArgumentError phases(v; mask = 1.0)

    # IDs of any type, 3D positions
    p3 = cellPhases([[(x[1], x[2], 0.0) for x in f] for f in FIELD.positions], v,
                    [string.(i) for i in FIELD.ids]; dt = DT, band = BAND)
    @test p3 isa PhaseResult{3,String}
    @test isequal(p3.phase, p.phase)
end

@testset "sign conventions" begin
    # noise-free wave: the phase is ωt − k·x, so φ = 0 at the peaks
    T, c = 180.0, 0.5
    still = agentField(; n = 200, L = 200.0, nframes = 200, rng = Xoshiro(2))
    v = planeWave(still.positions; dt = DT, period = T, speed = c, direction = 90, noise = 0)
    p = cellPhases(still.positions, v, still.ids; dt = DT, band = BAND)
    ω, k = 2π / T, 2π / T / c
    t = 100
    truth = [ω * (t - 1) * DT - k * x[2] for x in still.positions[t]]
    @test maximum(angdiff.(p.phase[t], truth)) < 0.05
    # a wave moving +y reads +90°
    w = phaseWaves(p; r = 30, frames = 50:10:150)
    @test rad2deg(circmean(w.direction)) ≈ 90 atol = 1
end

@testset "plane waves: speed and direction" begin
    # The prototype's validation table (µm/h), through real-like geometry.
    for (speed, θ) in ((10, -90), (15, -90), (25, 90), (40, 0)), r in (20.0, 30.0)
        T = 180.0
        λ = speed / 60 * T
        v = planeWave(FIELD.positions; dt = DT, period = T, speed = speed / 60, direction = θ,
                      rng = Xoshiro(3))
        w = phaseWaves(phases(v); r, frames = 1:4:200)
        good = w.r2 .>= 0.5
        if λ > 2r
            @test cleanShare(w) > 0.8
            @test 60median(filter(!isnan, w.speed[good])) ≈ speed rtol = 0.05
            @test rad2deg(angdiff(circmean(w.direction[good]), deg2rad(θ))) < 3
        else
            # The wrapping limit: below a wavelength of ~2r, at most half the fits come out clean.
            @test cleanShare(w) < 0.5
        end
    end
end

@testset "broadband signal: cycle length" begin
    # A pulse train puts most of its power in harmonics; only a correct band recovers its period.
    v = pulseTrain(FIELD.ids; dt = DT, period = 170, rng = Xoshiro(4))
    p = phases(v)
    c = cycleLengths(p)
    @test median(c.lengths) ≈ 170 atol = 3
    @test length(c.lengths) > 2000
    @test eltype(c.ids) == Int
    # the same pulses read through a band that excludes the fundamental give the 2nd harmonic
    q = cellPhases(FIELD.positions, v, FIELD.ids; dt = DT, band = (60, 100))
    @test median(cycleLengths(q).lengths) ≈ 85 atol = 3
    # a two-frequency signal: a weak harmonic does not move the cycle length
    two = [[exp(1 + cos(2π * (t - 1) * DT / 170) + 0.6cos(2π * (t - 1) * DT / 85 + 1)) for _ in f]
           for (t, f) in enumerate(FIELD.ids)]
    @test median(cycleLengths(phases(two)).lengths) ≈ 170 atol = 1
end

@testset "phase tracks the value" begin
    v = planeWave(FIELD.positions; dt = DT, period = 180, speed = 0.5, direction = 0, rng = Xoshiro(5))
    c = phaseValueCheck(phases(v))
    @test c.r2 > 0.8                     # noise SD 0.3 against a unit cosine
    @test abs(rad2deg(c.peak)) < 5
    @test length(c.centers) == 36
    @test c.median[18] > 1 && c.median[1] < -1          # peak at φ ≈ 0, trough at ±π
end

@testset "time-shift surrogates" begin
    v = planeWave(FIELD.positions; dt = DT, period = 180, speed = 0.5, direction = 0, rng = Xoshiro(6))
    p = phases(v)
    s = timeShift(p; rng = Xoshiro(7))
    for tr in p.tracks[1:20]
        @test sort([p.phase[t][r] for (t, r) in zip(tr.frames, tr.rows)]) ==
              sort([s.phase[t][r] for (t, r) in zip(tr.frames, tr.rows)])
    end
    @test isequal(isnan.(s.phase[100]), isnan.(p.phase[100]))
    @test s.phase[100] != p.phase[100]
end

@testset "nulls on uncoupled oscillators" begin
    v = oscillatorField(FIELD.positions, FIELD.ids; dt = DT, period = 180, rng = Xoshiro(8))
    p = phases(v)
    c = phaseCoherence(p; dr = 10, rmax = 60, frames = 1:4:200, rng = Xoshiro(9))
    @test all(abs.(c.coherence[2:end]) .< 0.05)
    @test all(abs.(c.excess[2:end]) .< 0.05)
    @test all(0 .< c.nullhigh[2:end] .- c.nulllow[2:end] .< 0.15)
    @test all(c.pvalue[2:end] .> 0.01)
    @test size(c.surrogates) == (7, 19)
    o = orderParameter(p; radii = [15, 30], frames = 1:4:200, rng = Xoshiro(9))
    @test all(abs.(o.excess) .< 0.03)
    @test all(o.R .> 0.1)                       # ≈ 1/√n, not 0
end

@testset "coherence on coupled oscillators" begin
    v = oscillatorField(FIELD.positions, FIELD.ids; dt = DT, period = 180,
                        correlationlength = 25.0, phasesd = 1.0, rng = Xoshiro(10))
    p = phases(v)
    c = phaseCoherence(p; dr = 10, rmax = 80, frames = 1:4:200, rng = Xoshiro(11))
    @test issorted(c.coherence[2:end]; rev = true)
    @test all(c.coherence[2:5] .> c.nullhigh[2:5])
    @test all(c.pvalue[2:5] .== 1 / 20)
    o = orderParameter(p; radii = [10, 20, 40], frames = 1:4:200, rng = Xoshiro(11))
    @test all(o.R .> o.nullhigh)
    @test issorted(o.R; rev = true)
end

@testset "masking" begin
    v = pulseTrain(FIELD.ids; dt = DT, period = 170, rng = Xoshiro(12))
    p = phases(v; mask = 0.2)
    n(p) = sum(φ -> count(!isnan, φ), p.phase)
    @test n(p) / n(phases(v)) ≈ 0.8 atol = 0.01
    @test !any(isnan, p.amplitude[100])
end

@testset "replicates" begin
    runs = [cellPhases(f.positions, planeWave(f.positions; dt = DT, period = 180, speed = 0.5,
                                              direction = 0, rng = Xoshiro(k)), f.ids; dt = DT, band = BAND)
            for (k, f) in enumerate(agentField(; n = 150, L = 150.0, nframes = 150, rng = Xoshiro(k))
                                    for k in 1:3)]
    c = phaseCoherence(runs; dr = 10, rmax = 40, frames = 1:4:150, nsurrogates = 5, rng = Xoshiro(1))
    @test c isa Replicates{PhaseCoherence}
    @test length(c.results) == 3
    @test c.pooled.counts == sum(r.counts for r in c.results)
    @test c.pooled.coherence ≈ sum(r.coherence .* r.counts for r in c.results) ./ c.pooled.counts
    spread = replicateSpread(c, :excess)
    @test spread.n == fill(3, 5)
    @test all(spread.std[2:end] .< 0.1)
    cl = cycleLengths(runs)
    @test length(cl.pooled.lengths) == sum(length(r.lengths) for r in cl.results)
    w = phaseWaves(runs; r = 30, frames = 1:10:150)
    @test length(w.pooled.r2) == sum(length(r.r2) for r in w.results)
    @test phaseValueCheck(runs).pooled.r2 > 0.8
    @test orderParameter(runs; radii = [20], nsurrogates = 2, frames = 1:10:150).pooled isa OrderParameter
end

@testset "plot data and recipes" begin
    v = planeWave(FIELD.positions; dt = DT, period = 180, speed = 0.5, direction = 90, rng = Xoshiro(13))
    p = phases(v)
    w = phaseWaves(p; r = 30, frames = 1:4:200)
    rose = roseData(w)
    @test sum(rose.share) ≈ 100
    @test length(rose.wedges) == 18
    @test rose.lo[argmax(rose.share)] == 80
    @test rose.meandirection ≈ 90 atol = 2
    @test rose.resultant > 0.95
    @test 60 * rose.speed[argmax(rose.share)] ≈ 30 rtol = 0.05
    @test rose.rings[1] == (sqrt(5.0), "5%")
    @test_throws ArgumentError roseData(w; bin = 7)

    c = phaseCoherence(p; dr = 10, rmax = 40, frames = 1:8:200, nsurrogates = 3)
    @test length(coherenceData(c).x) == 4
    o = orderParameter(p; radii = [10, 20], frames = 1:8:200, nsurrogates = 3)
    @test orderParameterData(o).R == o.R
    cl = cycleLengthData(cycleLengths(p))
    @test sum(cl.counts) <= length(cycleLengths(p).lengths)
    pv = phaseValueData(phaseValueCheck(p))
    @test all(sum(pv.density; dims = 2) .<= 1 + 1e-12)
    tr = traceData(p, 7)
    @test length(tr.time) == 200 && tr.time[2] == DT
    @test_throws ArgumentError traceData(p, -1)
    @test length(portraitData(p, 7).re) == count(!isnan, tr.phase)

    # recipes, without a plotting package (Plots defines is_key_supported for its backends)
    RecipesBase.is_key_supported(::Symbol) = true
    series(args...) = RecipesBase.apply_recipe(Dict{Symbol,Any}(), args...)
    @test length(series(w)) == 4 + 1 + 2 + 2           # axes, wedges, rings, arrow + halo
    @test length(series(c)) == 2
    @test length(series(o)) == 3
    @test length(series(cycleLengths(p))) == 2
    @test length(series(phaseValueCheck(p))) == 3
    @test length(series(p, :trace, 7)) == 2
    @test length(series(p, :portrait, 7)) == 1
    @test_throws ArgumentError series(p, :rose, 7)
end

include("physicell.jl")

end
