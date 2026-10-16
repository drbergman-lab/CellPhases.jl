# Readouts on a PhaseResult: cycle lengths, how well the phase tracks the value, phase coherence
# with distance and the local order parameter. Phase waves are in waves.jl.

_frames(p::PhaseResult, frames) = begin
    fr = collect(Int, frames)
    all(t -> 1 <= t <= nframes(p), fr) || throw(ArgumentError("frames must lie in 1:$(nframes(p))"))
    fr
end

#==================================================================================================#
# Cycle lengths

"""
    CycleLengths

The result of `CycleLengths(p::PhaseResult)`.

- `lengths`: time between successive upward crossings of φ = 0 (successive peaks).
- `ids`: the agent each cycle belongs to.
"""
struct CycleLengths{I}
    lengths::Vector{Float64}
    ids::Vector{I}
end

Base.show(io::IO, ::MIME"text/plain", c::CycleLengths) =
    print(io, "CycleLengths: $(length(c.lengths)) cycles from $(length(unique(c.ids))) agents",
          isempty(c.lengths) ? "" : ", median $(round(median(c.lengths); sigdigits = 4))")

"""
    CycleLengths(p::PhaseResult) -> CycleLengths

The time between successive peaks of each agent's band-passed signal: successive upward crossings
of φ = 0, each placed between frames by linear interpolation of the phase.

Only crossings within an unbroken stretch of phases count, so a masked frame or a gap never joins
two cycles into one.

# Examples
```jldoctest
julia> using Random

julia> rng = Xoshiro(1); xy = [(100.0rand(rng), 100.0rand(rng)) for _ in 1:100];

julia> positions, ids = fill(xy, 200), fill(1:100, 200);

julia> v = pulseTrain(ids; dt = 10, period = 170, rng = Xoshiro(2));

julia> c = CycleLengths(cellPhases(positions, v, ids; dt = 10, band = (100, 300)));

julia> round(median(c.lengths); digits = -1)
170.0
```
"""
function CycleLengths(p::PhaseResult{D,I}) where {D,I}
    lengths = Float64[]
    ids = I[]
    for tr in p.tracks
        crossings = Float64[]
        for k in 2:length(tr.frames)
            t0, t1 = tr.frames[k-1], tr.frames[k]
            a, b = p.phase[t0][tr.rows[k-1]], p.phase[t1][tr.rows[k]]
            if isnan(a) || isnan(b)
                _closerun!(lengths, ids, crossings, tr.id, p.dt)
                continue
            end
            a < 0 <= b && b - a < π && push!(crossings, t0 + (-a) / (b - a))
        end
        _closerun!(lengths, ids, crossings, tr.id, p.dt)
    end
    return CycleLengths{I}(lengths, ids)
end

function _closerun!(lengths, ids, crossings, id, dt)
    for d in diff(crossings)
        push!(lengths, d * dt)
        push!(ids, id)
    end
    empty!(crossings)
end

#==================================================================================================#
# Does the phase track the value?

"""
    PhaseValueCheck

The result of `PhaseValueCheck(p::PhaseResult)`.

- `phase`, `zscore`: every agent-frame's phase and its transformed value, z-scored within its
  track.
- `r2`: share of the variance of `zscore` explained by `cos φ` and `sin φ`.
- `peak`: the phase at which the fitted `zscore` peaks (0 when the value peaks with the band-passed
  signal).
- `centers`, `median`, `q25`, `q75`: `zscore` binned by phase.
"""
struct PhaseValueCheck
    phase::Vector{Float64}
    zscore::Vector{Float64}
    r2::Float64
    peak::Float64
    centers::Vector{Float64}
    median::Vector{Float64}
    q25::Vector{Float64}
    q75::Vector{Float64}
end

Base.show(io::IO, ::MIME"text/plain", c::PhaseValueCheck) =
    print(io, "PhaseValueCheck: $(length(c.zscore)) agent-frames; phase explains ",
          "$(round(100c.r2; digits = 1))% of the variance; peak at φ = $(round(rad2deg(c.peak); digits = 1))°")

"""
    PhaseValueCheck(p::PhaseResult; nbins = 36, minframes = 20) -> PhaseValueCheck

Check that the phase follows the value it was extracted from: regress each agent-frame's
transformed value, z-scored within its track, on `cos φ` and `sin φ`.

A high `r2` means the value is mostly a function of the phase. `peak` near 0 means the value peaks
where the band-passed signal does; a pulse-like value can peak elsewhere. Tracks with fewer than
`minframes` frames with a phase are left out.

# Examples
```jldoctest
julia> using Random

julia> rng = Xoshiro(1); xy = [(100.0rand(rng), 100.0rand(rng)) for _ in 1:100];

julia> positions, ids = fill(xy, 150), fill(1:100, 150);

julia> v = planeWave(positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> c = PhaseValueCheck(cellPhases(positions, v, ids; dt = 10, band = (100, 300)));

julia> c.r2 > 0.7, abs(rad2deg(c.peak)) < 5
(true, true)
```
"""
function PhaseValueCheck(p::PhaseResult; nbins::Integer = 36, minframes::Integer = 20)
    φall = Float64[]
    zscores = Float64[]
    for tr in p.tracks
        sel = [(t, r) for (t, r) in zip(tr.frames, tr.rows) if !isnan(p.phase[t][r])]
        length(sel) >= minframes || continue
        x = Float64.(p.transform.([p.values[t][r] for (t, r) in sel]))
        σ = std(x)
        σ > 0 || continue
        append!(zscores, (x .- mean(x)) ./ σ)
        append!(φall, [p.phase[t][r] for (t, r) in sel])
    end
    return _phasevalue(φall, zscores, nbins)
end

function _phasevalue(φ, zscore, nbins)
    edges = range(-π, π; length = nbins + 1)
    centers = collect((edges[1:end-1] .+ edges[2:end]) ./ 2)
    if length(zscore) < 3
        nanv = fill(NaN, nbins)
        return PhaseValueCheck(φ, zscore, NaN, NaN, centers, nanv, copy(nanv), copy(nanv))
    end
    X = [ones(length(φ)) cos.(φ) sin.(φ)]
    β = X \ zscore
    r2 = 1 - sum(abs2, zscore - X * β) / sum(abs2, zscore .- mean(zscore))
    bin = clamp.(floor.(Int, (φ .+ π) ./ (2π / nbins)) .+ 1, 1, nbins)
    stat(f) = [let s = zscore[bin .== k]; isempty(s) ? NaN : f(s) end for k in 1:nbins]
    return PhaseValueCheck(φ, zscore, r2, atan(β[3], β[2]), centers, stat(median),
                           stat(s -> quantile(s, 0.25)), stat(s -> quantile(s, 0.75)))
end

#==================================================================================================#
# Coherence with distance

"""
    PhaseCoherence

The result of `PhaseCoherence(p::PhaseResult)`. Per distance bin:

- `centers`: bin centres (position units); bin `k` holds pairs at distance within `dr/2` of
  `centers[k]`.
- `coherence`: mean `cos(φᵢ − φⱼ)` over same-frame pairs.
- `null`: the same, averaged over time-shift surrogates; `excess = coherence − null`.
- `nulllow`, `nullhigh`: 2.5% and 97.5% quantiles of the surrogate values.
- `pvalue`: one-sided, `(1 + #surrogates ≥ coherence) / (1 + nsurrogates)`.
- `counts`: pairs behind each value; `surrogates`: one column per surrogate.
"""
struct PhaseCoherence
    centers::Vector{Float64}
    coherence::Vector{Float64}
    null::Vector{Float64}
    excess::Vector{Float64}
    nulllow::Vector{Float64}
    nullhigh::Vector{Float64}
    pvalue::Vector{Float64}
    counts::Vector{Int}
    surrogates::Matrix{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", c::PhaseCoherence)
    println(io, "PhaseCoherence: $(length(c.centers)) distance bins to $(last(c.centers)), ",
            "$(size(c.surrogates, 2)) time-shift surrogates")
    print(io, "  excess coherence: ", round.(c.excess; digits = 2))
end

"""
    PhaseCoherence(p::PhaseResult; dr, rmax, frames = 1:nframes, nsurrogates = 19,
                   rng = Random.default_rng()) -> PhaseCoherence

How far apart agents can be and still be in step: the mean `cos(φᵢ − φⱼ)` over pairs of agents in
the same frame, binned by distance (bins of width `dr` centred on `0, dr, …, rmax`), with the same
quantity on [`timeShift`](@ref) surrogates as the null.

`frames` picks which frames to pool. Phases are band-passed, so an agent's phase moves by about
`2π·dt/period` per frame; when the band's shortest period spans many frames, consecutive frames
are strongly correlated and a sparser selection (`1:2:n`, `1:4:n`) costs less for little loss.
Compare two selections to check that on your data.

# Examples
```jldoctest
julia> using Random

julia> rng = Xoshiro(1); xy = [(150.0rand(rng), 150.0rand(rng)) for _ in 1:300];

julia> positions, ids = fill(xy, 150), fill(1:300, 150);

julia> v = planeWave(positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> p = cellPhases(positions, v, ids; dt = 10, band = (100, 300));

julia> c = PhaseCoherence(p; dr = 10, rmax = 40, frames = 1:4:150, rng = Xoshiro(2));

julia> c.excess[2] > 0.5, all(<(0.1), abs.(c.null))
(true, true)
```
"""
function PhaseCoherence(p::PhaseResult; dr::Real, rmax::Real, frames = 1:nframes(p),
                        nsurrogates::Integer = 19, rng::AbstractRNG = Random.default_rng())
    dr > 0 && rmax >= 0 || throw(ArgumentError("dr must be positive and rmax non-negative"))
    centers = collect(0:Float64(dr):Float64(rmax))
    nb = length(centers)
    sets = _phasesets(p, nsurrogates, rng)
    S = zeros(nb, length(sets))
    N = zeros(Int, nb, length(sets))
    for t in _frames(p, frames)
        ok = findall(!isnan, p.phase[t])
        length(ok) >= 2 || continue
        pts = p.positions[t][ok]
        cs = [_cossin(Φ[t][ok]) for Φ in sets]
        tree = KDTree(pts)
        for a in eachindex(pts), b in inrange(tree, pts[a], rmax + dr / 2)
            b > a || continue
            k = round(Int, norm(pts[b] - pts[a]) / dr) + 1
            k <= nb || continue
            for (s, (c, si)) in enumerate(cs)
                (isnan(c[a]) || isnan(c[b])) && continue
                S[k, s] += c[a] * c[b] + si[a] * si[b]
                N[k, s] += 1
            end
        end
    end
    M = S ./ N
    return _coherence(centers, M[:, 1], M[:, 2:end], N[:, 1])
end

function _coherence(centers, obs, sur, counts)
    null, lo, hi, pv = _nullstats(obs, sur)
    return PhaseCoherence(centers, obs, null, obs .- null, lo, hi, pv, counts, sur)
end

# The observed phases followed by `n` time-shift surrogates, as per-frame phase arrays.
_phasesets(p, n, rng) = pushfirst!([_shifted(p, rng)[1] for _ in 1:n], p.phase)

_cossin(φ) = (cos.(φ), sin.(φ))

function _nullstats(obs, sur)
    n = length(obs)
    null = fill(NaN, n); lo = fill(NaN, n); hi = fill(NaN, n); pv = fill(NaN, n)
    for k in 1:n
        s = filter(!isnan, sur[k, :])
        isempty(s) && continue
        null[k] = mean(s)
        lo[k], hi[k] = quantile(s, (0.025, 0.975))
        pv[k] = (1 + count(>=(obs[k]), s)) / (1 + length(s))
    end
    return null, lo, hi, pv
end

#==================================================================================================#
# Local order parameter

"""
    OrderParameter

The result of `OrderParameter(p::PhaseResult)`. Per radius:

- `radii`; `R`: the mean over agent-frames of the local order parameter `|mean e^{iφ}|` over the
  agents within that radius (centre included).
- `null`, `nulllow`, `nullhigh`, `pvalue`, `surrogates`: as in [`PhaseCoherence`](@ref), from
  time-shift surrogates; `excess = R − null`.
- `counts`: agent-frames behind each value.
"""
struct OrderParameter
    radii::Vector{Float64}
    R::Vector{Float64}
    null::Vector{Float64}
    excess::Vector{Float64}
    nulllow::Vector{Float64}
    nullhigh::Vector{Float64}
    pvalue::Vector{Float64}
    counts::Vector{Int}
    surrogates::Matrix{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", o::OrderParameter)
    println(io, "OrderParameter: radii $(o.radii), $(size(o.surrogates, 2)) time-shift surrogates")
    println(io, "  R:    ", round.(o.R; digits = 2))
    print(io, "  null: ", round.(o.null; digits = 2))
end

"""
    OrderParameter(p::PhaseResult; radii, frames = 1:nframes, minagents = 3, nsurrogates = 19,
                   rng = Random.default_rng()) -> OrderParameter

Local synchrony around each agent: `R = |mean e^{iφ}|` over the agents within radius `r` of it
(itself included), averaged over agents and frames, for each `r` in `radii`. Discs holding fewer
than `minagents` agents with a phase are skipped.

`R` is 1 when all neighbours share one phase. Unrelated phases give `R ≈ 1/√n` for `n` agents,
not 0, so compare `R` with `null`, computed the same way on [`timeShift`](@ref) surrogates.

# Examples
```jldoctest
julia> using Random

julia> rng = Xoshiro(1); xy = [(150.0rand(rng), 150.0rand(rng)) for _ in 1:300];

julia> positions, ids = fill(xy, 150), fill(1:300, 150);

julia> v = planeWave(positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> p = cellPhases(positions, v, ids; dt = 10, band = (100, 300));

julia> o = OrderParameter(p; radii = [10, 20, 40], frames = 1:4:150, rng = Xoshiro(2));

julia> all(o.R .> o.nullhigh)
true
```
"""
function OrderParameter(p::PhaseResult; radii, frames = 1:nframes(p), minagents::Integer = 3,
                        nsurrogates::Integer = 19, rng::AbstractRNG = Random.default_rng())
    radii = sort!(Float64.(collect(radii)))
    isempty(radii) || first(radii) >= 0 || throw(ArgumentError("radii must be non-negative"))
    nr = length(radii)
    sets = _phasesets(p, nsurrogates, rng)
    S = zeros(nr, length(sets))
    N = zeros(Int, nr, length(sets))
    for t in _frames(p, frames)
        ok = findall(!isnan, p.phase[t])
        isempty(ok) && continue
        pts = p.positions[t][ok]
        phasors = [cis.(Φ[t][ok]) for Φ in sets]
        tree = KDTree(pts)
        for a in eachindex(pts)
            nb = inrange(tree, pts[a], last(radii))
            d = [norm(pts[b] - pts[a]) for b in nb]
            for (s, e) in enumerate(phasors)
                isnan(e[a]) && continue
                for (j, r) in enumerate(radii)
                    acc = zero(ComplexF64)
                    n = 0
                    for (b, db) in zip(nb, d)
                        (db <= r && !isnan(e[b])) || continue
                        acc += e[b]
                        n += 1
                    end
                    n >= minagents || continue
                    S[j, s] += abs(acc) / n
                    N[j, s] += 1
                end
            end
        end
    end
    M = S ./ N
    obs = M[:, 1]
    sur = M[:, 2:end]
    null, lo, hi, pv = _nullstats(obs, sur)
    return OrderParameter(radii, obs, null, obs .- null, lo, hi, pv, N[:, 1], sur)
end
