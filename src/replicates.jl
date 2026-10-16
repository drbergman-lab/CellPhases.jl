# Readouts over replicate simulations: one result per replicate plus a pooled result, so conditions
# can be compared with the spread between replicates rather than within one run.

"""
    Replicates{T}

A readout run on each of several [`PhaseResult`](@ref)s (replicate simulations), returned by
the lowercase readout functions: `phaseCoherence(runs; …)`, `orderParameter`, `phaseWaves`,
`cycleLengths`, `phaseValueCheck`.

- `results::Vector{T}`: one result per replicate, in input order.
- `pooled::T`: the readout on all replicates together. Cycle lengths, phase-value pairs and wave
  fits are concatenated; coherence and order parameter are averaged with each replicate weighted by
  its pair (or agent-frame) counts, surrogate by surrogate.

Use [`replicateSpread`](@ref) for the between-replicate spread of any per-bin field.
"""
struct Replicates{T}
    results::Vector{T}
    pooled::T
end

Base.show(io::IO, m::MIME"text/plain", r::Replicates{T}) where {T} =
    (println(io, "Replicates{$(nameof(T))}: $(length(r.results)) replicates; pooled:"); show(io, m, r.pooled))

"""
    replicateSpread(r::Replicates, field::Symbol) -> (; mean, std, n)

Mean, standard deviation and number of replicates, per element, of the vector field `field`
(e.g. `:excess`, `:R`) across `r.results`. Replicates with `NaN` in an element are left out of
that element.

# Examples
```jldoctest
julia> using Random

julia> runs = map(1:3) do k
           rng = Xoshiro(k); xy = [(100rand(rng), 100rand(rng)) for _ in 1:150]
           positions, ids = fill(xy, 150), fill(1:150, 150)
           v = planeWave(positions; dt = 10, period = 180, speed = 0.5, direction = 0, rng)
           cellPhases(positions, v, ids; dt = 10, band = (100, 300))
       end;

julia> o = orderParameter(runs; radii = [20], frames = 1:4:150, nsurrogates = 4, rng = Xoshiro(9));

julia> s = replicateSpread(o, :R); s.n
1-element Vector{Int64}:
 3
```
"""
function replicateSpread(r::Replicates, field::Symbol)
    cols = [Float64.(getfield(x, field)) for x in r.results]
    n = length(first(cols))
    all(c -> length(c) == n, cols) || throw(DimensionMismatch("$field has different lengths across replicates"))
    vals = [filter(!isnan, [c[k] for c in cols]) for k in 1:n]
    return (mean = [isempty(v) ? NaN : mean(v) for v in vals],
            std = [length(v) < 2 ? NaN : std(v) for v in vals],
            n = length.(vals))
end

const _Runs = AbstractVector{<:PhaseResult}

"""
    cycleLengths(runs::AbstractVector{<:PhaseResult}) -> Replicates{CycleLengths}

[`CycleLengths`](@ref) for each of several replicate runs, plus a pooled result; see
[`Replicates`](@ref).
"""
function cycleLengths(ps::_Runs)
    rs = [CycleLengths(p) for p in ps]
    return Replicates(rs, CycleLengths(reduce(vcat, (r.lengths for r in rs)),
                                       reduce(vcat, (r.ids for r in rs))))
end

"""
    phaseValueCheck(runs::AbstractVector{<:PhaseResult}; kwargs...) -> Replicates{PhaseValueCheck}

[`PhaseValueCheck`](@ref) for each of several replicate runs, plus a pooled result; see
[`Replicates`](@ref). Keywords are those of `PhaseValueCheck(p::PhaseResult; …)`.
"""
function phaseValueCheck(ps::_Runs; nbins::Integer = 36, kwargs...)
    rs = [PhaseValueCheck(p; nbins, kwargs...) for p in ps]
    return Replicates(rs, _phasevalue(reduce(vcat, (r.phase for r in rs)),
                                      reduce(vcat, (r.zscore for r in rs)), nbins))
end

"""
    phaseWaves(runs::AbstractVector{<:PhaseResult}; kwargs...) -> Replicates{PhaseWaves}

[`PhaseWaves`](@ref) for each of several replicate runs, plus a pooled result; see
[`Replicates`](@ref). Keywords are those of `PhaseWaves(p::PhaseResult; …)`.
"""
function phaseWaves(ps::_Runs; kwargs...)
    rs = [PhaseWaves(p; kwargs...) for p in ps]
    cat(f) = reduce(vcat, (getfield(r, f) for r in rs))
    pooled = typeof(first(rs))(cat(:frame), cat(:id), cat(:position), cat(:gradient), cat(:r2),
                               cat(:direction), cat(:omega), cat(:speed), first(rs).r,
                               first(rs).minneighbours)
    return Replicates(rs, pooled)
end

"""
    phaseCoherence(runs::AbstractVector{<:PhaseResult}; kwargs...) -> Replicates{PhaseCoherence}

[`PhaseCoherence`](@ref) for each of several replicate runs, plus a pooled result; see
[`Replicates`](@ref). Keywords are those of `PhaseCoherence(p::PhaseResult; …)`.
"""
function phaseCoherence(ps::_Runs; kwargs...)
    rs = [PhaseCoherence(p; kwargs...) for p in ps]
    obs, sur, counts = _poolweighted(rs, :coherence)
    return Replicates(rs, _coherence(first(rs).centers, obs, sur, counts))
end

"""
    orderParameter(runs::AbstractVector{<:PhaseResult}; kwargs...) -> Replicates{OrderParameter}

[`OrderParameter`](@ref) for each of several replicate runs, plus a pooled result; see
[`Replicates`](@ref). Keywords are those of `OrderParameter(p::PhaseResult; …)`.
"""
function orderParameter(ps::_Runs; kwargs...)
    rs = [OrderParameter(p; kwargs...) for p in ps]
    obs, sur, counts = _poolweighted(rs, :R)
    null, lo, hi, pv = _nullstats(obs, sur)
    return Replicates(rs, OrderParameter(first(rs).radii, obs, null, obs .- null, lo, hi, pv,
                                         counts, sur))
end

# Count-weighted mean of the observed field and of each surrogate column across replicates.
function _poolweighted(rs, field)
    counts = reduce(+, (r.counts for r in rs))
    wmean(get) = begin
        s = sum(r -> let x = get(r); ifelse.(isnan.(x), 0.0, x .* r.counts) end, rs)
        s ./ counts
    end
    obs = wmean(r -> getfield(r, field))
    sur = wmean(r -> r.surrogates)
    return obs, sur, counts
end
