# The null: time-shift surrogates. Each track's phase series is rotated in time by its own random
# offset, which keeps every agent's rhythm and destroys the alignment between agents.

"""
    timeShift(p::PhaseResult; rng = Random.default_rng()) -> PhaseResult

A time-shift surrogate of `p`: each track's phase and amplitude series is circularly shifted by an
independent random number of frames. Positions, IDs and values are unchanged.

Every agent keeps its own rhythm and amplitude, while phases at the same moment in neighbouring
agents become unrelated. [`PhaseCoherence`](@ref) and [`OrderParameter`](@ref) use it for their
null bands and p-values; any readout can be run on it directly.

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 50, L = 50.0, nframes = 120, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300));

julia> s = timeShift(p; rng = Xoshiro(2));

julia> sort(filter(!isnan, s.phase[60])) == sort(filter(!isnan, p.phase[60]))
false
```
"""
function timeShift(p::PhaseResult{D,I}; rng::AbstractRNG = Random.default_rng()) where {D,I}
    Φ, A = _shifted(p, rng)
    return PhaseResult{D,I}(Φ, A, p.positions, p.ids, p.values, p.tracks, p.dt, p.band, p.order,
                            p.trim, p.minlength, p.transform, p.mask)
end

function _shifted(p::PhaseResult, rng::AbstractRNG)
    Φ = [fill(NaN, length(φ)) for φ in p.phase]
    A = [fill(NaN, length(a)) for a in p.amplitude]
    for tr in p.tracks
        m = length(tr.frames)
        s = m > 1 ? rand(rng, 1:m-1) : 0
        for (j, (t, r)) in enumerate(zip(tr.frames, tr.rows))
            k = mod1(j + s, m)
            t2, r2 = tr.frames[k], tr.rows[k]
            Φ[t][r] = p.phase[t2][r2]
            A[t][r] = p.amplitude[t2][r2]
        end
    end
    return Φ, A
end
