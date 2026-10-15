# Steps 1–3 of the pipeline: per-agent traces from per-frame input, the band-pass filter, and the
# Hilbert phase. Everything downstream reads a `PhaseResult`.

"""
    Track{I}

One run of consecutive frames of one agent that received a phase.

- `id::I`: the agent's ID.
- `frames::UnitRange{Int}`: the frames with a phase (the run with `trim` cut from each end).
- `rows::Vector{Int}`: the agent's row in each of those frames.

An agent that is missing from a frame, or changes into a cell type that was not selected, has its
trace split there; each piece long enough becomes its own `Track`.
"""
struct Track{I}
    id::I
    frames::UnitRange{Int}
    rows::Vector{Int}
end

"""
    PhaseResult{D,I}

The phase of every tracked agent in every frame, returned by [`cellPhases`](@ref).

All per-frame fields are aligned with the input: `phase[t][i]` belongs to the agent in row `i` of
frame `t`.

- `phase`: Hilbert phase in (−π, π]; 0 is the peak of the band-passed signal. `NaN` where the
  agent has no phase (track too short, inside the trimmed ends, or masked).
- `amplitude`: Hilbert amplitude of the band-passed, transformed signal; `NaN` where unavailable.
- `positions`, `ids`, `values`: the input, with positions as `SVector{D,Float64}` and values as
  given (before the transform).
- `tracks`: the [`Track`](@ref)s that received a phase.
- Settings: `dt`, `band`, `order`, `trim`, `minlength` (all times in the unit of `dt`),
  `transform`, `mask`.
"""
struct PhaseResult{D,I}
    phase::Vector{Vector{Float64}}
    amplitude::Vector{Vector{Float64}}
    positions::Vector{Vector{SVector{D,Float64}}}
    ids::Vector{Vector{I}}
    values::Vector{Vector{Float64}}
    tracks::Vector{Track{I}}
    dt::Float64
    band::Tuple{Float64,Float64}
    order::Int
    trim::Float64
    minlength::Float64
    transform::Any
    mask::Float64
end

function Base.show(io::IO, ::MIME"text/plain", p::PhaseResult{D}) where {D}
    nφ = sum(φ -> count(!isnan, φ), p.phase; init = 0)
    nall = sum(length, p.phase; init = 0)
    println(io, "PhaseResult: $(length(p.phase)) frames, dt = $(p.dt), $(D)D")
    println(io, "  tracks: $(length(p.tracks)) from $(length(unique(t.id for t in p.tracks))) agents")
    println(io, "  agent-frames with a phase: $nφ of $nall")
    print(io, "  band $(p.band[1])–$(p.band[2]), Butterworth order $(p.order), trim $(p.trim), ",
          "minlength $(p.minlength)", p.mask > 0 ? ", mask $(p.mask)" : "")
end

"""
    bandpassFilter(band; dt, order = 3)

The band-pass filter [`cellPhases`](@ref) applies: a Butterworth filter of the given `order`.

- `dt`: the frame interval.
- `band = (pmin, pmax)`: the shortest and longest **periods** to keep, in the same time unit as
  `dt`. Any unit works; it only has to be the same for both.

The conversion is explicit: the cutoff frequencies are `1/pmax` and `1/pmin` cycles per time unit
and the sampling rate is `fs = 1/dt` samples per time unit, and `fs` is always passed to
`DSP.digitalfilter`. Normalised frequencies are never used. With `dt = 10` and
`band = (100, 300)`: `fs = 0.1`, cutoffs `1/300` and `1/100` cycles per time unit.

The cutoffs are the half-power points of one pass (amplitude gain `1/√2`). `cellPhases` applies
the filter forward and backward (`filtfilt`), which squares the gain, so as applied the cutoffs
are half-**amplitude** points.

`pmin` must exceed `2dt`: the shortest period a frame interval can resolve is two frames.

# Examples
```jldoctest
julia> using DSP

julia> f = bandpassFilter((100, 300); dt = 10);

julia> periods = [100, 170, 300];   # same unit as dt

julia> H = freqresp(f, 2π .* (1 ./ periods) .* 10);   # radians per sample = 2π · (1/period) · dt

julia> round.(abs2.(H); digits = 3)   # amplitude gain as applied (forward and backward)
3-element Vector{Float64}:
 0.5
 1.0
 0.5
```
"""
function bandpassFilter(band; dt::Real, order::Integer = 3)
    length(band) == 2 || throw(ArgumentError("band must be (pmin, pmax): two periods, got $band"))
    pmin, pmax = Float64(band[1]), Float64(band[2])          # periods, in the unit of dt
    dt > 0 || throw(ArgumentError("dt (the frame interval) must be positive, got $dt"))
    0 < pmin < pmax || throw(ArgumentError(
        "band must be (pmin, pmax) periods with 0 < pmin < pmax, got $band"))
    pmin > 2dt || throw(ArgumentError(
        "the shortest period in the band (pmin = $pmin) must exceed two frames " *
        "(2dt = $(2dt)): shorter periods cannot be resolved at this frame interval"))
    order >= 1 || throw(ArgumentError("order must be at least 1, got $order"))
    flow, fhigh = 1 / pmax, 1 / pmin                          # cutoffs, cycles per time unit
    fs = 1 / dt                                               # sampling rate, samples per time unit
    return digitalfilter(Bandpass(flow, fhigh), Butterworth(order); fs = fs)
end

# Frames for the slowest pole of the filter to decay to 1/e: how far the filter's memory reaches.
_filtermemory(f) = -1 / log(maximum(abs, f.p))

"""
    cellPhases(positions, values, ids; dt, band, order = 3, trim = nothing, minlength = nothing,
               transform = log, mask = 0.0) -> PhaseResult

Give every tracked agent a phase for its oscillating value, frame by frame.

Inputs are per frame: `positions[t]` (points as tuples or vectors), `values[t]` and `ids[t]`,
one entry per agent present in frame `t`. Frames must be equally spaced, `dt` apart. All times
(`dt`, `band`, `trim`, `minlength`) are in one time unit of your choice.

For each agent, consecutive frames form a run; runs of at least `minlength` are
transformed (`log` by default), mean-subtracted, band-passed (zero-phase Butterworth of `order`,
periods `band = (pmin, pmax)`; see [`bandpassFilter`](@ref)), and
Hilbert-transformed. A time `trim` is dropped from each end of the run, where the filter has not
settled.

`band = (pmin, pmax)` are periods and `dt` the frame interval, in the same time unit; see
[`bandpassFilter`](@ref) for how they become cutoffs.

- `trim` defaults to the filter's memory (frames for its slowest pole to decay to 1/e), rounded up
  to whole frames.
- `minlength` defaults to `2trim + 2pmax`: two of the slowest cycles left after trimming. It
  cannot be set below `2pmax`.
- `mask = q` sets the phase to `NaN` in each run's frames whose amplitude is below that run's
  `q`-quantile. Off (`0`) by default; see the manual for when to use it.

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 200, L = 100.0, nframes = 150, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 90,
                     rng = Xoshiro(2));

julia> p = cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300));

julia> length(p.tracks), p.trim
(200, 180.0)
```
"""
function cellPhases(positions::AbstractVector, values::AbstractVector, ids::AbstractVector;
                    dt::Real, band, order::Integer = 3, trim::Union{Real,Nothing} = nothing,
                    minlength::Union{Real,Nothing} = nothing, transform = log, mask::Real = 0.0)
    nT = length(values)
    length(positions) == nT == length(ids) || throw(DimensionMismatch(
        "positions, values and ids must have one entry per frame; got $(length(positions)), $nT, $(length(ids))"))
    for t in 1:nT
        length(positions[t]) == length(values[t]) == length(ids[t]) || throw(DimensionMismatch(
            "frame $t has $(length(positions[t])) positions, $(length(values[t])) values and $(length(ids[t])) ids"))
    end
    0 <= mask < 1 || throw(ArgumentError("mask must be a quantile in [0, 1), got $mask"))
    filt = bandpassFilter(band; dt, order)
    pmax = Float64(Tuple(band)[2])
    trimframes = trim === nothing ? ceil(Int, _filtermemory(filt)) : round(Int, trim / dt)
    trimframes >= 0 || throw(ArgumentError("trim must be non-negative, got $trim"))
    minlength = minlength === nothing ? 2trimframes * dt + 2pmax : Float64(minlength)
    minlength >= 2pmax || throw(ArgumentError(
        "minlength ($minlength) must be at least two of the longest periods in the band " *
        "(2pmax = $(2pmax)); shorter tracks cannot carry that rhythm"))
    minframes = max(ceil(Int, minlength / dt - 1e-9), 2trimframes + 1, 12order + 1)
    nT * dt < minlength && @warn "the recording lasts $(nT * dt), shorter than minlength " *
        "($minlength = 2trim + 2pmax by default); no agent can receive a phase"

    D = _dims(positions)
    pos = [SVector{D,Float64}[_point(x, D) for x in positions[t]] for t in 1:nT]
    I = mapreduce(eltype, promote_type, ids; init = Union{})
    I === Union{} && (I = Int)
    idv = Vector{I}[convert(Vector{I}, collect(ids[t])) for t in 1:nT]
    vals = Vector{Float64}[Float64.(values[t]) for t in 1:nT]

    Φ = [fill(NaN, length(vals[t])) for t in 1:nT]
    A = [fill(NaN, length(vals[t])) for t in 1:nT]
    tracks = Track{I}[]
    for (id, (f, rows)) in _byagent(idv)
        for seg in _runs(f)
            length(seg) >= minframes || continue
            v = [vals[f[k]][rows[k]] for k in seg]
            x = _transform(transform, v, id)
            z = hilbert(filtfilt(filt, x .- mean(x)))
            keep = trimframes+1:length(seg)-trimframes
            amp = abs.(z[keep])
            cut = mask > 0 ? quantile(amp, mask) : -Inf
            for (j, k) in zip(keep, seg[keep])
                A[f[k]][rows[k]] = abs(z[j])
                Φ[f[k]][rows[k]] = abs(z[j]) < cut ? NaN : angle(z[j])
            end
            push!(tracks, Track{I}(id, f[seg[first(keep)]]:f[seg[last(keep)]], rows[seg[keep]]))
        end
    end
    sort!(tracks; by = tr -> (first(tr.frames), tr.rows[1]))
    return PhaseResult{D,I}(Φ, A, pos, idv, vals, tracks, Float64(dt), Float64.(Tuple(band)),
                            Int(order), trimframes * Float64(dt), Float64(minlength), transform,
                            Float64(mask))
end

function _dims(positions)
    for frame in positions, x in frame
        D = length(x)
        D >= 1 || break
        return D
    end
    throw(ArgumentError("no positions: every frame is empty"))
end

_point(x, D) = length(x) == D ? SVector{D,Float64}(Tuple(x)) :
    throw(DimensionMismatch("positions mix $(D)D points with $(length(x))D points"))

# id => (frames, rows), frames increasing.
function _byagent(ids::Vector{Vector{I}}) where {I}
    tracks = Dict{I,Tuple{Vector{Int},Vector{Int}}}()
    for t in eachindex(ids), (r, id) in enumerate(ids[t])
        f, rows = get!(() -> (Int[], Int[]), tracks, id)
        !isempty(f) && f[end] == t && throw(ArgumentError("ID $(repr(id)) appears twice in frame $t"))
        push!(f, t)
        push!(rows, r)
    end
    return tracks
end

# Index ranges of consecutive frames within the increasing frame list `f`.
function _runs(f)
    breaks = [0; findall(!=(1), diff(f)); length(f)]
    return [breaks[b]+1:breaks[b+1] for b in 1:length(breaks)-1]
end

function _transform(transform, v, id)
    transform === log && any(<=(0), v) && throw(ArgumentError(
        "agent $(repr(id)) has a value ≤ 0, which `log` cannot take; pass a different `transform` " *
        "(e.g. `identity`, or `x -> log(x + c)`)"))
    x = Float64.(transform.(v))
    all(isfinite, x) || throw(ArgumentError("the transform gave a non-finite value for agent $(repr(id))"))
    return x
end

nframes(p::PhaseResult) = length(p.phase)
