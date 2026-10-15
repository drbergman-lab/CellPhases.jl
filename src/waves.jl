# Local phase waves: a plane fitted to the neighbours' phase differences gives the local phase
# gradient; the wave travels down the gradient at speed ω/|∇φ|, with ω the agent's own phase rate.

"""
    PhaseWaves

Returned by [`phaseWaves`](@ref). One entry per agent-frame with enough neighbours for a fit:

- `frame`, `id`, `position`: where the fit was made.
- `gradient`: the fitted local phase gradient ∇φ (radians per position unit).
- `r2`: the fit's R²; low values mean no single plane wave explains the neighbours.
- `direction`: the direction of travel, the angle of −∇φ in radians (0 = +x, π/2 = +y).
- `omega`: the agent's own phase advance to the next frame, radians per minute; `NaN` when it has
  no phase in the next frame.
- `speed`: `omega / |∇φ|`, position units per minute, relative to the agents; `NaN` where `omega`
  is missing or not positive.
- `r`, `minneighbours`: the settings.

Report [`cleanShare`](@ref) next to any summary of directions or speeds.
"""
struct PhaseWaves{I}
    frame::Vector{Int}
    id::Vector{I}
    position::Vector{SVector{2,Float64}}
    gradient::Vector{SVector{2,Float64}}
    r2::Vector{Float64}
    direction::Vector{Float64}
    omega::Vector{Float64}
    speed::Vector{Float64}
    r::Float64
    minneighbours::Int
end

function Base.show(io::IO, ::MIME"text/plain", w::PhaseWaves)
    print(io, "PhaseWaves: $(length(w.r2)) local fits (r = $(w.r)); ",
          "$(round(100cleanShare(w); digits = 1))% with R² ≥ 0.5")
end

"""
    cleanShare(w::PhaseWaves; minr2 = 0.5)

The share of local fits with `r2 ≥ minr2`: how much of the tissue any direction or speed summary
speaks for.

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 300, L = 150.0, nframes = 150, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 90,
                     rng = Xoshiro(2));

julia> w = phaseWaves(cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300)); r = 30);

julia> cleanShare(w) > 0.9
true
```
"""
cleanShare(w::PhaseWaves; minr2::Real = 0.5) =
    isempty(w.r2) ? NaN : count(>=(minr2), w.r2) / length(w.r2)

"""
    phaseWaves(p::PhaseResult; r, minneighbours = 4, frames = 1:nframes) -> PhaseWaves

Fit a local plane wave around every agent with a phase: least squares of the wrapped phase
differences to the neighbours within `r` on their displacements gives ∇φ. The wave travels along
−∇φ at speed `ω/|∇φ|`, where `ω` is the agent's own phase advance to the next frame. 2D only.

Speeds are relative to the agents (ω follows the agent). Wrapped differences are only unambiguous
while neighbours are less than half a wavelength apart, so waves shorter than about `2r` are
misread; choose `r` well below the shortest wavelength of interest.

# Examples
```jldoctest
julia> using Random, Statistics

julia> field = agentField(; n = 300, L = 150.0, nframes = 150, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 90,
                     rng = Xoshiro(2));

julia> w = phaseWaves(cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300)); r = 30);

julia> good = w.r2 .>= 0.5;

julia> round(rad2deg(angle(mean(cis.(w.direction[good])))); digits = -1)
90.0

julia> round(median(filter(!isnan, w.speed[good])); digits = 1)
0.5
```
"""
function phaseWaves(p::PhaseResult{D,I}; r::Real, minneighbours::Integer = 4,
                    frames = 1:nframes(p)) where {D,I}
    D == 2 || throw(ArgumentError("phaseWaves needs 2D positions; got $(D)D"))
    r > 0 || throw(ArgumentError("r must be positive, got $r"))
    minneighbours >= 2 || throw(ArgumentError("minneighbours must be at least 2 for a plane fit"))
    next = [zeros(Int, length(φ)) for φ in p.phase]
    for tr in p.tracks, k in 1:length(tr.frames)-1
        next[tr.frames[k]][tr.rows[k]] = tr.rows[k+1]
    end
    out = PhaseWaves{I}(Int[], I[], SVector{2,Float64}[], SVector{2,Float64}[], Float64[],
                        Float64[], Float64[], Float64[], Float64(r), Int(minneighbours))
    for t in _frames(p, frames)
        ok = findall(!isnan, p.phase[t])
        length(ok) > minneighbours || continue
        pts = p.positions[t][ok]
        φ = p.phase[t][ok]
        tree = KDTree(pts)
        for a in eachindex(pts)
            nb = filter(!=(a), inrange(tree, pts[a], r))
            length(nb) >= minneighbours || continue
            g, q = _planefit(pts, φ, a, nb)
            all(isfinite, g) || continue
            ω = NaN
            if t < nframes(p) && next[t][ok[a]] > 0
                φn = p.phase[t+1][next[t][ok[a]]]
                isnan(φn) || (ω = angle(cis(φn - φ[a])) / p.dt)
            end
            push!(out.frame, t)
            push!(out.id, p.ids[t][ok[a]])
            push!(out.position, pts[a])
            push!(out.gradient, g)
            push!(out.r2, q)
            push!(out.direction, atan(-g[2], -g[1]))
            push!(out.omega, ω)
            push!(out.speed, ω > 0 ? ω / norm(g) : NaN)
        end
    end
    return out
end

function _planefit(pts, φ, a, nb)
    X = Matrix{Float64}(undef, length(nb), 2)
    y = Vector{Float64}(undef, length(nb))
    for (k, b) in enumerate(nb)
        d = pts[b] - pts[a]
        X[k, 1], X[k, 2] = d[1], d[2]
        y[k] = angle(cis(φ[b] - φ[a]))
    end
    g = X \ y
    q = 1 - sum(abs2, y - X * g) / max(sum(abs2, y), eps())
    return SVector{2,Float64}(g[1], g[2]), q
end
