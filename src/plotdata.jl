# Plot-ready data from the readouts: plain numbers, tuples and polygons, no plotting package.
# The recipes (recipes.jl, ext/) draw from these, and so can any hand-built figure.

"""
    roseData(w::PhaseWaves; bin = 20, minr2 = 0.5, rings = (5, 10)) -> NamedTuple

Everything needed to draw the directions of local phase waves as a rose diagram, from the fits
with `r2 ≥ minr2`. Directions are compass angles in degrees, 0 = right (+x), 90 = up (+y).

Wedges are `bin` degrees wide; a wedge's length is `√share` (share in percent), so its area is
proportional to the number of fits. Fields:

- `wedges`: one closed polygon per wedge, a vector of `(x, y)` tuples starting at the origin.
- `lo`, `share`: each wedge's start angle (degrees) and its share of fits (%).
- `speed`: median speed per wedge (position units per time unit; `NaN` for an empty wedge).
- `rings`: `(radius, label)` reference circles at the shares in `rings`; `labelangle`: the
  centre of the shortest wedge, where ring labels least overlap the wedges.
- `rlim`: radius of the outermost ring or wedge.
- `meandirection` (degrees), `resultant` (mean resultant length R), `arrow`: the tip of an arrow
  from the origin of length `R·rlim` in the mean direction.
- `compass`: `(x, y, label, halign, valign)` for the labels right, up, left, down.
- `up`, `down`: shares (%) within ±45° of up and of down; `clean`: [`cleanShare`](@ref) (%);
  `n`: number of fits drawn.

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 300, L = 150.0, nframes = 150, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 90,
                     rng = Xoshiro(2));

julia> w = PhaseWaves(cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300)); r = 30);

julia> rose = roseData(w);

julia> rose.lo[argmax(rose.share)], round(rose.meandirection; digits = -1), rose.up > 95
(80.0, 90.0, true)
```
"""
function roseData(w::PhaseWaves; bin::Real = 20, minr2::Real = 0.5, rings = (5, 10))
    360 % bin == 0 || throw(ArgumentError("bin must divide 360, got $bin"))
    good = w.r2 .>= minr2
    deg = rad2deg.(w.direction[good])
    speed = w.speed[good]
    lo = collect(-180.0:bin:180-bin)
    inwedge(l) = l .<= deg .< l + bin
    share = [isempty(deg) ? 0.0 : 100 * count(inwedge(l)) / length(deg) for l in lo]
    wspeed = [let s = filter(!isnan, speed[inwedge(l)]); isempty(s) ? NaN : median(s) end for l in lo]
    wedges = [_wedge(l, l + bin, sqrt(s)) for (l, s) in zip(lo, share)]
    rlim = sqrt(max(maximum(share; init = 0.0), maximum(rings; init = 0) * 1.25))
    m = isempty(deg) ? complex(NaN) : mean(cis.(deg2rad.(deg)))
    R, θm = abs(m), rad2deg(angle(m))
    c = rlim * 1.04
    compass = [(c, 0.0, "right", :left, :center), (0.0, c, "up", :center, :bottom),
               (-c, 0.0, "left", :right, :center), (0.0, -c, "down", :center, :top)]
    within(θ) = isempty(deg) ? NaN : 100 * count(d -> abs(mod(d - θ + 180, 360) - 180) < 45, deg) / length(deg)
    return (wedges = wedges, lo = lo, share = share, speed = wspeed,
            rings = [(sqrt(Float64(s)), "$(s)%") for s in rings],
            labelangle = lo[argmin(share)] + bin / 2, rlim = rlim,
            meandirection = θm, resultant = R, arrow = (R * rlim * cosd(θm), R * rlim * sind(θm)),
            compass = compass, up = within(90), down = within(-90),
            clean = 100 * cleanShare(w; minr2), n = length(deg))
end

function _wedge(lo, hi, r; n = 12)
    pts = [(0.0, 0.0); [(r * cosd(θ), r * sind(θ)) for θ in range(lo, hi; length = n)]]
    return push!(pts, (0.0, 0.0))
end

_circle(r; n = 181) = [(r * cosd(θ), r * sind(θ)) for θ in range(0, 360; length = n)]

"""
    coherenceData(c::PhaseCoherence; skipzero = true) -> (; x, y, low, high)

Excess phase coherence against distance, with the surrogate 95% band shifted by the same null
(`low`, `high`), so 0 is "no more in step than chance". `skipzero` drops the first bin, which holds
only the closest pairs.

# Examples
```jldoctest
julia> c = PhaseCoherence([0.0, 10.0], [0.9, 0.5], [0.0, 0.1], [0.9, 0.4], [-0.1, 0.0],
                          [0.1, 0.2], [0.05, 0.05], [10, 100], zeros(2, 0));

julia> coherenceData(c)
(x = [10.0], y = [0.4], low = [-0.1], high = [0.1])
```
"""
function coherenceData(c::PhaseCoherence; skipzero::Bool = true)
    k = skipzero ? (2:length(c.centers)) : eachindex(c.centers)
    return (x = c.centers[k], y = c.excess[k], low = c.nulllow[k] .- c.null[k],
            high = c.nullhigh[k] .- c.null[k])
end

"""
    orderParameterData(o::OrderParameter) -> (; radii, R, null, low, high)

The local order parameter and its surrogate null with the 95% band, against radius.

# Examples
```jldoctest
julia> o = OrderParameter([10.0], [0.8], [0.3], [0.5], [0.25], [0.35], [0.05], [100], zeros(1, 0));

julia> orderParameterData(o).R
1-element Vector{Float64}:
 0.8
```
"""
orderParameterData(o::OrderParameter) =
    (radii = o.radii, R = o.R, null = o.null, low = o.nulllow, high = o.nullhigh)

"""
    cycleLengthData(c::CycleLengths; edges = nothing) -> (; edges, counts, median)

A histogram of cycle lengths. `edges` defaults to 30 bins over the central 98% of the lengths.

# Examples
```jldoctest
julia> cycleLengthData(CycleLengths([150.0, 170.0, 175.0, 260.0], [1, 1, 2, 2]);
                       edges = 100:50:300)
(edges = [100.0, 150.0, 200.0, 250.0, 300.0], counts = [0, 3, 0, 1], median = 172.5)
```
"""
function cycleLengthData(c::CycleLengths; edges = nothing)
    L = c.lengths
    if edges === nothing
        lo, hi = isempty(L) ? (0.0, 1.0) : quantile(L, (0.01, 0.99))
        hi > lo || (hi = lo + 1)
        edges = range(lo, hi; length = 31)
    end
    e = Float64.(collect(edges))
    counts = [count(x -> e[k] <= x < e[k+1], L) for k in 1:length(e)-1]
    return (edges = e, counts = counts, median = isempty(L) ? NaN : median(L))
end

"""
    phaseValueData(v::PhaseValueCheck; zscoreedges = range(-2.5, 3.5; length = 41))
        -> (; phaseedges, zscoreedges, density, centers, median, q25, q75, r2, peak)

The value against phase as a density: `density[i, j]` is the share of agent-frames in phase bin
`i` whose z-scored value lies in `zscoreedges[j:j+1]` (each phase column sums to at most 1), plus the
binned median and quartiles from [`PhaseValueCheck`](@ref).

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 50, L = 50.0, nframes = 120, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> d = phaseValueData(PhaseValueCheck(cellPhases(field.positions, v, field.ids; dt = 10,
                                                     band = (100, 300))));

julia> size(d.density), round(maximum(sum(d.density; dims = 2)); digits = 2)   # each phase column sums to ≤ 1
((36, 40), 1.0)
```
"""
function phaseValueData(v::PhaseValueCheck; zscoreedges = range(-2.5, 3.5; length = 41))
    nb = length(v.centers)
    pe = collect(range(-π, π; length = nb + 1))
    ze = Float64.(collect(zscoreedges))
    H = zeros(nb, length(ze) - 1)
    colsum = zeros(Int, nb)
    for (φ, zs) in zip(v.phase, v.zscore)
        i = clamp(floor(Int, (φ + π) / (2π / nb)) + 1, 1, nb)
        colsum[i] += 1
        j = searchsortedlast(ze, zs)
        1 <= j < length(ze) && (H[i, j] += 1)
    end
    density = H ./ max.(colsum, 1)
    return (phaseedges = pe, zscoreedges = ze, density = density, centers = v.centers, median = v.median,
            q25 = v.q25, q75 = v.q75, r2 = v.r2, peak = v.peak)
end

"""
    traceData(p::PhaseResult, id) -> (; time, value, phase, amplitude, signal)

One agent's trace over every frame it appears in: `time` (from the first frame, in the unit of `dt`), the input
`value`, its `phase` and `amplitude` (`NaN` outside its tracks), and the band-passed `signal`
`amplitude·cos(phase)`.

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 5, L = 50.0, nframes = 120, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> tr = traceData(cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300)), 3);

julia> length(tr.time), count(!isnan, tr.phase)
(120, 84)
```
"""
function traceData(p::PhaseResult, id)
    sel = [(t, r) for t in eachindex(p.ids) for r in findall(==(id), p.ids[t])]
    isempty(sel) && throw(ArgumentError("no agent with ID $(repr(id))"))
    φ = [p.phase[t][r] for (t, r) in sel]
    A = [p.amplitude[t][r] for (t, r) in sel]
    return (time = [(t - 1) * p.dt for (t, _) in sel], value = [p.values[t][r] for (t, r) in sel],
            phase = φ, amplitude = A, signal = A .* cos.(φ))
end

"""
    portraitData(p::PhaseResult, id) -> (; time, re, im)

One agent's analytic signal in the complex plane, over the frames where it has a phase: `re` is
the band-passed signal and `im` its Hilbert quadrature. The point circles the origin
counter-clockwise once per cycle, and φ is its angle.

# Examples
```jldoctest
julia> using Random

julia> field = agentField(; n = 5, L = 50.0, nframes = 120, rng = Xoshiro(1));

julia> v = planeWave(field.positions; dt = 10, period = 180, speed = 0.5, direction = 0,
                     rng = Xoshiro(2));

julia> length(portraitData(cellPhases(field.positions, v, field.ids; dt = 10, band = (100, 300)), 3).re)
84
```
"""
function portraitData(p::PhaseResult, id)
    tr = traceData(p, id)
    k = findall(!isnan, tr.phase)
    return (time = tr.time[k], re = tr.signal[k], im = tr.amplitude[k] .* sin.(tr.phase[k]))
end
