# Synthetic values with known answers, on positions and IDs the caller supplies (typically a real
# run's): plane waves, pulse trains and an oscillator field. Every generator returns per-frame
# vectors in the shape `cellPhases` takes. Values are positive, so `transform = log` also applies.

"""
    planeWave(positions; dt, period, speed, direction, noise = 0.3, baseline = 1.0,
              rng = Random.default_rng())

Values of a plane wave `exp(baseline + cos(ωt − k·x) + noise·ε)` at the given per-frame positions
(2D): period `period`, travelling at `speed` position units per time unit in `direction`
degrees (0 = +x, 90 = +y). Frame `t` is at time `(t − 1)·dt`; ε is standard normal noise.

The phase [`cellPhases`](@ref) returns is `ωt − k·x`, so φ = 0 at the wave's peaks.

# Examples
```jldoctest
julia> v = planeWave([[(0.0, 0.0)], [(0.0, 0.0)]]; dt = 45, period = 180, speed = 1,
                     direction = 0, noise = 0);

julia> round.(log.(only.(v)); digits = 3)
2-element Vector{Float64}:
 2.0
 1.0
```
"""
function planeWave(positions; dt::Real, period::Real, speed::Real, direction::Real,
                   noise::Real = 0.3, baseline::Real = 1.0, rng::AbstractRNG = Random.default_rng())
    ω = 2π / period
    k = ω / speed .* (cosd(direction), sind(direction))
    return [[exp(baseline + cos(ω * (t - 1) * dt - k[1] * x[1] - k[2] * x[2]) + noise * randn(rng))
             for x in positions[t]] for t in eachindex(positions)]
end

"""
    pulseTrain(ids; dt, period, width = period / 10, jitter = 0.05, baseline = 0.1, noise = 0.2,
               rng = Random.default_rng())

Pulse-like values for the agents in per-frame `ids`: each agent fires a Gaussian pulse of
standard deviation `width` once per cycle, on top of `baseline`, with multiplicative
log-normal noise of SD `noise`. Each agent's period is `period·(1 + jitter·ε)` and its first
pulse is at a random time, so agents are not in step.

A pulse train is broadband: most of its power lies in harmonics of the period, which is what makes
it a real test of the band-pass filter.

# Examples
```jldoctest
julia> using Random

julia> v = pulseTrain([[1], [1], [1]]; dt = 10, period = 170, rng = Xoshiro(1));

julia> length(v), all(>(0), only.(v))
(3, true)
```
"""
function pulseTrain(ids; dt::Real, period::Real, width::Real = period / 10, jitter::Real = 0.05,
                    baseline::Real = 0.1, noise::Real = 0.2, rng::AbstractRNG = Random.default_rng())
    I = eltype(eltype(ids))
    agent = Dict{I,Tuple{Float64,Float64}}()   # id => (period, time of a peak)
    return map(eachindex(ids)) do t
        map(ids[t]) do id
            T, t0 = get!(() -> (T = period * (1 + jitter * randn(rng)); (T, T * rand(rng))), agent, id)
            Δ = mod((t - 1) * dt - t0 + T / 2, T) - T / 2
            (baseline + exp(-Δ^2 / (2width^2))) * exp(noise * randn(rng))
        end
    end
end

"""
    oscillatorField(positions, ids; dt, period, correlationlength = 0.0, phasesd = π,
                    diffusion = 0.0, noise = 0.3, baseline = 1.0, rng = Random.default_rng())

Oscillators sharing one `period` whose phase offsets are spatially correlated:
`exp(baseline + cos(ωt + θᵢ + ψᵢ(t)) + noise·ε)`.

- `θᵢ` is set from the agent's first position. With `correlationlength = 0` it is uniform and
  independent between agents (no coherence at any distance). Otherwise it is a Gaussian random
  field with standard deviation `phasesd` and covariance `exp(−d²/2ℓ²)`, `ℓ = correlationlength`,
  so the expected phase coherence at distance `d` is `exp(−phasesd²(1 − exp(−d²/2ℓ²)))`.
- `ψᵢ(t)` is a random walk with variance `diffusion` rad² per time unit: phase noise that slowly
  breaks the alignment.

# Examples
```jldoctest
julia> using Random

julia> rng = Xoshiro(1); xy = [(50.0rand(rng), 50.0rand(rng)) for _ in 1:20];

julia> positions, ids = fill(xy, 4), fill(1:20, 4);

julia> v = oscillatorField(positions, ids; dt = 10, period = 180,
                           correlationlength = 20.0, phasesd = 1.0, rng = Xoshiro(2));

julia> length(v[4])
20
```
"""
function oscillatorField(positions, ids; dt::Real, period::Real, correlationlength::Real = 0.0,
                         phasesd::Real = π, diffusion::Real = 0.0, noise::Real = 0.3,
                         baseline::Real = 1.0, rng::AbstractRNG = Random.default_rng())
    ω = 2π / period
    K = 400                               # random Fourier features for the Gaussian field
    wave = [SVector(randn(rng), randn(rng)) ./ max(correlationlength, eps()) for _ in 1:K]
    shift = 2π .* rand(rng, K)
    field(x) = phasesd * sqrt(2 / K) * sum(cos(wave[j][1] * x[1] + wave[j][2] * x[2] + shift[j]) for j in 1:K)
    I = eltype(eltype(ids))
    state = Dict{I,Tuple{Float64,Float64,Int}}()   # id => (θ, ψ, last frame)
    return map(eachindex(ids)) do t
        map(eachindex(ids[t])) do i
            id = ids[t][i]
            θ, ψ, last = get!(state, id) do
                (correlationlength > 0 ? field(positions[t][i]) : 2π * rand(rng), 0.0, t)
            end
            diffusion > 0 && t > last && (ψ += sqrt(diffusion * (t - last) * dt) * randn(rng))
            state[id] = (θ, ψ, t)
            exp(baseline + cos(ω * (t - 1) * dt + θ + ψ) + noise * randn(rng))
        end
    end
end
