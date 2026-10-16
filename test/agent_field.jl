# Toy agents for tests and docs examples: positions and IDs with no real run behind them. Not
# package API; load with
#     include(joinpath(pkgdir(CellPhases), "test", "agent_field.jl"))

using Random: Random, AbstractRNG

"""
    agentField(; n, L, nframes, step = 0.0, rng = Random.default_rng()) -> (; positions, ids)

`n` tracked agents placed uniformly at random in the square `[0, L]²`, present in all `nframes`
frames with IDs `1:n`. With `step > 0` each agent takes a Gaussian random-walk step of that
standard deviation per frame, reflected at the walls.
"""
function agentField(; n::Integer, L::Real, nframes::Integer, step::Real = 0.0,
                    rng::AbstractRNG = Random.default_rng())
    x = [(L * rand(rng), L * rand(rng)) for _ in 1:n]
    positions = Vector{Vector{NTuple{2,Float64}}}(undef, nframes)
    for t in 1:nframes
        positions[t] = copy(x)
        step > 0 || continue
        x = [_reflect.(xi .+ step .* (randn(rng), randn(rng)), L) for xi in x]
    end
    return (positions = positions, ids = [collect(1:n) for _ in 1:nframes])
end

_reflect(c, L) = (c = mod(c, 2L); c > L ? 2L - c : c)
