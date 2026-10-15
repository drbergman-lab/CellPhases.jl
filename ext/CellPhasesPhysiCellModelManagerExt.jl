# PhysiCellModelManager support, loaded when PhysiCellModelManager and PhysiCellOutput are present.
#
# A thin adapter: PCMM resolves simulations to their output folders; the PhysiCellOutput extension
# does the rest. Monads and samplings give one PhaseResult per replicate, which every readout takes
# as a vector (see `Replicates`).

module CellPhasesPhysiCellModelManagerExt

using CellPhases
using PhysiCellModelManager: Simulation, Monad, Sampling, simulationIDs, monadIDs
using PhysiCellOutput: PhysiCellSequence

"""
    cellPhases(simulation::Simulation; value, cell_types = nothing, kwargs...) -> PhaseResult
    cellPhases(monad::Monad; kwargs...) -> Vector{PhaseResult}
    cellPhases(sampling::Sampling; kwargs...) -> Dict{Int,Vector{PhaseResult}}

Phases of a cell variable in PhysiCellModelManager simulations. Requires
`using PhysiCellModelManager`.

A `Simulation` is loaded as a `PhysiCellSequence` and every keyword is forwarded to
[`cellPhases(::PhysiCellSequence)`](@ref CellPhases.cellPhases(::PhysiCellOutput.PhysiCellSequence)).
A `Monad` gives one result per replicate simulation; a `Sampling` gives those vectors keyed by
monad ID. Pass the vector to any readout to get per-replicate results and a pooled one.

# Examples
```julia
using PhysiCellModelManager, CellPhases
runs = cellPhases(Monad(1); value = :notch, cell_types = "rgc", band = (100, 300))
coherence = phaseCoherence(runs; dr = 5, rmax = 80)     # Replicates: .results and .pooled
replicateSpread(coherence, :excess)
```
"""
CellPhases.cellPhases(simulation::Simulation; kwargs...) = _phases(simulation.id; kwargs...)

CellPhases.cellPhases(monad::Monad; kwargs...) = [_phases(id; kwargs...) for id in simulationIDs(monad)]

CellPhases.cellPhases(sampling::Sampling; kwargs...) =
    Dict(id => cellPhases(Monad(id); kwargs...) for id in monadIDs(sampling))

_phases(simulation_id::Integer; kwargs...) = cellPhases(PhysiCellSequence(simulation_id); kwargs...)

end
