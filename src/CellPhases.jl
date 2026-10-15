"""
    CellPhases

Per-agent phase analysis of oscillating values in agent-based model output: band-pass each tracked
agent's trace, take its Hilbert phase, and measure how the phases relate in space (coherence with
distance, local order parameter, local phase waves).

Entry point: [`cellPhases`](@ref). Readouts: [`cycleLengths`](@ref), [`phaseValueCheck`](@ref),
[`phaseCoherence`](@ref), [`orderParameter`](@ref), [`phaseWaves`](@ref).
"""
module CellPhases

using DSP: DSP, Bandpass, Butterworth, digitalfilter, filtfilt, hilbert
using LinearAlgebra: norm
using NearestNeighbors: KDTree, inrange
using Random: Random, AbstractRNG, Xoshiro
using RecipesBase: RecipesBase, @recipe, @series
using StaticArrays: SVector
using Statistics: mean, median, quantile, std

export cellPhases, PhaseResult, Track, bandpassFilter, timeShift
export cycleLengths, CycleLengths, phaseValueCheck, PhaseValueCheck
export phaseCoherence, PhaseCoherence, orderParameter, OrderParameter
export phaseWaves, PhaseWaves, cleanShare
export Replicates, replicateSpread
export agentField, planeWave, pulseTrain, oscillatorField
export roseData, coherenceData, orderParameterData, cycleLengthData, phaseValueData,
       portraitData, traceData
export roseplot, roseplot!, coherenceplot, coherenceplot!, orderparameterplot,
       orderparameterplot!, cyclelengthplot, cyclelengthplot!, phasevalueplot, phasevalueplot!,
       portraitplot, portraitplot!

include("phases.jl")
include("surrogates.jl")
include("readouts.jl")
include("waves.jl")
include("replicates.jl")
include("synthetic.jl")
include("plotdata.jl")
include("recipes.jl")
include("makie.jl")

end
