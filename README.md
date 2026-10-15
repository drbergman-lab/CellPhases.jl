# CellPhases.jl

A phase for every tracked agent's oscillating value in agent-based model output, and how those
phases relate in space: coherence with distance, local synchrony, and local phase waves with
their direction and speed.

```julia
using CellPhases

# per-frame positions, values and IDs, frames dt apart
p = cellPhases(positions, values, ids; dt = 10, band = (100, 300))   # band: periods, same time unit as dt

cycleLengths(p)                                  # time between successive peaks
phaseCoherence(p; dr = 5, rmax = 80)             # ⟨cos Δφ⟩ by distance, with a time-shift null
orderParameter(p; radii = [10, 20, 40])          # local synchrony R, with its null
w = phaseWaves(p; r = 30)                        # local wave direction and speed
roseData(w)                                      # plot-ready rose diagram
```

With `using PhysiCellOutput` it reads PhysiCell output directly, and with
`using PhysiCellModelManager` it takes simulations and monads (one result per replicate; every
readout accepts the vector):

```julia
p = cellPhases(PhysiCellSequence("output"); value = :notch, cell_types = "rgc", band = (100, 300))
```

`plot(result)` draws each readout with Plots.jl; with a Makie backend loaded there are
`roseplot`, `coherenceplot` and friends, and Makie's `plot(result)` works too. Synthetic data with known answers
(`planeWave`, `pulseTrain`, `oscillatorField`) is exported for trying the methods.

Build the docs locally with `julia --project=docs docs/make.jl`.
