# PhysiCell and replicates

With PhysiCellOutput loaded, `cellPhases` takes a whole output folder; with PhysiCellModelManager
loaded, it takes simulations, monads and samplings.

!!! tierbrief
    Name the cell variable with `value` and pick cell types with `cell_types`. `dt` is read from
    the snapshot times.

!!! tierfull
    Restricting to the cell types that carry the signal matters: a variable that is constant in
    most cell types adds agents with no rhythm. A cell that changes into an unselected type and
    back has a gap there, and its trace is split. Dead cells are dropped unless
    `include_dead = true`.

```julia
using PhysiCellOutput, CellPhases

sequence = PhysiCellSequence("output")
p = cellPhases(sequence; value = :notch, cell_types = "rgc", band = (100, 300))
```

!!! tierdev
    `ext/CellPhasesPhysiCellOutputExt.jl` calls PhysiCellOutput's `cellDataSequence`, which
    groups each cell's values by ID over time, lays the rows out per frame on the snapshot time
    grid (`_frames`, rows in ID order), and calls the array method. Snapshot times must lie on
    one grid; a missing snapshot is an empty frame. The extension test writes PhysiCell-format
    output at test time (`test/physicell_output.jl`) and checks the result against the array
    method on the same data.

!!! tierjournal "2026-10-15 — PhysiCell input through cellDataSequence"
    Decided: the extension reads through PhysiCellOutput's `cellDataSequence` with
    `include_cell_type_name = true` and filters cell types itself.
    Why: `cellDataSequence` already groups each cell's values by ID across snapshots and drops
    dead cells; a second reader would duplicate it.

## Replicates and conditions

!!! tierbrief
    Every readout also takes a vector of `PhaseResult`s and returns a [`Replicates`](@ref): one
    result per replicate in `results` and a `pooled` one.

!!! tierfull
    One simulation's coherence has a null band but no measure of how much it varies between
    runs. To compare conditions, for example wild type against a knockout, run each condition's
    replicates and compare the between-replicate spread (`replicateSpread`) rather than the
    within-run null. Pooling concatenates cycle lengths, phase-value pairs and wave fits, and
    averages coherence and the order parameter weighted by each replicate's pair counts.

```julia
using PhysiCellModelManager, CellPhases

wt = cellPhases(Monad(1); value = :notch, cell_types = "rgc", band = (100, 300))
ko = cellPhases(Monad(2); value = :notch, cell_types = "rgc", band = (100, 300))
cwt = phaseCoherence(wt; dr = 5, rmax = 80, frames = 1:2:600)
cko = phaseCoherence(ko; dr = 5, rmax = 80, frames = 1:2:600)
replicateSpread(cwt, :excess), replicateSpread(cko, :excess)
```

!!! tierdev
    `ext/CellPhasesPhysiCellModelManagerExt.jl` resolves a `Simulation` by ID to a
    `PhysiCellSequence` and forwards; `Monad` maps over `simulationIDs`, `Sampling` over
    `monadIDs`. It is not loaded in the test environment (PCMM is heavy and needs a database);
    the code it forwards to is tested.
