# PhysiCell output support, loaded when PhysiCellOutput is present.
#
# A thin adapter: PhysiCellOutput's `cellDataSequence` groups each cell's values by ID over time;
# this lays them out per frame on the snapshot time grid and calls the array method. No second
# PhysiCell reader.

module CellPhasesPhysiCellOutputExt

using CellPhases
using PhysiCellOutput: PhysiCellSequence, cellDataSequence
using StaticArrays: SVector

"""
    cellPhases(sequence::PhysiCellSequence; value, cell_types = nothing, include_dead = false,
               dims = 2, kwargs...) -> PhaseResult

Phases of a cell variable across the snapshots of one PhysiCell output folder. Requires
`using PhysiCellOutput`.

- `value`: the cell variable, e.g. `:notch` or any custom data label.
- `cell_types`: cell type names to keep (a name or a collection of names); `nothing` keeps all.
  A cell outside the selection in some snapshots (it changed type) has a gap there, which splits
  its trace.
- `include_dead`: keep cells PhysiCell marks as dead.
- `dims`: use `position_1`, `position_2` (and `position_3` with `dims = 3`).

`dt` is the snapshot interval, read from the snapshot times; snapshots must lie on one time grid
(a missing snapshot is a gap). The remaining keywords (`band`, `order`, `trim`, `minlength`,
`transform`, `mask`) are those of the array method.

# Examples
```julia
using PhysiCellOutput, CellPhases
sequence = PhysiCellSequence("output")
p = cellPhases(sequence; value = :notch, cell_types = "rgc", band = (100, 300))
```
"""
function CellPhases.cellPhases(sequence::PhysiCellSequence; value::Union{Symbol,AbstractString},
                               cell_types = nothing, include_dead::Bool = false,
                               dims::Integer = 2, kwargs...)
    haskey(kwargs, :dt) && throw(ArgumentError("dt is read from the snapshot times; do not pass it"))
    frames = _frames(sequence; value, cell_types, include_dead, dims)
    return cellPhases(frames.positions, frames.values, frames.ids; dt = frames.dt, kwargs...)
end

# Per-frame positions, values and IDs on the snapshot time grid, with the frame interval.
function _frames(sequence::PhysiCellSequence; value, cell_types = nothing,
                 include_dead::Bool = false, dims::Integer = 2)
    1 <= dims <= 3 || throw(ArgumentError("dims must be 1, 2 or 3 for PhysiCell output, got $dims"))
    isempty(sequence.snapshots) && throw(ArgumentError("the sequence has no snapshots"))
    keep = _celltypes(cell_types, sequence)
    times = sort!(unique([s.time for s in sequence.snapshots]))
    length(times) >= 2 || throw(ArgumentError("at least two snapshots are needed for a time step"))
    dt = minimum(diff(times))
    frame(t) = begin
        k = (t - times[1]) / dt
        isapprox(k, round(k); atol = 1e-6) || throw(ArgumentError(
            "snapshot time $t is not on the grid $(times[1]) + k·$dt; snapshots must be equally spaced"))
        round(Int, k) + 1
    end
    nT = frame(times[end])

    label = String(value)
    agents = cellDataSequence(sequence, [label, "position"]; include_dead,
                              include_cell_type_name = true)
    positions = [SVector{dims,Float64}[] for _ in 1:nT]
    values = [Float64[] for _ in 1:nT]
    ids = [Int[] for _ in 1:nT]
    for (id, data) in agents
        v = getproperty(data, Symbol(label))
        pos = data.position
        for (k, t) in enumerate(data.time)
            keep === nothing || data.cell_type_name[k] in keep || continue
            f = frame(t)
            push!(positions[f], SVector{dims,Float64}(ntuple(j -> pos[k, j], dims)))
            push!(values[f], v[k])
            push!(ids[f], id.id)
        end
    end
    # Rows in ID order, so the layout does not depend on dictionary order.
    for f in 1:nT
        o = sortperm(ids[f])
        positions[f], values[f], ids[f] = positions[f][o], values[f][o], ids[f][o]
    end
    return (positions = positions, values = values, ids = ids, dt = Float64(dt))
end

_celltypes(::Nothing, sequence) = nothing
_celltypes(name::Union{AbstractString,Symbol}, sequence) = _celltypes([name], sequence)
function _celltypes(names, sequence)
    keep = Set(String.(names))
    known = Set(values(sequence.cell_type_to_name_dict))
    unknown = setdiff(keep, known)
    isempty(unknown) || throw(ArgumentError(
        "unknown cell type(s) $(join(sort!(collect(unknown)), ", ")); this output has " *
        join(sort!(collect(known)), ", ")))
    return keep
end

end
