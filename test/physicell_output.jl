# Writes synthetic data in the PhysiCell output format that PhysiCellOutput reads: `initial.xml`
# (cell types and column labels), one `outputNNNNNNNN.xml` per frame (time) and its `_cells.mat`
# (one column per cell, one row per label). Only what the reader needs is written. Same layout as
# SpaceTimeCorrelations' test/physicell_output.jl.

using MAT: matwrite

const PC_TYPES = ["rgc", "layer"]                       # cell type IDs 0 and 1
const PC_LABELS = [("ID", 1), ("position", 3), ("cell_type", 1), ("dead", 1), ("notch", 1)]

function _physicellxml(path, time, matfile)
    types = join(["<type ID=\"$(i - 1)\" type=\"$(i - 1)\">$t</type>" for (i, t) in enumerate(PC_TYPES)])
    k = 0
    labels = join(map(PC_LABELS) do (name, size)
        s = "<label index=\"$k\" size=\"$size\" units=\"none\">$name</label>"
        k += size
        s
    end)
    write(path, """
        <?xml version="1.0"?>
        <MultiCellDS version="2" type="snapshot/simulation">
          <metadata>
            <current_time units="min">$time</current_time>
            <current_runtime units="sec">0.0</current_runtime>
          </metadata>
          <cellular_information><cell_populations><cell_population type="individual"><custom>
            <simplified_data type="matlab" source="PhysiCell">
              <cell_types>$types</cell_types>
              <labels>$labels</labels>
              <filename>$matfile</filename>
            </simplified_data>
          </custom></cell_population></cell_populations></cellular_information>
        </MultiCellDS>
        """)
end

"""
    writePhysiCellOutput(dir, rgc; others = 30, dead = Dict(), layer = Dict(), dt = 10.0)

Write `rgc` (per-frame `positions`, `values`, `ids`) as rgc cells carrying `notch`, plus `others`
fixed cells of type "layer" with constant notch 0.1. `dead[t]` lists rgc IDs to mark dead in frame
`t`; `layer[t]` lists rgc IDs written as type "layer" in frame `t`. Returns `dir`.
"""
function writePhysiCellOutput(dir, rgc; others = 30, dead = Dict{Int,Vector{Int}}(),
                              layer = Dict{Int,Vector{Int}}(), dt = 10.0)
    mkpath(dir)
    _physicellxml(joinpath(dir, "initial.xml"), 0.0, "initial_cells.mat")
    nrow = sum(last, PC_LABELS)
    for t in eachindex(rgc.positions)
        n = length(rgc.positions[t])
        A = zeros(nrow, n + others)
        for i in 1:n
            id = rgc.ids[t][i]
            A[:, i] = [id, rgc.positions[t][i][1], rgc.positions[t][i][2], 0.0,
                       id in get(layer, t, Int[]) ? 1.0 : 0.0,
                       id in get(dead, t, Int[]) ? 1.0 : 0.0, rgc.values[t][i]]
        end
        for j in 1:others   # a fixed row of "layer" cells along the top
            A[:, n + j] = [100_000 + j, 2.0 * j, 500.0, 0.0, 1.0, 0.0, 0.1]
        end
        base = "output" * lpad(t - 1, 8, '0')
        _physicellxml(joinpath(dir, base * ".xml"), (t - 1) * dt, base * "_cells.mat")
        matwrite(joinpath(dir, base * "_cells.mat"), Dict("cells" => A))
    end
    return dir
end
