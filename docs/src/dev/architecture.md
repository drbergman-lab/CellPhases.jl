# Architecture

## Module map

| File | Holds |
|:--|:--|
| `src/CellPhases.jl` | Module, imports, exports, includes. |
| `src/phases.jl` | Steps 1–3: `Track`, `PhaseResult`, `bandpassFilter`, `cellPhases` (traces, filter, Hilbert phase, trim, mask). |
| `src/time_shift.jl` | `timeShift` and `_shifted`, the time-shift null. |
| `src/readouts.jl` | `CycleLengths`, `PhaseValueCheck`, `PhaseCoherence`, `OrderParameter`, their result types, and the null statistics. |
| `src/waves.jl` | `PhaseWaves`, `PhaseWaves`, `cleanShare`, the local plane fit. |
| `src/replicates.jl` | `Replicates`, `replicateSpread`, and the readout methods for vectors of results with their pooling. |
| `src/synthetic.jl` | `planeWave`, `pulseTrain`, `oscillatorField`: values with known answers on given agents. |
| `src/plotdata.jl` | `roseData` and the other plot-data functions. |
| `src/recipes.jl` | RecipesBase recipes, drawn from the plot-data functions. |
| `src/makie.jl` | Declarations and docstrings of the Makie recipe functions (`roseplot` …). |
| `ext/CellPhasesMakieExt.jl` | The Makie recipes, their preferred axes, and `plottype` for each result. |
| `ext/CellPhasesPhysiCellOutputExt.jl` | The `PhysiCellSequence` method: `cellDataSequence`, per-frame layout, then the array method. |
| `ext/CellPhasesPhysiCellModelManagerExt.jl` | `Simulation`, `Monad`, `Sampling` methods: resolve to output folders and forward. |
| `test/runtests.jl` | Filter, tracks, sign conventions, plane-wave validation, broadband test, nulls, replicates, plot data and recipes. |
| `test/agent_field.jl` | `agentField`, toy agents for tests and docs examples (not API). |
| `test/makie.jl` | The Makie extension, with Makie alone (no backend). |
| `test/physicell.jl`, `test/physicell_output.jl` | The PhysiCellOutput extension, on synthetic data written in the PhysiCell output format at test time. |

## One call, start to finish

`cellPhases(positions, values, ids; dt, band)` then `PhaseWaves(p; r)`:

1. `cellPhases` checks lengths, builds the filter (`bandpassFilter`, always with `fs = 1/dt`),
   and derives `trim` (`_filtermemory`) and `minlength` if not given.
2. Positions become `SVector{D,Float64}`; `_byagent` gathers each ID's `(frames, rows)`;
   `_runs` splits at gaps.
3. Each long-enough run is transformed (`_transform`), mean-subtracted, `filtfilt`-ed and
   `hilbert`-ed. Interior frames receive phase and amplitude, written back at `[frame][row]`; the
   run becomes a `Track`.
4. `PhaseWaves` maps each `(frame, row)` to the agent's row in the next frame through the tracks,
   builds a `KDTree` per frame, fits a plane per agent (`_planefit`) and pushes one entry per fit.

`PhaseCoherence` and `OrderParameter` draw all surrogate phase sets first (`_phasesets`) and
accumulate observed and surrogate sums in one pass over neighbour pairs; `_nullstats` turns the
surrogate columns into null, band and p-value.

## Invariants

- **Per-frame fields stay aligned with the input.** `phase[t][i]` is the agent in row `i` of
  input frame `t`; nothing reorders rows. The PhysiCell extension sorts rows by ID once, before
  the array method.
- **One time unit, the user's.** Band, `trim`, `minlength` and `dt` share whatever unit the
  input uses; no unit is assumed. Speeds are position units per that unit. `fs` is always
  passed to DSP.
- **Scales have no defaults:** `band`, `dr`, `rmax`, `radii`, `r`.
- **Sign conventions:** φ = 0 at peaks of the band-passed signal; direction is the angle of −∇φ;
  a wave moving +y reads +90°.
- **Tracks, not IDs, are the unit of a series.** Surrogates shift tracks; cycle lengths and
  z-scores are per track.
- **The null is the time-shift surrogate.** Never a per-frame shuffle for bands or p-values.
- **Extensions are thin adapters** ending in the array method; no second PhysiCell reader.
- **Plot data holds the decisions; recipes only draw.**

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Each area is its own `@testset` in `test/runtests.jl`; the shared geometry `FIELD` is defined at
the top, so a single testset can be pasted into a REPL after those few lines. The PhysiCell test
writes its input to a temporary directory with MAT.jl.

## Adding to it

A new readout goes in `src/readouts.jl` (or its own file) with a result struct and a `show`
method, a method for `AbstractVector{<:PhaseResult}` in `src/replicates.jl` that defines how it
pools, a plot-data function in `src/plotdata.jl` if it is drawn, and a recipe that draws from
that function. A new input form should end in a call to the array method of `cellPhases`.
