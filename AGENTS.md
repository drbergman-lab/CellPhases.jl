# AGENTS.md — CellPhases.jl

Per-agent phase analysis of oscillating values in agent-based model output. The entry point is
`cellPhases` in [`src/phases.jl`](src/phases.jl): per-frame positions, values and IDs in, a
`PhaseResult` (Hilbert phase and amplitude per agent per frame) out. Readouts are constructors
on a `PhaseResult`: `CycleLengths`, `PhaseValueCheck`, `PhaseCoherence`, `OrderParameter`
([`src/readouts.jl`](src/readouts.jl)) and `PhaseWaves` ([`src/waves.jl`](src/waves.jl)). Their
lowercase functions (`phaseCoherence(runs; …)` …) take a vector of results from replicate runs
and return a `Replicates` ([`src/replicates.jl`](src/replicates.jl)).
Plot-ready data is in [`src/plotdata.jl`](src/plotdata.jl), Plots recipes in
[`src/recipes.jl`](src/recipes.jl), Makie recipes in
[`ext/CellPhasesMakieExt.jl`](ext/CellPhasesMakieExt.jl) (functions declared in
[`src/makie.jl`](src/makie.jl)), synthetic values with known answers (on caller-supplied
agents, usually a real run's) in [`src/synthetic.jl`](src/synthetic.jl). Toy agents
(`agentField`) are test-only, in [`test/agent_field.jl`](test/agent_field.jl); docs examples
`include` it.

Start with [`docs/src/dev/architecture.md`](docs/src/dev/architecture.md): module map, the data
flow through one call, invariants. Design decisions are `!!! tierjournal` blocks on the manual
pages, collected into `docs/src/dev/journal.md`. The original brief is
`~/.julia/dev/SpaceTimeCorrelations/sandbox/CellPhases_BRIEF.md`; the prototype it was built from
is in `~/.julia/dev/SpaceTimeCorrelations/sandbox/phase/`.

## Commands

Always run Julia with a project flag; never edit `Manifest.toml` or add a dependency without asking.

```bash
julia --project=. -e 'using Pkg; Pkg.test()'                  # the test suite
julia --project=docs -e 'using Pkg; Pkg.instantiate()'         # once; docs/Project.toml sources ".."
julia --project=docs docs/make.jl                             # build docs into docs/build, runs doctests
julia --project=docs -e 'include("docs/journal.jl")'          # regenerate docs/src/dev/journal.md alone
```

The docs build needs at least one commit (Documenter reads `HEAD`), not a remote.

## Invariants

- **Physical units, one time unit, never assumed.** `band = (pmin, pmax)` are periods and `dt`
  the frame interval, both required, in one time unit the user chooses; `trim`, `minlength`,
  cycle lengths and speeds (position units per time unit) follow it. No unit is named anywhere.
  `bandpassFilter` always passes `fs = 1/dt` to DSP: without it DSP reads cutoffs as fractions
  of Nyquist (this happened once and produced a spurious result). The `filter` testset measures
  the realised response at `pmin`, `pmax` and mid-band for four frame intervals; keep it.
- **Scales have no defaults:** `band`, `dr`, `rmax`, `radii`, `r`.
- **Per-frame fields stay aligned with the input rows.**
- **Sign conventions:** φ = 0 at peaks; direction = angle of −∇φ; a wave moving +y reads +90°.
- **The null is the time-shift surrogate** (`timeShift`), per track. A per-frame shuffle gives a
  far too narrow band; do not use it for bands or p-values.
- **Report the clean-fit share** (`cleanShare`) with any direction or speed summary.
- **Runtime dependencies:** DSP, LinearAlgebra, NearestNeighbors, Random, RecipesBase,
  StaticArrays, Statistics. Makie, PhysiCellOutput and PhysiCellModelManager are `[weakdeps]`
  reached through `ext/`. MAT, Makie and PhysiCellOutput are test-only. MakieCore is not used:
  Makie 0.24 no longer depends on it.
- **Extensions are thin adapters** ending in the array method. No second PhysiCell reader:
  PhysiCellOutput's `cellDataSequence` loads the data.
- **No simulation data in the repo.** PhysiCell test input is generated at test time by
  `test/physicell_output.jl`.
- **PhysiCellModelManager stays out of the test and docs environments.**
- **Plot data holds the decisions; recipes only draw**, in both Plots and Makie. No packaged
  multi-panel figure.
- **Out of scope:** wavelet ridge analysis, cluster/patch tracking, the pooled space-time
  correlation (that is SpaceTimeCorrelations.jl).

## Conventions

- Functions `camelCase` (`cellPhases`), internal helpers `_lowercase`, types `PascalCase`, files
  `snake_case.jl`. A readout of one run is a constructor of its result type
  (`PhaseCoherence(p; …)`); the lowercase name is reserved for the replicate form.
- Every exported name has a docstring with a runnable `jldoctest`. `checkdocs = :all`.
- Docs follow the tiered-docs house style: `!!! tierbrief` / `tierfull` / `tierdev` /
  `tierjournal`, shallowest first within a section, never wrapping a code block.
  `docs/src/dev/journal.md` is generated; edit the blocks on the pages.

## Sharp edges

- Wavelengths shorter than about `2r` are misread by `PhaseWaves`; the clean-fit share collapses
  there. Tested in "plane waves: speed and direction".
- `OrderParameter` under the null is about `1/√n`, not 0; always read it against `null`.
- An agent with a gap has two tracks; surrogates and z-scores treat them independently.
- Doctests run with `using CellPhases, Statistics` (see `docs/make.jl`).
