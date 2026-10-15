using CellPhases
using Documenter

# Loaded so the PhysiCellOutput extension is present and its docstring is rendered and checked.
# The PhysiCellModelManager extension is left out of the build to keep PCMM out of the docs
# environment; its docstring lives in ext/ for the REPL.
using PhysiCellOutput
const CellPhasesPhysiCellOutputExt = Base.get_extension(CellPhases, :CellPhasesPhysiCellOutputExt)

ENV["GKSwstype"] = "100"   # GR without a display, for the plotting page

# Regenerate docs/src/dev/journal.md from the `!!! tierjournal` blocks on the pages,
# and check that every page's tier blocks run shallowest first.
include("journal.jl")

DocMeta.setdocmeta!(CellPhases, :DocTestSetup, :(using CellPhases, Statistics); recursive = true)

makedocs(;
    modules = [CellPhases, CellPhasesPhysiCellOutputExt],
    authors = "Daniel Bergman <danielrbergman@gmail.com> and contributors",
    sitename = "CellPhases.jl",
    # Named explicitly so the docs build before the repository has a remote.
    repo = Remotes.GitHub("drbergman-lab", "CellPhases.jl"),
    format = Documenter.HTML(;
        canonical = "https://drbergman-lab.github.io/CellPhases.jl",
        edit_link = "main",
        assets = ["assets/tiers.css", "assets/tiers.js"],
        size_threshold = 400 * 2^10,
    ),
    checkdocs = :all,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "man/phases.md",
            "man/checks.md",
            "man/synchrony.md",
            "man/waves.md",
            "man/physicell.md",
            "man/plotting.md",
            "man/synthetic.md",
            "man/pitfalls.md",
        ],
        "API reference" => "reference.md",
        "Developers" => [
            "Architecture" => "dev/architecture.md",
            "Journal" => "dev/journal.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/drbergman-lab/CellPhases.jl",
    devbranch = "main",
)
