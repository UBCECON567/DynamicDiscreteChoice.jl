using DynamicDiscreteChoice
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(DynamicDiscreteChoice, :DocTestSetup, :(using DynamicDiscreteChoice); recursive=true)


bib = CitationBibliography("ddc.bib")

makedocs(bib;
    modules=[DynamicDiscreteChoice],
    authors="Paul Schrimpf <paul.schrimpf@ubc.ca> and contributors",
    repo="https://github.com/UBCECON567/DynamicDiscreteChoice.jl/blob/{commit}{path}#{line}",
    sitename="DynamicDiscreteChoice.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ubcecon567.github.io/DynamicDiscreteChoice.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/UBCECON567/DynamicDiscreteChoice.jl",
    devbranch="main",
)


# add pluto html as in https://github.com/JuliaManifolds/Manopt.jl/blob/00aa4a8f833b4052bfa97ca8bcf7956c3286639c/docs/make.jl#L62-L64
# or https://github.com/fonsp/Pluto.jl/discussions/1345
