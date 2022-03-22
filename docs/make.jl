using DynamicDiscreteChoice
using Documenter, DocumenterCitations

DocMeta.setdocmeta!(DynamicDiscreteChoice, :DocTestSetup, :(using DynamicDiscreteChoice); recursive=true)


bib = CitationBibliography("ddc.bib")

makedocs(bib;
    modules=[DynamicDiscreteChoice],
    authors="Paul Schrimpf <paul.schrimpf@gmail.com> and contributors",
    repo="https://github.com/ubcecon567/DynamicDiscreteChoice.jl/blob/{commit}{path}#{line}",
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
    repo="github.com/ubcecon567/DynamicDiscreteChoice.jl",
    devbranch="main",
)
