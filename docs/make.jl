using DynamicDiscreteChoice
using Documenter, DocumenterCitations, DemoCards

DocMeta.setdocmeta!(DynamicDiscreteChoice, :DocTestSetup, :(using DynamicDiscreteChoice); recursive=true)


# generate demo files
demopage, postprocess_cb, demo_assets = makedemos("demos") # this is the relative path to docs/
# if there are generated css assets, pass it to Documenter.HTML
assets = []
isnothing(demo_assets) || (push!(assets, demo_assets))


bib = CitationBibliography("ddc.bib")

makedocs(bib;
    modules=[DynamicDiscreteChoice],
    authors="Paul Schrimpf <paul.schrimpf@ubc.ca> and contributors",
    repo="https://github.com/UBCECON567/DynamicDiscreteChoice.jl/blob/{commit}{path}#{line}",
    sitename="DynamicDiscreteChoice.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://UBCECON567.github.io/DynamicDiscreteChoice.jl",
        assets=assets,
    ),
    pages=[
        "Home" => "index.md",
        "Method" => "method.md",
        "Example" => "rr-example.md",
        demopage,
    ],
)

# for democards
postprocess_cb()

deploydocs(;
    repo="github.com/UBCECON567/DynamicDiscreteChoice.jl",
    devbranch="main",
    push_preview=true,
)


# add pluto html as in https://github.com/JuliaManifolds/Manopt.jl/blob/00aa4a8f833b4052bfa97ca8bcf7956c3286639c/docs/make.jl#L62-L64
# or https://github.com/fonsp/Pluto.jl/discussions/1345
