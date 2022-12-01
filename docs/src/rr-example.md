# Example: [rr1995](@cite)

[rr1995](@cite) use a dynamic discrete choice model to analyze how Three Mile Island and the associated regulatory changes affected nuclear power plant operation. 

## Load Data

```@example
using ZipFile, DataFrames, Downloads, CSV, CategoricalArrays
dataurl="http://qed.econ.queensu.ca/jae/1995-v10.S/rust-rothwell/rr-data.zip"
datazip="rr-data.zip"
if !isfile(datazip)
    Downloads.download(dataurl, datazip)
end

nukesum = let 
	zfr = ZipFile.Reader(datazip) 
	df = filter(x->x.name=="nukesum.asc",zfr.files)[1] |>
		x->CSV.File(x, header=false, delim=" ", ignorerepeated=true) |>
		DataFrame 
	rename!(df, ["name","firm",
	"start month","start day","start year",
    "thermal capacity","nameplate rating",
    "net power","arch engineer","steam system",
    "builder","turbine","state"])
	df
end

describe(nukesum)
```

We can write documentation with code examples in Documenter style md files, but if we want to include lots of code, it might be more convenient to instead use
[Literate.jl](https://github.com/fredrikekre/Literate.jl), perhaps in combination with [DemoCards.jl](https://github.com/JuliaDocs/DemoCards.jl). 

It's possible to combine other literate programming formats with documenter. I used [Weave.jl](https://github.com/JunoLab/Weave.jl) with Documenter in [GMMInference.jl](https://github.com/schrimpf/GMMInference.jl). 

[Pluto.jl can also be integrated with Documenter](https://github.com/fonsp/Pluto.jl/discussions/1345).

