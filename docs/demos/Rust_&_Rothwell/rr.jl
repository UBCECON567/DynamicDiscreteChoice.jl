# ---
# title: Reproducing Rust & Rothwell (1995)
# date: 2022
# author: Paul Schrimpf
# ---

using ZipFile, DataFrames, Downloads, CSV, CategoricalArrays

# # Data 
# The data for this paper is available from the *Journal of Applied Econometrics* data archive. 
# It's supplied as a zip archive containing text files. Each text file named nameX.mon contains monthly data from a plant, reactor number X. There's also nukesum.asc with time-variant information about each reactor and plant. 
#
# First, we donwload the zip archive.
dataurl="http://qed.econ.queensu.ca/jae/1995-v10.S/rust-rothwell/rr-data.zip"
datazip="rr-data.zip"
if !isfile(datazip)
    Downloads.download(dataurl, datazip)
end 

# Now, we read the time invariant information.

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
        
# And the time-varying information.

nuket = let 
    zfr = ZipFile.Reader(datazip)
    df = vcat([DataFrame(CSV.File(f, header=false, delim=" ", ignorerepeated=true))
    for f ∈ zfr.files if occursin(r"\.mon$", f.name)]...)
    rename!(df, ["name","ID","steam system","year","month",
    "vintage","age","hours refuel","hours planned outage",
    "hours forced outage","hours total",
    "scram in","scram out","number forced outages"])
    df[!,"steam system"] = recode(df[!,"steam system"], 
    1=>"Babcock Wilcox",2=>"Combustion", 3=>"GE", 4=>"Westinghouse", 
    5=>"Other1",6=>"Other2") 
    df[!,"vintage"] = recode(df[!,"vintage"], 0=>"preTMI", 1=>"postTMI")
    df
end;
    
describe(nuket)
    
