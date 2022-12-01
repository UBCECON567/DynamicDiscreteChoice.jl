module DynamicDiscreteChoice

import NamedArrays
import DataStructures: OrderedDict
import Base.getindex, Base.*, Base.show
import Distributions
import NLsolve
import ShiftedArrays: lag
import LinearAlgebra: I
import ProgressMeter: @showprogress
import Statistics: quantile

include("markovchain.jl") 
include("ddc.jl")  # defines types and functions for representing and simulating model
include("estimate.jl")
export MarkovChain, getindex, *, show, rand,
  emax, value, choicep, estimate, ControlledMarkovChain, eltype

end
