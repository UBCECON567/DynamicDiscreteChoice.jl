struct MarkovChain{StateType, M}
  states::OrderedDict{StateType, Int64}
  P::M
  function MarkovChain(states::OrderedDict{StateType,Int64},P) where StateType
    @assert length(states)==size(P,1)
    @assert size(P,1)==size(P,2)
    @assert all(sum(P, dims=2) .≈ one(eltype(P)))
    return(new{StateType,typeof(P)}(states,P))
  end
end

function MarkovChain(states::AbstractVector{StateType}, P) where StateType
  return(MarkovChain(OrderedDict( s => i for (i,s) ∈ enumerate(states)), P))
end

function getindex(c::MarkovChain{StateType,T}, s1::StateType, s2::StateType) where{T, StateType <: Union{Tuple, Symbol, String} }
  return(c.P[c.states[s1], c.states[s2]])
end

"""
   *(c1::MarkovChain, c2::MarkovChain)

Assuming `c1` and `c2` are independent, create a new MarkovChain
representing the evolution of the combination of states in `c1` and
`c2`.
"""
function *(c1::MarkovChain, c2::MarkovChain)
  states = OrderedDict((s1.first, s2.first) => i for
                         (i, (s1, s2)) ∈ enumerate(Iterators.product(c1.states, c2.states)) )
  P = similar(c1.P, length(states), length(states))
  for ((s1old, s2old),iold) ∈ states
    for ((s1new, s2new),inew) ∈ states
      P[inew,iold] = c1[s1new, s1old]*c2[s2new, s2old]
    end
  end
  return(MarkovChain(states,P))
end

function show(io::IO, c::MarkovChain)
  sstates = OrderedDict( "$s" => i for (s,i) ∈ c.states)
  print(io, "hello")
  show(io, NamedArrays.NamedArray(Matrix(c.P),(sstates,sstates), ("old","new")))
end

expectations(x, c::MarkovChain) = c.P*x

import Random: rand
rand(c::MarkovChain, s::Integer) =  rand(Distributions.Categorical(c.P[s,:]))
rand(c::MarkovChain, s::Union{String,Symbol,Tuple}) =  rand(Distributions.Categorical(c.P[c.states[s],:]))
