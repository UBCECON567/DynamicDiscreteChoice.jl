module DynamicDiscreteChoice

import NamedArrays
import DataStructures: OrderedDict
import Base.getindex, Base.*, Base.show
import Distributions
import NLsolve

include("markovchain.jl")
export MarkovChain, getindex, *, show, rand,
  emax, value, choicep


struct DDC{StateType, ActionType, PayType, T, TranType, Dist}
  states::OrderedDict{StateType,Int64}
  actions::OrderedDict{ActionType,Int64}
  payoffs::PayType
  discount::T
  transition::TranType
  Fϵ::Dist
end


"""
  emax(v, d=Distributions.Gumbel())

Returns E[maxᵢ v[i] + ϵ[i] ] where ϵ[i] are i.i.d. d
"""
function emax(v, d::Distributions.Gumbel=Distributions.Gumbel())
  maxv = maximum(v)
  return(maxv + d.μ + d.θ*log(sum(v->exp((v.-maxv)/d.θ),v)) + d.θ*Base.MathConstants.γ)
end



function choicespecificvalue!(v, V,ddc)
  v .= ddc.payoffs
  for a ∈ axes(v,1)
    v[a,:] .+= ddc.discount*expectations(V,ddc.transition(a))
  end
  return(v)
end

function exante_bellman(V,v, ddc)
  choicespecificvalue!(v, V,ddc)
  out = deepcopy(V)
  for i ∈ eachindex(V)
    out[i] = emax(v[:,i], ddc.Fϵ)
  end
  return(out)
end

"""
   value(ddc; kwargs...)

Compute the value function for the dynamic discrete choice model `ddc`.

`NLsolve.fixedpoint` is used to compute the value functions. `kwargs` are
options passed to `NLsolve.fixedpoint`.

Returns a name tuple containing the value function at each state, `V`,
the choice specific value function, `v`, and the output of
`NLsolve.fixedpoint`, `solver_output`.
"""
function value(ddc; kwargs...)
  V = zeros(length(ddc.states))
  v = deepcopy(ddc.payoffs)
  res = NLsolve.fixedpoint(V->exante_bellman(V,v,ddc), V;
                           kwargs...)
  V = res.zero
  choicespecificvalue!(v,V,ddc)
  return(V=V, v=v, solver_output=res)
end

function choicep(v, ddc)
  @assert ddc.Fϵ == Distributions.Gumbel()
  p = deepcopy(v)
  for s ∈ axes(v,2)
    vmax = maximum(v[:,s])
    p[:,s] .= exp.(v[:,s].-vmax)./sum(exp.(v[:,s].-vmax))
  end
  return(p)
end

"""
   simulate(T, ddc)

Simulates dynamic discrete choice problem for `T` periods.

Returns a named tuple consisting of states, state indices, actions, and action indices.
"""
function simulate(T, ddc)
  V, v, _ = value(ddc)
  return(simulate(T, ddc, v))
end

function simulate(T, ddc, v)
  state, si = rand(ddc.states)
  states = Vector{typeof(state)}(undef, T)
  idxs = Vector{typeof(si)}(undef, T)
  action, ai = rand(ddc.actions)
  actions = Vector{typeof(action)}(undef, T)
  idxa = Vector{typeof(ai)}(undef, T)
  for t ∈ 1:T
    idxs[t] = si
    ϵ = rand(Distributions.Gumbel(), length(ddc.actions))
    ai = argmax(v[:,si] + ϵ)
    idxa[t] = ai
    si = rand(ddc.transition(ai),si)
  end
  for a ∈ ddc.actions
    actions[idxa.==a[2]] .= a[1]
  end
  for s ∈ ddc.states
    for i in eachindex(states)
      if idxs[i]==s[2]
        states[i] = s[1]
      end
    end
  end

  return(states=states, state_i = idxs, actions=actions, action_i = idxa)
end

end
