module DynamicDiscreteChoice

import NamedArrays
import DataStructures: OrderedDict
import Base.getindex, Base.*, Base.show
import Distributions
import NLsolve
import ShiftedArrays: lag
import LinearAlgebra: I

include("markovchain.jl")
export MarkovChain, getindex, *, show, rand,
  emax, value, choicep, estimate


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

"""
  emaxp(p, i0=1, d=Distributions.Gumbel())

Returns E[maxᵢ v[i] - v[i0] + ϵ[i] ] where ϵ[i] are i.i.d. d
and p[i] = P(i = argmaxⱼ v[j] - v[i0] + ϵ[j])
"""
function emaxp(p, i0=1, d::Distributions.Gumbel=Distributions.Gumbel())
  @assert d.θ ≈ 1.0
  @assert d.μ ≈ 0.0
  return(Base.MathConstants.γ - log(p[i0]))
end


"""
    choicespecificvalue!(v, V,ddc)

Given a dynamic discrete choice model, `ddc`, and exante value
function `V`, sets `v` equal to the choice specific value functions.

```math
v[a,s] = u[a,s] + δ E[V(s') | a, s]
```

"""
function choicespecificvalue!(v, V,ddc)
  v .= ddc.payoffs
  for a ∈ axes(v,1)
    v[a,:] .+= ddc.discount*expectations(V,ddc.transition(a))
  end
  return(v)
end


"""
    exante_bellman(V,v, ddc)

Compute the exante (before ϵ is known) Bellman operator for dynamic
discrete choice model `ddc`. i.e. returns
```math
Ṽ[s] = E[maxₐu[a,s] + ϵ[a] + δ E[V[s'] | a, s] ]
```

"""
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

"""
    choicep(v, ddc)

Given choice specific value functions, returns conditional choice probabilities.
"""
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

"""
   simulate(T, ddc, v)

Simulates dynamic discrete choice problem for `T` periods with choice
specific value functions `v`.

Returns a named tuple consisting of states, state indices, actions,
and action indices.
"""

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


"""
    estimate(action_data, state_data, discount; zero_action=action_data[1],
                  Fϵ=Distributions.Gumbel(),
                  actions = OrderedDict(a => i for (i,a) ∈ enumerate(unique(action_data))),
                  states = OrderedDict(s => i for (i,s) ∈ enumerate(unique(state_data))))

Given a vector of observed actions, `action_data`, and states,
`state_data`, and discount rate, computes an estimate of the payoffs
of the associated dynamic discrete choice model.

Normalizes the payoff of `zero_action` to 0 in all states.

`action_data` and `state_data` should be single time-series
realizations of the game.  Conditional sample means are used to
estimate conditional choice probabilities and transition
probabilities. Use `estimate(p, ddc::DDC;
zero_action=first(ddc.actions)[1])` if you want to use other estimates
for choice and transition probabilities.

"""
function estimate(action_data, state_data, discount; zero_action=action_data[1],
                  Fϵ=Distributions.Gumbel(),
                  actions = OrderedDict(a => i for (i,a) ∈ enumerate(unique(action_data))),
                  states = OrderedDict(s => i for (i,s) ∈ enumerate(unique(state_data))))
  payoffs = NamedArrays.NamedArray(zeros(length(actions), length(states)),
                                   (actions, states), ("action", "state"))
  trans = let
    tp = NamedArrays.NamedArray( zeros(length(states), length(states), length(actions)) ,
                                 (states, states, actions), ("new", "old", "action") )
    cnt = NamedArrays.NamedArray( zeros(length(states), length(actions)) ,
                                  (states, actions), ("old", "action") )
    for t ∈ 2:length(state_data)
      tp[state_data[t],state_data[t-1],action_data[t-1]] += 1
      cnt[state_data[t-1],action_data[t-1]] += 1
    end
    for (sold, action) ∈ Iterators.product(keys(states), keys(actions))
      tp[:,sold, action] ./= max(cnt[sold,action],one(eltype(cnt)))
    end
    (a)->MarkovChain(states, tp[:,:,a]')
  end
  ddc = DDC(states, actions, payoffs, discount, trans, Fϵ)
  choicep = deepcopy(payoffs)
  choicep .= zero(eltype(choicep))
  cnts = NamedArrays.NamedArray(zeros(length(states)), (states,),("state",))
  for (a,s) ∈ zip(action_data, state_data)
    choicep[a,s] += 1
    cnts[s] += 1
  end
  for s ∈ keys(states)
    choicep[:,s] ./= max(cnts[s], one(eltype(cnts)))
  end
  pay, v = estimate(choicep, ddc, zero_action=zero_action)
  return(payoffs=pay, v=v, trans=trans, choicep=choicep)
end


"""
    hotzmiller(p, a0, d::Distributions.Gumbel=Distributions.Gumbel())

Given choice probabilities, return differences of choice specific value functions.
"""
function hotzmiller(p, a0, d::Distributions.Gumbel=Distributions.Gumbel())
  @assert d.θ ≈ 1.0
  @assert d.μ ≈ 0.0
  return(log.(p) .- log.(p[a0]))
end

"""
    estimate(p, ddc::DDC; zero_action=first(ddc.actions)[1])

Estimate payoffs of dynamic discrete choice model. `p` should be
estimates of the conditional choice probabilities. `ddc.transition(a)`
should return a `MarkovChain` representing the transitions when action
`a` is chosen.
"""
function estimate(p, ddc::DDC; zero_action=first(ddc.actions)[1])

  payoffs = deepcopy(ddc.payoffs)
  δ = ddc.discount

  v = similar(payoffs)
  # recover v[zero_action, :]
  # v[0, s] = 0 + discount * E[max v[a,s] + ϵ[a] | 0, s]
  #         = 0 + discount * E[v[0,s] | 0, s] +
  #             + discount * E[max v[a,s] - v[0,s] + ϵ[a] | 0, s]
  #         = 0 + discount * E[v[0,s] | 0, s] +
  #             + discount * E[max log(P[a,s]) - log(P[0,s]) + ϵ[a] | 0, s]
  #         = 0 + discount * E[v[0,s] | 0, s] +
  #             + discount * E[ log (sum(exp(log(P[a,s]) - log(P[0,s])) + γ | 0, s]
  em = emaxp.(eachcol(p), zero_action, ddc.Fϵ)
  v[zero_action, :] = (I - δ*ddc.transition(zero_action).P) \ (δ*expectations(em, ddc.transition(zero_action)))
  for s ∈ keys(ddc.states)
    @views v[:,s] .= hotzmiller(p[:,s], zero_action) .+ v[zero_action,s]
  end
  for a ∈ keys(ddc.actions)
    payoffs[a,:] .= v[a,:] .- δ*expectations(emax.(eachcol(v), ddc.Fϵ), ddc.transition(a))
  end
  return(payoffs=payoffs, v=v)
end

end
