module DynamicDiscreteChoice

import POMDPs, Distributions

ActionType=Symbol
StateType = Tuple{:Symbol, :Symbol, :Int64, :Int64, Float64}
struct NuclearMDP{Td,T} <: POMDPs.MDP{StateType, ActionType}
  # fields should be the parameters of the problem
  discount::Td
  ϕc::T
  ϕr::T
  ϕrf::T
  ϕd::T
  ϕa::Vector{T}
  ϕu0f::T
  ϕu100f::T
  ϕseason::Vector{T}
  maxduration::Int64
  states::Vector{StateType}
  stateindices::Dict{StateType,Int64}
  actions::Vector{ActionType}
  actionindices::Dict{ActionType,Int64}
end

function NuclearMDP(discount, ϕc, ϕr, ϕd, ϕa, ϕu0f, ϕu100f, ϕseason, maxduration, nepsilon)
  previous_spell = [:problem, :refuelling, :operating]
  signal = [:none, :force_outage, :problem]
  duration = 0:dp.maxduration
  season = 1:length(dp.ϕseason)
  ϵgrid = (range(0,1,nepsilon+1).-1/(2*nepsilon))[2:end]
  states = [(r,f,d,m,ϵ) for (r,f,d,m) in
              product(previous_spell, signal, duration, season,ϵgrid)]
  stateindices = Dict(s=>i for (i,s) in enumerate(states))
  g=range(1,100,length(dp.ϕa)-2)
  actions = [:close, :refuel, :run0,
             (Symbol(("run{g[i]}-{g[i+1]}")) for i in 1:(length(g)-1))..., :run100]
  actionindices = Dict(a=>i for (i,a) in enumerate(actions))
  return(NuclearMDP(discount, ϕc, ϕr, ϕd, ϕa, ϕu0f, ϕu100f, ϕseason,
                    states, stateindices, actions, actionindices))
end

POMDPs.stateindex(dp::NuclearMDP, s::StateType) = dp.stateindices[s]
POMDPs.states(dp::NuclearMDP) = dp.states

POMDPs.actionindex(dp::NuclearMDP, a::ActionType) = dp.actionindices[a]
function POMDPs.actions(dp::NuclearMDP, s::StateType)
  r,signal,_,_,_ = s
  if r==:operating
    if signal!=:problem
      return(dp.actions)
    else #if signal==:problem
      return([:close, :run0]) # enter major problem spell
    end
  elseif r==:refuelling
    if signal==:problem
      return([:close,:refuel])
    else
      return(dp.actions)
    end
  else # if r==:problem
    if signal==:problem
      return([:close, :run0])
    else
      return(dp.actions)
    end
  end
end

function POMDPs.transition(dp::NuclearMDP, s, a)
  if a==:
end



end
