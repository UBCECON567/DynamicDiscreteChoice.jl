"""
    estimate_series(action_data, state_data, discount; zero_action=action_data[1],
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
function estimate_series(action_data, state_data, discount; zero_action=action_data[1],
                        Fϵ=Distributions.Gumbel(),
                        actions = OrderedDict(a => i for (i,a) ∈ enumerate(unique(action_data))),
                        states = OrderedDict(s => i for (i,s) ∈ enumerate(unique(state_data))))

    zero_action ∈ keys(actions) || error("Action \"$zero_action\" not in keys(actions)=$(keys(actions))")    
    payoffs = NamedArrays.NamedArray(zeros(length(actions), length(states)),(actions, states), ("action", "state"))
    trans, tp = let
        tp = NamedArrays.NamedArray( zeros(length(states), length(states), length(actions)) ,
            (states, states, actions), ("new", "old", "action") )
            cnt = NamedArrays.NamedArray( zeros(length(states), length(actions)) ,
            (states, actions), ("old", "action") )
        for t ∈ 2:length(state_data)
            if !(ismissing(state_data[t]) || ismissing(state_data[t-1]) || ismissing(action_data[t-1]))
                tp[state_data[t],state_data[t-1],action_data[t-1]] += 1
                cnt[state_data[t-1],action_data[t-1]] += 1
            end
        end
        for (sold, action) ∈ Iterators.product(keys(states), keys(actions))
            tp[:,sold, action] ./= max(cnt[sold,action],one(eltype(cnt)))
        end
        trans = ControlledMarkovChain(OrderedDict(a=>MarkovChain(states,tp[:,:,a]') for a ∈ keys(actions)))
        trans, tp
    end
    ddc = DDC(states, actions, payoffs, discount, trans, Fϵ)
    choicep = deepcopy(payoffs)
    choicep .= zero(eltype(choicep))
    cnts = NamedArrays.NamedArray(zeros(length(states)), (states,),("state",))
    for (a,s) ∈ zip(action_data, state_data)
        if !ismissing(a) && !ismissing(s)
            choicep[a,s] += 1
            cnts[s] += 1
        end
    end
    for s ∈ keys(states)
        choicep[:,s] ./= max(cnts[s], one(eltype(cnts)))
    end
    pay, v = estimate(choicep, trans, ddc, zero_action=zero_action)
    return(payoffs=pay, v=v, choicep=choicep, transarray=tp, ddc=ddc)
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
    estimate(p, transition, ddc::DDC; zero_action=first(ddc.actions)[1])
    
Estimate payoffs of dynamic discrete choice model. `p` should be
estimates of the conditional choice probabilities. `transition`
should be a "states" by "states" by "actions" array, with `transition[new,old,action]=P(s=new|s=old,a=action)`.
"""
function estimate(p, tp::AbstractArray, ddc::DDC; zero_action=first(ddc.actions)[1])
    trans = ControlledMarkovChain(OrderedDict(a=>MarkovChain(ddc.states,tp[:,:,a]') for a ∈ keys(ddc.actions)))
    estimate(p, trans, ddc, zero_action=zero_action)
end 
    
    
"""
    estimate(p, transition, ddc::DDC; zero_action=first(ddc.actions)[1])
    
Estimate payoffs of dynamic discrete choice model. `p` should be
estimates of the conditional choice probabilities. `transition`
should be a `ControlledMarkovChain` representing the transitions when action
`a` is chosen.
"""
function estimate(p, transition::ControlledMarkovChain, ddc::DDC; zero_action=first(ddc.actions)[1])
    payoffs = similar(p, promote_type(eltype(p), eltype(transition)), size(p))
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
    v[zero_action, :] = (I - δ*transition(zero_action).P) \ (δ*expectations(em, transition(zero_action)))
    for s ∈ keys(ddc.states)
        @views v[:,s] .= hotzmiller(p[:,s], zero_action) .+ v[zero_action,s]
    end
    for a ∈ keys(ddc.actions)
        payoffs[a,:] .= v[a,:] .- δ*expectations(emax.(eachcol(v), ddc.Fϵ), transition(a))
    end
    return(payoffs=payoffs, v=v)
end

"""
    bootstrap_series(payoffs,v,T,ddc; B=999)

Use the parametric bootstrap for inference for a dynamic discrete choice model. 
Data is assumed to come from a single time series of length `T`.

# Inputs:
- `payoffs` estimated payoffs
- `v` estimated choice specific value functions
- `T` length of data 
- `ddc` DDC model containing estimated transitions in `ddc.transition`
- `B` number of bootstrap replications
"""
function bootstrap_series(payoffs,v,T, ddc; B=999)
    bpayoffs = Vector{typeof(payoffs)}(undef,B)
    @showprogress "Bootstrapping ..." for b ∈ eachindex(bpayoffs)
        state_data, _, action_data, _ = simulate(T,ddc,v)
        bpayoffs[b], _ = estimate_series(action_data, state_data, ddc.discount; zero_action=first(ddc.actions)[1],
                                        Fϵ=ddc.Fϵ,
                                        actions = ddc.actions,
                                        states = ddc.states)
    end
    return(bpayoffs)
end


"""
    bootstrap_table(θ,θb;coverage=0.95)

Given a matrix of parameter estimates, θ, and a vector of matrices of 
bootstrap replicants, θb, returns an array with rows alternating between
 rows of θ and rows of tuples containing confidence intervals with 
 `coverage` coverage probability.
"""
function bootstrap_table(θ,θb;coverage=0.95)
    out = similar(θ, Union{eltype(θ), Tuple{eltype(θ), eltype(θ)}}, (size(θ,1)*2, size(θ,2)))
    bootstrap_table!(out,θ,θb; coverage=coverage)
    return(out)
end 

function bootstrap_table!(out, θ,θb;coverage=0.95)
    for i ∈ eachindex(θ)        
        bs = [b[i] - θ[i] for b ∈ θb]
        ci = quantile(bs, [(1-coverage)/2, 1-(1-coverage)/2])
        out[2i-1] = θ[i]
        out[2i] = (ci[1], ci[2]) .+ θ[i]
    end
    return(out)
end

function bootstrap_table(θ::NamedArrays.NamedArray,θb;coverage=0.95)
    out = similar(θ, Union{eltype(θ), Tuple{eltype(θ), eltype(θ)}}, (size(θ,1)*2, size(θ,2)))
    bootstrap_table!(out,θ,θb; coverage=coverage)
    NamedArrays.setnames!(out, vcat([[n,"$n $(Int(round(100*coverage)))% CI"] for n in names(θ)[1]]...), 1)
    return(out)
end
