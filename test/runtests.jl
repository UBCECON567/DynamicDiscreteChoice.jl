using .DynamicDiscreteChoice

using Test, SparseArrays, Statistics, NamedArrays
import DataStructures: OrderedDict
import Distributions

@testset "MarkovChain" begin

  c1 = DynamicDiscreteChoice.MarkovChain( ["lo","hi"], sparse([0.9 0.1;  0.2 0.8]) )
  c2 = DynamicDiscreteChoice.MarkovChain( [:a,:b, :c], sparse([0.9 0.1 0.0;
                                                               0.1 0.7  0.2;
                                                               0.0 0.4  0.6] ))

  c3=c1*c2
  @test issparse(c3.P)

end


@testset "emax" begin
  dists = [Distributions.Gumbel(), Distributions.Gumbel(2.5, 1.0),
           Distributions.Gumbel(0.5, 0.5), Distributions.Gumbel(0, 2.0)]
  vs = [zeros(4), [1., 10., 100.], rand(10)]
  S = 100_000
  simemax(v,d) = mean(maximum(v + rand(d,size(v))) for s in 1:S)
  for (d, v) ∈ Iterators.product(dists, vs)
    @test abs(DynamicDiscreteChoice.emax(v,d) - simemax(v,d))/sqrt((d.θ^2*π^2/6)/S) < 3
  end
end


@testset "DDC" begin
  actions = ["in","out"]
  endostates = ["inlast","outlast"]
  exostates = ["lo","hi"]
  exochain = DynamicDiscreteChoice.MarkovChain(exostates, [0.8 0.2; 0.3 0.7])
  trans(a::Int) = trans(actions[a])
  function trans(a)
    if (a=="in")
      endochain = DynamicDiscreteChoice.MarkovChain(endostates, [1.0 0.0; 1.0 0.0])
    else
      endochain = DynamicDiscreteChoice.MarkovChain(endostates, [0.0 1.0; 0.0 1.0])
    end
    return(endochain*exochain)
  end
  states = trans("in").states
  actiondict = OrderedDict(a => i for (i,a) ∈ enumerate(actions))
  payoffs = NamedArray(zeros(2,4), (actiondict,states), ("action","state"))
  payoffs["out",:] .= 0.0
  payoffs["in",("inlast","lo")] .= -0.1
  payoffs["in",("inlast","hi")] .= 0.2
  payoffs["in",("outlast","lo")] .= -0.5
  payoffs["in",("outlast","hi")] .= 0.1
  discount = 0.5
  Fϵ = Distributions.Gumbel()
  ddc = DynamicDiscreteChoice.DDC(states, actiondict, payoffs, discount, trans, Fϵ)

  res = DynamicDiscreteChoice.value(ddc, show_trace=true, method=:anderson, m=0)
  p = DynamicDiscreteChoice.choicep(res.v, ddc)

  @testset "value" begin
    # simulated average value should approximately equal the compute res.V
    v = res.v
    T = 10_000
    pay = zeros(T)
    states = Vector{Int}(undef, T)
    let s = 1
      for t ∈ eachindex(pay)
        ϵ = rand(Distributions.Gumbel(), length(ddc.actions))
        a = argmax(v[:,s] + ϵ)
        pay[t] = ddc.payoffs[a,s] + ϵ[a]
        states[t] = s
        s = rand(ddc.transition(a),s)
      end
    end
    simV = [
      mean( sum(pay[t:end].*ddc.discount.^(0:(T-t))) for t ∈ findall(states.==s) )
    for s in 1:length(ddc.states) ]
    # these should be close, but need very large T to be really close
    display(hcat(simV, res.V))

    sim = DynamicDiscreteChoice.simulate(T, ddc, v)
    for s ∈ keys(ddc.states)
      for a ∈ keys(ddc.actions)
        @show p[a,s], mean(sim.actions[findall(==(s),sim.states)].==a)
        @test isapprox(p[a,s], mean(sim.actions[findall(==(s),sim.states)].==a),
                       atol=3*sqrt(p[a,s]*(1-p[a,s])/sum(==(s), sim.states)))
      end
    end
  end
  @testset "estimate" begin
    T = 20_000
    sim = DynamicDiscreteChoice.simulate(T, ddc, res.v)
    est = DynamicDiscreteChoice.estimate(sim.actions, sim.states, ddc.discount, zero_action="out", states=ddc.states, actions=ddc.actions)
  end
end
