using DynamicDiscreteChoice

using Test, SparseArrays, Statistics, NamedArrays
import DataStructures: OrderedDict
import Random, Distributions, ForwardDiff

@testset "MarkovChain" begin

  c1 = DynamicDiscreteChoice.MarkovChain( ["lo","hi"], sparse([0.9 0.1;  0.2 0.8]) )
  c2 = DynamicDiscreteChoice.MarkovChain( [:a,:b, :c], sparse([0.9 0.1 0.0;
                                                               0.1 0.7  0.2;
                                                               0.0 0.4  0.6] ))
  c3=c1*c2
  @test issparse(c3.P)
end

@testset "ControlledMarkovChain" begin
  c1 = MarkovChain( ["lo","hi"], sparse([0.9 0.1;  0.2 0.8]) )
  endostates = ["outlast","inlast"]
  trans = ControlledMarkovChain( OrderedDict(
    "out"=> DynamicDiscreteChoice.MarkovChain(endostates, [0.0 1.0; 0.0 1.0])*c1,
    "in" => DynamicDiscreteChoice.MarkovChain(endostates, [1.0 0.0; 1.0 0.0])*c1
  ))
  @test isa(trans("in"),MarkovChain)
  @test isa(trans(1),MarkovChain)
  @test trans(2) === trans("in")
  @test trans(1) === trans("out")
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
  actions = ["out","in"]
  endostates = ["outlast","inlast"]
  exostates = ["lo","hi"]
  exochain = DynamicDiscreteChoice.MarkovChain(exostates, [0.8 0.2; 0.3 0.7])
  trans = ControlledMarkovChain( OrderedDict(
    "out"=> DynamicDiscreteChoice.MarkovChain(endostates, [0.0 1.0; 0.0 1.0])*exochain,
    "in" => DynamicDiscreteChoice.MarkovChain(endostates, [1.0 0.0; 1.0 0.0])*exochain
  ))
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

  res = DynamicDiscreteChoice.value(ddc, show_trace=true, method=:anderson, m=1)

  # Value function from past calculation
  Vexpect = reshape([1.2322159644235953, 0.6615604023460782, 1.2322159644235953, 1.0615604023460781, 1.2866470931936551, 1.3436909484346855, 1.2866470931936551, 1.4436909484346854],size(res.v))
  @test Vexpect ≈ res.v

  p = DynamicDiscreteChoice.choicep(res.v, ddc)
  pexpect = reshape([0.6389144289123461, 0.3610855710876540, 0.5425606483367234, 0.4574393516632768, 0.4857429020315101, 0.5142570979684898, 0.4608195280318177, 0.5391804719681823], size(p))
  @test pexpect ≈ p

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
  @testset "Estimate" begin
    T = 20_000
    Random.seed!(8760)
    sim = DynamicDiscreteChoice.simulate(T, ddc, res.v)
    est = DynamicDiscreteChoice.estimate_series(sim.actions, sim.states, ddc.discount, zero_action="out", states=ddc.states, actions=ddc.actions)
    payoffsexpect = [0.0000000000000002, -0.5052696682597060, 0.0000000000000002, -0.0911048412094677, 0.0000000000000000, 0.1256519392163353, 0.0000000000000000, 0.1596524777057666]
    @test vec(est.payoffs) ≈ payoffsexpect

    # different interfaces to estimate
    @test est.payoffs ≈ DynamicDiscreteChoice.estimate(est.choicep, est.transarray, est.ddc, zero_action="out").payoffs
    @test est.payoffs ≈ DynamicDiscreteChoice.estimate(est.choicep, est.ddc.transition, est.ddc, zero_action="out").payoffs



    # Test compatibility with ForwardDiff
    vec2named = (x,na)->NamedArray(reshape(x,size(na)), tuple(names(na)...), tuple(dimnames(na)...))
    @test est.choicep ≈ vec2named(vec(est.choicep),est.choicep)
    estgivenp = x->DynamicDiscreteChoice.estimate(vec2named(x,est.choicep),est.transarray,est.ddc, zero_action="out").payoffs
    @test isa(ForwardDiff.jacobian(estgivenp, vec(est.choicep)), AbstractMatrix)
    estgivent = x->DynamicDiscreteChoice.estimate(est.choicep,vec2named(x,est.transarray),est.ddc, zero_action="out").payoffs
    @test isa(ForwardDiff.jacobian(estgivent, vec(est.transarray)), AbstractMatrix)

    @testset "Bootstrap" begin
      B = 99
      bpayoffs = DynamicDiscreteChoice.bootstrap_series(est.payoffs,est.v,T,est.ddc,B=B)
      btable = DynamicDiscreteChoice.bootstrap_table(est.payoffs, bpayoffs)
      @test isa(btable, NamedArrays.NamedArray)
      for j in axes(payoffs)[2]
        @test payoffs[2,j] > btable[4,j][1]
        @test payoffs[2,j] < btable[4,j][2]
      end
    end
  end
end
