include("InformativeObservation.jl")
using .InformativeObservation, TOML, LinearAlgebra, Test
c=config();s=settings(); rows=NamedTuple[]
@testset "Piecewise integration at clipping boundaries" begin
 for (lambda,rep) in [(2,2),(0,33),(1,30),(2,13)]
  f=TOML.parsefile(joinpath(studyroot,"outputs","confirmatory_gh21","lambda$lambda","rep"*lpad(rep,4,'0')*".toml"))
  d=simulate(c,f["seed"],lambda,f["alpha"])
  fits=[fit(d,c;nodes=n) for n in (16,32)]
  x=fits[1]["theta"];q=gl(16)
  F,G,_=objective(x,d,c,q)
  ng=[begin delta=zeros(4);delta[j]=1e-5;(objective(x+delta,d,c,q;derivatives=false)[1]-objective(x-delta,d,c,q;derivatives=false)[1])/2e-5 end for j in 1:4]
  H=exact_hessian(x,d,c,q);H2=exact_hessian(x,d,c,q;step=5e-5)
  de=maximum(abs.(fits[1]["estimate"]-fits[2]["estimate"]));ds=maximum(abs.(fits[1]["se"]-fits[2]["se"]))
  @test all(z->z["valid"],fits)
  @test de<1e-4
  @test ds<1e-4
  @test maximum(abs.(G-ng))<1e-4
  @test maximum(abs.(H-H2))<.002
  push!(rows,(lambda=lambda,replication=rep,estimate_difference=de,se_difference=ds,gradient_error=maximum(abs.(G-ng)),hessian_step_error=maximum(abs.(H-H2))))
  println(rows[end]);flush(stdout)
 end
 d=simulate(c,1234,0.,log(.81/.19);nseq=20)
 empty=Dataset(d.z0,d.plan,d.action,falses(size(d.y)),fill(Int8(-1),size(d.y)),0.,log(.81/.19))
 x=[.2,-.5,.3,1.]
 @test abs(objective(x,empty,c,gl(16))[1])<1e-12
 @test objective(x,d,c,gl(16);aware=true)[1]-objective(x,d,c,gl(16))[1]≈-sum(d.observed)*log(.81)-sum(.!d.observed)*log(.19) atol=1e-10
end
write_csv(joinpath(studyroot,"validation","piecewise_diagnostic.csv"),rows)
