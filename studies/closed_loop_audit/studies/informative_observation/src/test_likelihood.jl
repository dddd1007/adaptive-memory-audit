include("InformativeObservation.jl")
using .InformativeObservation, LinearAlgebra, Test, TOML
c=config(); s=settings(); q=gh(21); x=[.2,-.5,.3,log(.7/.3)]
dest=joinpath(studyroot,"validation"); mkpath(dest)
@testset "Conditional likelihood and DGP" begin
    @test sum(q.weights)≈1 atol=1e-14
    @test sum(q.weights.*q.nodes)≈0 atol=1e-13
    @test sum(q.weights.*q.nodes.^2)≈1 atol=1e-13
    @test sum(q.weights.*q.nodes.^4)≈3 atol=1e-12
    prior=c["population"]["initial_state_variance"]; noise=c["population"]["baseline_error_variance"]
    @test prior/(prior+noise)≈c["estimation"]["posterior_mean_multiplier"]
    @test prior*noise/(prior+noise)≈c["estimation"]["posterior_variance"]
    alpha=log(.81/.19)
    d=simulate(c,c["master_seed"]+c["validation_seed_offset"]+9999,0.0,alpha;nseq=20)
    empty=Dataset(d.z0,d.plan,d.action,falses(size(d.y)),fill(Int8(-1),size(d.y)),0.0,alpha)
    @test abs(objective(x,empty,c,q)[1])<1e-12
    @test maximum(abs.(objective(x,empty,c,q)[2]))<1e-12
    @test all(d.y[.!d.observed].==-1)
    @test all(z->z in (0,1),d.y[d.observed])
    same=simulate(c,c["master_seed"]+c["validation_seed_offset"]+9999,0.0,alpha;nseq=20)
    @test d.y==same.y && d.action==same.action && d.z0==same.z0
    constant=-sum(d.observed)*log(.81)-sum(.!d.observed)*log(.19)
    for theta in (x,[-.1,-.9,.5,1.6],[.6,-.2,-.1,-.3])
        f,g,H=objective(theta,d,c,q)
        af,ag,aH=objective(theta,d,c,q;aware=true)
        @test af-f≈constant atol=1e-10
        @test ag≈g atol=1e-10
        @test aH≈H atol=1e-9
    end
    checks=NamedTuple[]
    for lambda in (0.0,1.0,2.0), aware in (false,true)
        data=simulate(c,c["master_seed"]+c["validation_seed_offset"]+9998,lambda,alpha;nseq=30)
        f,g,H=objective(x,data,c,q;aware)
        eps=1e-5; eye=Matrix{Float64}(I,4,4)
        numeric_g=[(objective(x+eps*eye[:,j],data,c,q;aware,derivatives=false)[1]-objective(x-eps*eye[:,j],data,c,q;aware,derivatives=false)[1])/(2eps) for j in 1:4]
        numeric_H=hcat([(objective(x+eps*eye[:,j],data,c,q;aware)[2]-objective(x-eps*eye[:,j],data,c,q;aware)[2])/(2eps) for j in 1:4]...)
        gd=maximum(abs.(g-numeric_g)); hd=maximum(abs.(H-numeric_H))
        @test gd<s["gradient_check_atol"]
        @test hd<s["hessian_check_atol"]
        push!(checks,(lambda=lambda,aware=aware,gradient_error=gd,hessian_error=hd))
    end
    write_csv(joinpath(dest,"derivative_checks.csv"),checks)
    # Check deterministic skip recursion against a one-node, one-sequence calculation.
    d1=Dataset([1.0],zeros(1,2),[0.0 0.0],BitMatrix([false true]),Int8[-1 1],0.0,alpha)
    one=(nodes=[0.0],weights=[1.0]); phi=logistic(x[4])
    expected=-log(logistic(c["response"]["intercept"]+x[2]*phi*.9))
    @test objective(x,d1,c,one)[1]≈expected atol=1e-12
    # Clipped paths are also differentiated consistently away from the boundary.
    dc=Dataset([30.0],zeros(1,3),zeros(1,3),trues(1,3),zeros(Int8,1,3),1.0,alpha)
    f,g,H=objective(x,dc,c,q;aware=true)
    @test isfinite(f)&&all(isfinite,g)&&all(isfinite,H)
end
