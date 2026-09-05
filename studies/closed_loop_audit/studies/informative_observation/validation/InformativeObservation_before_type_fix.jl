module InformativeObservation
using LinearAlgebra, Random, Statistics, TOML, SHA, Printf, Dates
export Dataset, config, settings, gh, simulate, objective, fit, write_dataset, read_dataset,
       write_toml, write_csv, logistic, studyroot, filehash
const studyroot = normpath(joinpath(@__DIR__, ".."))
config() = TOML.parsefile(joinpath(studyroot,"config","design.toml"))
settings() = TOML.parsefile(joinpath(studyroot,"config","execution.toml"))
filehash(path) = bytes2hex(sha256(read(path)))
logistic(x) = x>=0 ? inv(1+exp(-x)) : exp(x)/(1+exp(x))
softplus(x) = max(x,0)+log1p(exp(-abs(x)))
logbern(y,x) = y==1 ? -softplus(-x) : -softplus(x)
struct Dataset
    z0::Vector{Float64}
    plan::Matrix{Float64}
    action::Matrix{Float64}
    observed::BitMatrix
    y::Matrix{Int8}
    lambda::Float64
    alpha::Float64
end
function gh(n)
    e=eigen(SymTridiagonal(zeros(n),sqrt.(collect(1.0:n-1))))
    (nodes=e.values,weights=vec(e.vectors[1,:].^2))
end
function simulate(c,seed,lambda,alpha;nseq=c["population"]["learners"]*c["population"]["items_per_learner"])
    pop=c["population"]; p=c["platform"]; l=c["learner"]; r=c["response"]
    tmax=Int(pop["opportunities"]); rng=Xoshiro(seed)
    pi,pa,pv,pn,en,br,bi0,bi1=Float64.((p["plan_intercept"],p["adaptivity"],p["previous_failure_coefficient"],p["plan_noise_sd"],p["enactment_noise_sd"],p["belief_retention"],p["belief_success_increment"],p["belief_failure_increment"]))
    ri,rd,rs,rv=Float64.((r["intercept"],r["delay_coefficient"],r["state_coefficient"],r["previous_failure_coefficient"]))
    phi,mi0,mi1=Float64.((l["retention"],l["success_increment"],l["failure_increment"]))
    mlo,mhi=Float64.(l["state_clip"]); blo,bhi=Float64.(p["belief_clip"])
    # All random arrays are drawn before any observation-dependent branch.
    initial=randn(rng,nseq).*sqrt(pop["initial_state_variance"])
    z0=initial+randn(rng,nseq).*sqrt(pop["baseline_error_variance"])
    ep=randn(rng,nseq,tmax); ea=randn(rng,nseq,tmax)
    ur=rand(rng,nseq,tmax); uy=rand(rng,nseq,tmax)
    plan=zeros(nseq,tmax); action=similar(plan); obs=falses(nseq,tmax); y=fill(Int8(-1),nseq,tmax)
    M=copy(initial); B=copy(z0); V=zeros(nseq)
    @inbounds for t in 1:tmax, i in 1:nseq
        plan[i,t]=pi+pa*B[i]+pv*V[i]+pn*ep[i,t]
        action[i,t]=plan[i,t]+en*ea[i,t]
        obs[i,t]=ur[i,t]<logistic(alpha+lambda*M[i])
        if obs[i,t]
            eta=ri+rd*action[i,t]+rs*M[i]+rv*V[i]
            y[i,t]=uy[i,t]<logistic(eta)
            M[i]=clamp(phi*M[i]+(y[i,t]==0 ? mi0 : mi1),mlo,mhi)
            B[i]=clamp(br*B[i]+(y[i,t]==0 ? bi0 : bi1),blo,bhi)
            V[i]=y[i,t]
        else
            M[i]=clamp(phi*M[i],mlo,mhi)
            B[i]=clamp(br*B[i],blo,bhi)
        end
    end
    Dataset(z0,plan,action,obs,y,Float64(lambda),Float64(alpha))
end

"""Conditional marginal NLL with analytic first/second derivatives.
Integrates M1|Z0, advances M on every opportunity, and keeps all-zero sequences.
The optional selection term is inside the sequence integral; alpha/lambda are fixed.
Derivatives of a clipped state are zero strictly outside its unclipped bounds.
"""
function objective(theta,d::Dataset,c,q;aware=false,derivatives=true)
    bd,bs,bv,gamma=theta
    phi=logistic(gamma); dp=phi*(1-phi); ddp=dp*(1-2phi)
    l=c["learner"]; intercept=Float64(c["response"]["intercept"])
    lo,hi=Float64.(l["state_clip"]); success=Float64(l["success_increment"]); failure=Float64(l["failure_increment"])
    pm=Float64(c["estimation"]["posterior_mean_multiplier"]); sd=sqrt(Float64(c["estimation"]["posterior_variance"]))
    n,tmax=size(d.y); nq=length(q.nodes)
    ll=zeros(nq); scores=zeros(4,nq); hs=zeros(4,4,nq)
    G=zeros(4); H=zeros(4,4); value=0.0
    @inbounds for i in 1:n
        fill!(scores,0); fill!(hs,0)
        for k in 1:nq
            m=pm*d.z0[i]+sd*q.nodes[k]; dm=0.0; ddm=0.0; v=0.0
            lk=log(q.weights[k])
            for t in 1:tmax
                observed=d.observed[i,t]
                if aware
                    s=d.alpha+d.lambda*m; ps=logistic(s)
                    lk+=logbern(observed,s)
                    if derivatives
                        a=d.lambda*dm
                        scores[4,k]+=(observed-ps)*a
                        hs[4,4,k]+=(observed-ps)*d.lambda*ddm-ps*(1-ps)*a*a
                    end
                end
                if observed
                    yt=d.y[i,t]; a1=d.action[i,t]; eta=intercept+bd*a1+bs*m+bv*v
                    py=logistic(eta); lk+=logbern(yt,eta)
                    if derivatives
                        a=(a1,m,v,bs*dm); residual=yt-py; w=py*(1-py)
                        for u in 1:4
                            scores[u,k]+=residual*a[u]
                            for z in 1:4
                                hs[u,z,k]-=w*a[u]*a[z]
                            end
                        end
                        hs[2,4,k]+=residual*dm; hs[4,2,k]+=residual*dm
                        hs[4,4,k]+=residual*bs*ddm
                    end
                    inc=yt==0 ? success : failure; v=yt
                else
                    inc=0.0
                end
                raw=phi*m+inc
                newdm=dp*m+phi*dm; newddm=ddp*m+2dp*dm+phi*ddm
                if lo<raw<hi
                    m=raw; dm=newdm; ddm=newddm
                else
                    m=clamp(raw,lo,hi); dm=0.0; ddm=0.0
                end
            end
            ll[k]=lk
        end
        mx=maximum(ll); denom=sum(exp(x-mx) for x in ll)
        value-=mx+log(denom)
        if derivatives
            sg=zeros(4)
            for k in 1:nq
                w=exp(ll[k]-mx)/denom
                for u in 1:4
                    sg[u]+=w*scores[u,k]
                    for z in 1:4
                        H[u,z]-=w*(hs[u,z,k]+scores[u,k]*scores[z,k])
                    end
                end
            end
            G-=sg
            H+=sg*sg'
        end
    end
    value,G,H
end
function optimize_start(x0,d,c,q,s;aware=false)
    x=Float64.(x0); iterations=0; reason="maximum_iterations"
    for iter in 1:s["max_iterations"]
        iterations=iter
        f,g,H=objective(x,d,c,q;aware)
        if !(isfinite(f)&&all(isfinite,g)&&all(isfinite,H)); reason="nonfinite"; break end
        if maximum(abs.(g))<s["gradient_tolerance"]; reason="gradient"; break end
        ev=eigvals(Symmetric(H)); shift=max(0.0,1e-4-minimum(ev))
        step=-(Symmetric(H)+shift*I)\g
        if norm(step)>3; step*=3/norm(step) end
        slope=dot(g,step); scale=1.0; accepted=false
        for bt in 1:35
            candidate=x+scale*step
            nf=objective(candidate,d,c,q;aware,derivatives=false)[1]
            if isfinite(nf)&&nf<=f+1e-4*scale*slope
                x=candidate; accepted=true; break
            end
            scale*=0.5
        end
        if !accepted; reason="line_search"; break end
    end
    f,g,H=objective(x,d,c,q;aware)
    (theta=x,nll=f,gradient=g,hessian=H,iterations=iterations,reason=reason)
end
function fit(d,c,s=settings();nodes=c["estimation"]["gh_nodes"],aware=false)
    q=gh(nodes); started=time_ns()
    fits=[optimize_start(x,d,c,q,s;aware) for x in s["starts"]]
    valid_objectives=[isfinite(x.nll) ? x.nll : Inf for x in fits]
    best=fits[argmin(valid_objectives)]
    x=best.theta; H=best.hessian; phi=logistic(x[4]); estimates=[x[1:3];phi]
    ev=eigvals(Symmetric(H)); condition=cond(H); grad=maximum(abs.(best.gradient))
    boundary=abs(x[4])>=s["gamma_boundary"]||maximum(abs.(x[1:3]))>=s["parameter_boundary"]
    valid=all(isfinite,x)&&isfinite(best.nll)&&grad<s["gradient_tolerance"]&&minimum(ev)>0&&condition<s["maximum_hessian_condition"]&&!boundary
    status=valid ? "valid" : boundary ? "boundary" : minimum(ev)<=0 ? "hessian_not_positive" : condition>=s["maximum_hessian_condition"] ? "ill_conditioned" : "not_converged"
    se=fill(NaN,4); lower=fill(NaN,4); upper=fill(NaN,4); se_theta=fill(NaN,4)
    if valid
        covariance=inv(Symmetric(H))
        se_theta=sqrt.(diag(covariance)); se=copy(se_theta); se[4]*=phi*(1-phi)
        lower=estimates-1.959963984540054se; upper=estimates+1.959963984540054se
        lower[4]=logistic(x[4]-1.959963984540054se_theta[4]); upper[4]=logistic(x[4]+1.959963984540054se_theta[4])
    end
    Dict("theta"=>x,"estimate"=>estimates,"se"=>se,"se_theta"=>se_theta,"lower"=>lower,"upper"=>upper,
        "nll"=>best.nll,"gradient_max"=>grad,"hessian_min_eigenvalue"=>minimum(ev),"hessian_condition"=>condition,
        "valid"=>valid,"status"=>status,"boundary"=>boundary,"nodes"=>nodes,"aware"=>aware,
        "observed_rate"=>mean(d.observed),"all_missing_sequences"=>count(i->!any(d.observed[i,:]),axes(d.y,1)),
        "iterations"=>best.iterations,"stop_reason"=>best.reason,"start_nll"=>[f.nll for f in fits],
        "start_gradient_max"=>[maximum(abs.(f.gradient)) for f in fits],"elapsed_seconds"=>(time_ns()-started)/1e9)
end
function write_toml(path,obj)
    mkpath(dirname(path)); tmp=path*".tmp-"*string(getpid())
    open(tmp,"w") do io; TOML.print(io,obj;sorted=true) end
    mv(tmp,path;force=true)
end
function write_csv(path,rows)
    isempty(rows)&&error("No rows to write")
    mkpath(dirname(path)); names=keys(first(rows))
    open(path,"w") do io
        println(io,join(string.(names),","))
        for row in rows
            println(io,join((replace(string(getproperty(row,k)),","=>";") for k in names),","))
        end
    end
end
function write_dataset(path,d)
    rows=[(sequence=i,opportunity=t,z0=d.z0[i],plan=d.plan[i,t],action=d.action[i,t],observed=Int(d.observed[i,t]),y=d.y[i,t],lambda=d.lambda,alpha=d.alpha) for t in axes(d.y,2) for i in axes(d.y,1)]
    write_csv(path,rows)
end
function read_dataset(path)
    lines=readlines(path); rows=[parse.(Float64,split(x,',')) for x in lines[2:end]]
    n=Int(maximum(x[1] for x in rows)); t=Int(maximum(x[2] for x in rows))
    z0=zeros(n); p=zeros(n,t); a=zeros(n,t); r=falses(n,t); y=fill(Int8(-1),n,t)
    for x in rows
        i,j=Int.(x[1:2]); z0[i]=x[3]; p[i,j]=x[4]; a[i,j]=x[5]; r[i,j]=x[6]==1; y[i,j]=Int8(x[7])
    end
    Dataset(z0,p,a,r,y,rows[1][8],rows[1][9])
end
end
