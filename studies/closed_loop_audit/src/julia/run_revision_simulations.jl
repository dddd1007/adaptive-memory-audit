#!/usr/bin/env julia

"""
Authoritative known-truth simulations for the adaptive-selection manuscript.

The data-generating processes, estimators, Monte Carlo loops, and dynamic-
parameter recovery study are implemented in Julia.  R is used only as an
independent refitting and CR2/Satterthwaite validation layer on frozen data.
No Python result is consumed by this script or by the manuscript tables.
"""

include(joinpath(@__DIR__, "BRMCore.jl"))
using .BRMCore
using LinearAlgebra
using Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = ensure_directory(joinpath(ROOT, "outputs_julia"))
const VALID = ensure_directory(joinpath(ROOT, "validation", "frozen_data"))
const MASTER_SEED = 20260903

clip1(x, lo, hi) = min(max(Float64(x), Float64(lo)), Float64(hi))
normal01(rng) = rand_normal!(rng)
bernoulli(rng, p) = rand_bernoulli!(rng, clip1(p, 0.0, 1.0))

function covariance1(x, y)
    n = min(length(x), length(y))
    mx, my = mean1(x), mean1(y)
    sum((Float64(x[i]) - mx) * (Float64(y[i]) - my) for i in 1:n) / (n - 1)
end

function append_column!(target::Vector{Float64}, values)
    for value in values
        push!(target, Float64(value))
    end
end

function append_column!(target::Vector{Int}, values)
    for value in values
        push!(target, Int(value))
    end
end

function design_matrix(columns...)
    n = length(columns[1])
    X = zeros(Float64, n, length(columns))
    for j in eachindex(columns), i in 1:n
        X[i, j] = Float64(columns[j][i])
    end
    X
end

function logistic_full(X::Matrix{Float64}, y::Vector{Float64}; max_iter=80)
    n, p0 = size(X)
    X1 = hcat(ones(Float64, n), X)
    p = p0 + 1
    beta = zeros(Float64, p)
    converged = false
    for _ in 1:max_iter
        eta = X1 * beta
        prob = sigmoid.(clamp.(eta, -35.0, 35.0))
        w = max.(prob .* (1 .- prob), 1.0e-7)
        info = X1' * (w .* X1) + 1.0e-10I
        score = X1' * (y - prob)
        step = info \ score
        beta += step
        if maximum(abs.(step)) < 1.0e-9
            converged = true
            break
        end
    end
    prob = sigmoid.(clamp.(X1 * beta, -35.0, 35.0))
    w = max.(prob .* (1 .- prob), 1.0e-7)
    bread = inv(Symmetric(X1' * (w .* X1) + 1.0e-10I))
    beta, prob, Matrix(bread), converged, X1, w
end

"""CR3 cluster linearization with a t(G-1) reference."""
function logistic_cr3(X::Matrix{Float64}, y::Vector{Float64}, clusters::Vector{Int})
    beta, prob, bread, converged, X1, w = logistic_full(X, y)
    Z = sqrt.(w) .* X1
    pearson = (y - prob) ./ sqrt.(w)
    meat = zeros(Float64, size(X1, 2), size(X1, 2))
    for g in sort(unique(clusters))
        idx = findall(==(g), clusters)
        Zg = Z[idx, :]
        rg = pearson[idx]
        small = Matrix(I, size(X1, 2), size(X1, 2)) - bread * (Zg' * Zg)
        adjusted = rg + Zg * (pinv(small; rtol=1.0e-10) * bread * (Zg' * rg))
        ug = Zg' * adjusted
        meat += ug * ug'
    end
    vcov = bread * meat * bread
    se = sqrt.(max.(diag(vcov), 0.0))
    beta, se, converged
end

"""CR2/Satterthwaite inference for a logistic M-estimator.

The adjustment is applied in Pearson-residual coordinates with the inverse
Fisher information as the working-model bread.  Frozen-data refits are checked
against clubSandwich::vcovCR(type = "CR2").
"""
function logistic_cr2_satt(X::Matrix{Float64}, y::Vector{Float64}, clusters::Vector{Int}; focal=2)
    beta, prob, bread, converged, X1, w = logistic_full(X, y)
    Z = sqrt.(w) .* X1
    pearson = (y - prob) ./ sqrt.(w)
    meat = zeros(Float64, size(X1,2), size(X1,2))
    c = zeros(Float64,size(X1,2)); c[focal]=1.0
    qnorm2 = Float64[]; hrows = Vector{Vector{Float64}}()
    evals,evecs=eigen(Symmetric(bread))
    Bhalf=evecs*Diagonal(sqrt.(max.(evals,0.0)))*evecs'
    for g in sort(unique(clusters))
        idx=findall(==(g),clusters); Zg=Z[idx,:]
        U=Zg*Bhalf
        vals,vecs=eigen(Symmetric(U'*U))
        keep=findall(x->x>1.0e-12,vals)
        function apply_A(v)
            isempty(keep) && return v
            lam=vals[keep]; V=vecs[:,keep]
            Q=U*(V*Diagonal(1 ./ sqrt.(lam)))
            factors=1 ./ sqrt.(max.(1 .- lam,1.0e-10)) .- 1
            v+Q*(factors.*(Q'*v))
        end
        adj=apply_A(pearson[idx]); ug=Zg'*adj; meat+=ug*ug'
        q=apply_A(Zg*bread*c); push!(qnorm2,dot(q,q))
        push!(hrows,vec(q' * Zg * Bhalf))
    end
    vcov=bread*meat*bread
    se=sqrt.(max.(diag(vcov),0.0))
    H=reduce(vcat,(reshape(h,1,:) for h in hrows))
    P=Diagonal(qnorm2)-H*H'
    df=tr(P)^2/sum(P.^2)
    beta,se,df,converged
end

"""Accurate 97.5% t critical approximation for df above four."""
function tcrit975(df::Real)
    v = max(Float64(df), 4.01)
    z = 1.959963984540054
    z + (z^3 + z) / (4v) + (5z^5 + 16z^3 + 3z) / (96v^2) +
        (3z^7 + 19z^5 + 17z^3 - 15z) / (384v^3)
end

"""CR2 and Satterthwaite inference for the centered assignment coefficient."""
function ols_cr2_satt(X::Matrix{Float64}, y::Vector{Float64}, clusters::Vector{Int}; focal=2)
    Binv = inv(Symmetric(X' * X))
    beta = Binv * (X' * y)
    resid = y - X * beta
    c = zeros(Float64, size(X, 2)); c[focal] = 1.0
    meat = zeros(Float64, size(X, 2), size(X, 2))
    qnorm2 = Float64[]; hrows = Vector{Vector{Float64}}()
    evals, evecs = eigen(Symmetric(Binv))
    Bhalf = evecs * Diagonal(sqrt.(max.(evals, 0.0))) * evecs'
    for g in sort(unique(clusters))
        idx = findall(==(g), clusters)
        Xg = X[idx, :]
        U = Xg * Bhalf
        vals, vecs = eigen(Symmetric(U' * U))
        keep = findall(x -> x > 1.0e-12, vals)
        function apply_A(v)
            isempty(keep) && return v
            lam = vals[keep]
            V = vecs[:, keep]
            Q = U * (V * Diagonal(1 ./ sqrt.(lam)))
            factors = 1 ./ sqrt.(max.(1 .- lam, 1.0e-10)) .- 1
            v + Q * (factors .* (Q' * v))
        end
        adj = apply_A(resid[idx])
        ug = Xg' * adj
        meat += ug * ug'
        q = apply_A(Xg * Binv * c)
        push!(qnorm2, dot(q, q))
        push!(hrows,vec(q' * Xg * Bhalf))
    end
    vcov = Binv * meat * Binv
    se = sqrt(max(vcov[focal, focal], 0.0))
    H=reduce(vcat,(reshape(h,1,:) for h in hrows))
    P=Diagonal(qnorm2)-H*H'
    df=tr(P)^2/sum(P.^2)
    beta[focal], se, df
end

function simulate_spaced(seed::Int, n_people::Int, adaptivity::Float64;
                         reliability=0.90, context_strength=0.0,
                         randomized_delta=0.0, compliance=1.0,
                         oracle_policy=false, context_proxy_r2=-1.0,
                         force_randomize=false)
    rng = DeterministicRNG(seed)
    n_items = 100
    n_pairs = n_people * n_items
    ability = [0.65 * normal01(rng) for _ in 1:n_people]
    difficulty = [0.65 * normal01(rng) for _ in 1:n_items]
    memory = zeros(Float64, n_pairs)
    person = zeros(Int, n_pairs)
    for pair in 1:n_pairs
        i = div(pair - 1, n_items) + 1
        j = mod(pair - 1, n_items) + 1
        person[pair] = i
        memory[pair] = ability[i] - difficulty[j] + 0.35 * normal01(rng)
    end
    vm = var(memory)
    noise_sd = sqrt(vm * (1 - reliability) / reliability)
    belief = [memory[i] + noise_sd * normal01(rng) for i in 1:n_pairs]
    previous = zeros(Float64, n_pairs)
    data = Dict{Symbol,Any}(
        :person=>Int[], :failure=>Float64[], :previous=>Float64[],
        :memory=>Float64[], :belief=>Float64[], :plan=>Float64[],
        :actual=>Float64[], :context=>Float64[], :proxy=>Float64[],
        :assign=>Float64[], :true_itt=>Float64[])
    for _ in 1:2
        policy_state = oracle_policy ? memory : belief
        plan = [0.75 + adaptivity * policy_state[i] - 0.55 * previous[i] +
                0.20 * normal01(rng) for i in 1:n_pairs]
        context = [normal01(rng) for _ in 1:n_pairs]
        proxy = context_proxy_r2 < 0 ? zeros(n_pairs) :
                [sqrt(context_proxy_r2) * context[i] +
                 sqrt(1 - context_proxy_r2) * normal01(rng) for i in 1:n_pairs]
        natural = [0.50 * normal01(rng) + 0.40 * context_strength * context[i]
                   for i in 1:n_pairs]
        randomized = randomized_delta > 0 || force_randomize
        assign = randomized ? [bernoulli(rng, 0.5) for _ in 1:n_pairs] : zeros(n_pairs)
        comply = randomized ? [bernoulli(rng, compliance) for _ in 1:n_pairs] : zeros(n_pairs)
        actual = [plan[i] + natural[i] + randomized_delta * (2assign[i] - 1) * comply[i]
                  for i in 1:n_pairs]
        eta = [-2.40 + 0.25actual[i] - 0.55memory[i] + 0.35previous[i] +
               0.50context_strength * context[i] for i in 1:n_pairs]
        failure = [bernoulli(rng, sigmoid(eta[i])) for i in 1:n_pairs]
        true_itt = fill(NaN, n_pairs)
        if randomized
            for i in 1:n_pairs
                p_hi = sigmoid(-2.40 + 0.25 * (plan[i] + natural[i] + randomized_delta) -
                               0.55memory[i] + 0.35previous[i] + 0.50context_strength*context[i])
                p_lo = sigmoid(-2.40 + 0.25 * (plan[i] + natural[i] - randomized_delta) -
                               0.55memory[i] + 0.35previous[i] + 0.50context_strength*context[i])
                true_itt[i] = compliance * (p_hi - p_lo)
            end
        end
        for (key, values) in ((:person,person),(:failure,failure),(:previous,previous),
                              (:memory,memory),(:belief,belief),(:plan,plan),(:actual,actual),
                              (:context,context),(:proxy,proxy),(:assign,assign),(:true_itt,true_itt))
            append_column!(data[key], values)
        end
        success = 1 .- failure
        memory = [0.75memory[i] + 0.45success[i] - 0.65failure[i] + 0.20normal01(rng)
                  for i in 1:n_pairs]
        belief = 0.75 .* belief + 0.45 .* success - 0.65 .* failure
        previous = failure
    end
    data
end

function simulate_tutor(seed::Int, adaptivity::Float64; n_people=60, n_skills=6,
                        n_opps=8, randomized_delta=0.0, compliance=1.0,
                        force_randomize=false, family="RL", calibrated_initial=false)
    rng = DeterministicRNG(seed)
    n_pairs = n_people * n_skills
    ability = [0.80normal01(rng) for _ in 1:n_people]
    challenge = [0.45normal01(rng) for _ in 1:n_skills]
    theta = zeros(Float64, n_pairs); person = zeros(Int, n_pairs)
    for pair in 1:n_pairs
        i = div(pair - 1, n_skills) + 1; j = mod(pair - 1, n_skills) + 1
        person[pair] = i
        theta[pair] = 0.15 + ability[i] - challenge[j] + 0.25normal01(rng)
    end
    baseline = calibrated_initial ? copy(theta) :
               [theta[i] + 0.20normal01(rng) for i in 1:n_pairs]
    belief = [theta[i] + 0.20normal01(rng) for i in 1:n_pairs]
    recon = copy(baseline)
    variance = fill(0.80, n_pairs)
    previous = fill(0.50, n_pairs)
    data = Dict{Symbol,Any}(
        :person=>Int[], :success=>Float64[], :previous=>Float64[],
        :theta=>Float64[], :belief=>Float64[], :recon=>Float64[],
        :plan=>Float64[], :actual=>Float64[], :assign=>Float64[],
        :true_itt=>Float64[], :opportunity=>Float64[])
    ymat = zeros(Float64, n_pairs, n_opps)
    dmat = zeros(Float64, n_pairs, n_opps)
    for opp in 1:n_opps
        plan = [clip1(adaptivity*belief[i] - 0.35 + 0.10normal01(rng), -2.25, 2.25)
                for i in 1:n_pairs]
        base = [round(clip1(plan[i] + 0.25normal01(rng), -2.25, 2.25)/0.25)*0.25
                for i in 1:n_pairs]
        randomized = randomized_delta > 0 || force_randomize
        assign = randomized ? [bernoulli(rng, 0.5) for _ in 1:n_pairs] : zeros(n_pairs)
        comply = randomized ? [bernoulli(rng, compliance) for _ in 1:n_pairs] : zeros(n_pairs)
        actual = [base[i] + randomized_delta*(2assign[i]-1)*comply[i] for i in 1:n_pairs]
        prob = [sigmoid(0.35 + theta[i] - actual[i]) for i in 1:n_pairs]
        success = [bernoulli(rng, prob[i]) for i in 1:n_pairs]
        true_itt = fill(NaN, n_pairs)
        if randomized
            for i in 1:n_pairs
                p_h = sigmoid(0.35 + theta[i] - (base[i] + randomized_delta))
                p_e = sigmoid(0.35 + theta[i] - (base[i] - randomized_delta))
                true_itt[i] = compliance * ((1-p_h) - (1-p_e))
            end
        end
        for (key, values) in ((:person,person),(:success,success),(:previous,previous),
                              (:theta,theta),(:belief,belief),(:recon,recon),(:plan,plan),
                              (:actual,actual),(:assign,assign),(:true_itt,true_itt),
                              (:opportunity,fill(Float64(opp),n_pairs)))
            append_column!(data[key], values)
        end
        ymat[:,opp] = success; dmat[:,opp] = actual
        if family == "RL"
            theta = clamp.(theta + 0.30 .* (success - prob) .+ 0.04, -4.5, 4.5)
        else
            predvar = variance .+ 0.15
            post = 1 ./ (1 ./ predvar + max.(prob .* (1 .- prob), 1.0e-5))
            theta = clamp.(theta + post .* (success - prob) .+ 0.04, -4.5, 4.5)
            variance = post
        end
        pq = sigmoid.(0.35 .+ belief - actual)
        belief = clamp.(belief + 0.25 .* (success-pq) .+ 0.04, -4.5, 4.5)
        pr = sigmoid.(0.35 .+ recon - actual)
        recon = clamp.(recon + 0.30 .* (success-pr) .+ 0.04, -4.5, 4.5)
        previous = success
    end
    data, (baseline=baseline, ymat=ymat, dmat=dmat)
end

"""Exact dynamic benchmark for the spaced-retrieval architecture.

The recorded baseline `z0` is the learner state at the first decision point and
is available to both the platform and analyst.  The platform belief is updated
with fixed, distinct parameters and drives selection.  The learner state is
deterministic conditional on `z0` and observed responses, which makes the
conditional likelihood correctly specified in the clean benchmark.

When `context_strength > 0`, an unrecorded context changes enacted delay and
retrieval.  When `observation_strength > 0`, the same context also changes the
probability that the opportunity is observed.  These two switches distinguish
outcome-model omission with complete observation from a nonignorable
observation mechanism.
"""
function simulate_spaced_dynamic(seed::Int, adaptivity::Float64;
                                 n_people=30, n_items=20, n_opps=8,
                                 context_strength=0.0,
                                 observation_strength=0.0)
    rng = DeterministicRNG(seed)
    n_pairs = n_people * n_items
    ability = [0.65 * normal01(rng) for _ in 1:n_people]
    difficulty = [0.65 * normal01(rng) for _ in 1:n_items]
    memory = zeros(Float64, n_pairs)
    for pair in 1:n_pairs
        i = div(pair - 1, n_items) + 1
        j = mod(pair - 1, n_items) + 1
        memory[pair] = ability[i] - difficulty[j] + 0.35 * normal01(rng)
    end
    z0 = copy(memory)
    belief = copy(z0)
    previous = zeros(Float64, n_pairs)
    ymat = fill(NaN, n_pairs, n_opps)
    amat = fill(NaN, n_pairs, n_opps)
    observed = falses(n_pairs, n_opps)
    for opp in 1:n_opps
        plan = [0.75 + adaptivity * belief[i] - 0.55 * previous[i] +
                0.20 * normal01(rng) for i in 1:n_pairs]
        context = [normal01(rng) for _ in 1:n_pairs]
        actual = [plan[i] + 0.50 * normal01(rng) +
                  0.40 * context_strength * context[i] for i in 1:n_pairs]
        response_prob = [sigmoid(-2.40 + 0.25 * actual[i] - 0.55 * memory[i] +
                                 0.35 * previous[i] +
                                 0.50 * context_strength * context[i])
                         for i in 1:n_pairs]
        for i in 1:n_pairs
            keep = observation_strength <= 0 ? true :
                   bernoulli(rng, sigmoid(1.40 + observation_strength *
                                          (0.80 * context[i] + 0.40 * memory[i]))) == 1
            if keep
                failure = Float64(bernoulli(rng, response_prob[i]))
                observed[i, opp] = true
                ymat[i, opp] = failure
                amat[i, opp] = actual[i]
                success = 1.0 - failure
                memory[i] = clip1(0.75 * memory[i] + 0.45 * success -
                                  0.65 * failure, -4.5, 4.5)
                belief[i] = clip1(0.70 * belief[i] + 0.40 * success -
                                  0.60 * failure, -4.5, 4.5)
                previous[i] = failure
            end
        end
    end
    (z0=z0, ymat=ymat, amat=amat, observed=observed,
     observed_rate=count(observed) / length(observed))
end

function observational_fits(data, architecture::String)
    if architecture == "spaced"
        action_z, scale = standardize(data[:actual]); plan_z,_=standardize(data[:plan])
        belief_z,_=standardize(data[:belief]); state_z,_=standardize(data[:memory])
        outcome=data[:failure]; previous=data[:previous]; truth=0.25scale
        specs=(("marginal",design_matrix(action_z)),
               ("recorded_belief",design_matrix(action_z,belief_z,previous)),
               ("policy_log",design_matrix(action_z,plan_z,previous)),
               ("oracle_state",design_matrix(action_z,state_z,previous)))
        state=data[:memory]
    else
        action_z, scale = standardize(data[:actual]); plan_z,_=standardize(data[:plan])
        belief_z,_=standardize(data[:belief]); recon_z,_=standardize(data[:recon])
        state_z,_=standardize(data[:theta]); outcome=1 .- data[:success]
        previous=data[:previous]; truth=scale
        specs=(("marginal",design_matrix(action_z)),
               ("recorded_belief",design_matrix(action_z,belief_z,previous)),
               ("reconstructed_history",design_matrix(action_z,recon_z,previous)),
               ("policy_log",design_matrix(action_z,plan_z,recon_z,previous)),
               ("oracle_state",design_matrix(action_z,state_z,previous)))
        state=data[:theta]
    end
    rows=NamedTuple[]
    G=length(unique(data[:person])); crit=tcrit975(G-1)
    for (label,X) in specs
        beta,se,df,ok=logistic_cr2_satt(X,outcome,data[:person])
        est=beta[2]; s=se[2]; crit=tcrit975(df)
        push!(rows,(estimator=label,truth=truth,estimate=est,se=s,
                    covered=(est-crit*s<=truth<=est+crit*s),sign_error=est<0,
                    df=df,converged=ok))
    end
    cov_as=covariance1(data[:actual],state); var_a=var(data[:actual])
    if architecture == "spaced"
        cov_ah=covariance1(data[:actual],data[:previous])
        prediction=scale*(0.25 - 0.55cov_as/var_a + 0.35cov_ah/var_a)
    else
        prediction=scale*(1-cov_as/var_a)
    end
    rows,(truth=truth,prediction=prediction,outcome_rate=mean1(outcome),
          plan_actual_corr=correlation(data[:plan],data[:actual]),
          action_state_corr=correlation(data[:actual],state))
end

function context_fit(data; include_proxy=false)
    a,scale=standardize(data[:actual]); p,_=standardize(data[:plan])
    X=include_proxy ? begin q,_=standardize(data[:proxy]); design_matrix(a,p,data[:previous],q) end :
                      design_matrix(a,p,data[:previous])
    beta,se,df,ok=logistic_cr2_satt(X,data[:failure],data[:person])
    truth=0.25scale; crit=tcrit975(df)
    (truth=truth,estimate=beta[2],se=se[2],covered=(beta[2]-crit*se[2]<=truth<=beta[2]+crit*se[2]),converged=ok)
end

function wcls_fit(data, outcome)
    A=data[:assign]; centered=A .- 0.5
    plan,_=standardize(data[:plan]); state,_=standardize(data[:belief])
    X=hcat(ones(length(A)),centered,plan,state,data[:previous])
    est,se,df=ols_cr2_satt(X,outcome,data[:person],focal=2)
    truth=mean1(data[:true_itt]); crit=tcrit975(df)
    (truth=truth,estimate=est,se=se,df=df,covered=(est-crit*se<=truth<=est+crit*se),
     detected=(est-crit*se>0),rejected=(abs(est/se)>crit))
end

function dynamic_nll(par, baseline, ymat, dmat)
    b0,bstate,bdiff,gamma=par
    alpha=0.80sigmoid(gamma)
    state=copy(baseline); value=0.0
    for t in 1:size(ymat,2), i in 1:size(ymat,1)
        eta=clip1(b0+bstate*state[i]+bdiff*dmat[i,t],-30,30)
        p=sigmoid(eta); y=ymat[i,t]
        value -= y*log(max(p,1e-12))+(1-y)*log(max(1-p,1e-12))
        state[i]=clip1(state[i]+alpha*(y-p)+0.04,-4.5,4.5)
    end
    value
end

function spaced_dynamic_nll(par, z0, ymat, amat, observed)
    bdelay,bstate,bprevious,gamma=par
    retention=sigmoid(gamma)
    state=copy(z0); previous=zeros(Float64,length(z0)); value=0.0
    for t in 1:size(ymat,2), i in 1:size(ymat,1)
        observed[i,t] || continue
        eta=clip1(-2.40+bdelay*amat[i,t]+bstate*state[i]+bprevious*previous[i],-30,30)
        p=sigmoid(eta); y=ymat[i,t]
        value -= y*log(max(p,1e-12))+(1-y)*log(max(1-p,1e-12))
        success=1-y
        state[i]=clip1(retention*state[i]+0.45*success-0.65*y,-4.5,4.5)
        previous[i]=y
    end
    value
end

function nelder_mead(f, x0; step=0.12, maxiter=350, tol=1e-7)
    n=length(x0); simplex=[copy(x0) for _ in 1:n+1]
    for i in 1:n simplex[i+1][i]+=step end
    values=[f(x) for x in simplex]
    for _ in 1:maxiter
        order=sortperm(values); simplex=simplex[order]; values=values[order]
        maximum(abs.(values[2:end] .- values[1])) < tol && break
        centroid=sum(simplex[1:n])/n
        xr=centroid+(centroid-simplex[end]); fr=f(xr)
        if fr < values[1]
            xe=centroid+2(xr-centroid); fe=f(xe)
            simplex[end],values[end]=(fe<fr ? xe : xr),(fe<fr ? fe : fr)
        elseif fr < values[n]
            simplex[end],values[end]=xr,fr
        else
            xc=centroid+0.5(simplex[end]-centroid); fc=f(xc)
            if fc < values[end]
                simplex[end],values[end]=xc,fc
            else
                for i in 2:n+1
                    simplex[i]=simplex[1]+0.5(simplex[i]-simplex[1]); values[i]=f(simplex[i])
                end
            end
        end
    end
    order=sortperm(values); simplex[order[1]],values[order[1]]
end

function numeric_hessian(f,x; h=2.0e-4)
    n=length(x); H=zeros(n,n); fx=f(x)
    for i in 1:n
        ei=zeros(n); ei[i]=h
        H[i,i]=(f(x+ei)-2fx+f(x-ei))/h^2
        for j in i+1:n
            ej=zeros(n); ej[j]=h
            H[i,j]=(f(x+ei+ej)-f(x+ei-ej)-f(x-ei+ej)+f(x-ei-ej))/(4h^2)
            H[j,i]=H[i,j]
        end
    end
    H
end

function tutor_dynamic_parameter_recovery(seed, adaptivity)
    # The learner state starts at an observed calibrated baseline.  The tutor
    # belief is separate and drives selection through observed history.
    data,meta=simulate_tutor(seed,adaptivity,n_people=30,n_skills=6,n_opps=8,
                             calibrated_initial=true)
    f=x->dynamic_nll(x,meta.baseline,meta.ymat,meta.dmat)
    x0=[0.30,0.95,-0.95,log(0.30/(0.80-0.30))]
    fit,_=nelder_mead(f,x0)
    H=numeric_hessian(f,fit)
    cov=pinv(H;rtol=1e-9)
    alpha=0.80sigmoid(fit[4]); dalpha=alpha*(1-alpha/0.80)
    estimates=[fit[1],fit[2],fit[3],alpha]
    ses=[sqrt(max(cov[1,1],0)),sqrt(max(cov[2,2],0)),sqrt(max(cov[3,3],0)),
         abs(dalpha)*sqrt(max(cov[4,4],0))]
    truths=[0.35,1.0,-1.0,0.30]
    labels=["intercept","state_slope","difficulty_slope","update_gain"]
    [(parameter=labels[i],truth=truths[i],estimate=estimates[i],se=ses[i],
      covered=(estimates[i]-1.96ses[i]<=truths[i]<=estimates[i]+1.96ses[i]),
      observed_rate=1.0) for i in 1:4]
end

function spaced_dynamic_parameter_recovery(seed, adaptivity;
                                           context_strength=0.0,
                                           observation_strength=0.0)
    meta=simulate_spaced_dynamic(seed,adaptivity,n_people=30,n_items=20,n_opps=8,
                                 context_strength=context_strength,
                                 observation_strength=observation_strength)
    f=x->spaced_dynamic_nll(x,meta.z0,meta.ymat,meta.amat,meta.observed)
    x0=[0.20,-0.50,0.30,log(0.70/(1-0.70))]
    fit,_=nelder_mead(f,x0,step=0.10,maxiter=450)
    H=numeric_hessian(f,fit)
    cov=pinv(H;rtol=1e-9)
    retention=sigmoid(fit[4]); dret=retention*(1-retention)
    estimates=[fit[1],fit[2],fit[3],retention]
    ses=[sqrt(max(cov[1,1],0)),sqrt(max(cov[2,2],0)),sqrt(max(cov[3,3],0)),
         abs(dret)*sqrt(max(cov[4,4],0))]
    truths=[0.25,-0.55,0.35,0.75]
    labels=["delay_slope","state_slope","previous_failure","retention"]
    [(parameter=labels[i],truth=truths[i],estimate=estimates[i],se=ses[i],
      covered=(estimates[i]-1.96ses[i]<=truths[i]<=estimates[i]+1.96ses[i]),
      observed_rate=meta.observed_rate) for i in 1:4]
end

function state_history(initial,ymat,dmat,family,parameter; misspecified=false)
    state=copy(initial); variance=fill(0.80,length(initial)); hist=zeros(size(ymat))
    intercept=misspecified ? 0.60 : 0.35; slope=misspecified ? 0.80 : 1.0
    for t in 1:size(ymat,2)
        hist[:,t]=state
        p=sigmoid.(intercept .+ slope.*state .- dmat[:,t])
        if family=="RL"
            state=clamp.(state+parameter.*(ymat[:,t]-p).+0.04,-4.5,4.5)
        else
            pred=variance.+parameter
            post=1 ./ (1 ./ pred + max.(p.*(1 .- p),1e-5))
            state=clamp.(state+post.*(ymat[:,t]-p).+0.04,-4.5,4.5); variance=post
        end
    end
    hist
end

function recovery_choice(meta; misspecified=false)
    intercept=misspecified ? 0.60 : 0.35; slope=misspecified ? 0.80 : 1.0
    grids=Dict("RL"=>[.10,.20,.30,.40,.50],"Bayesian"=>[.03,.08,.15,.25,.40])
    best=Dict{String,Float64}()
    for family in keys(grids)
        values=Float64[]
        for par in grids[family]
            h=state_history(meta.baseline,meta.ymat,meta.dmat,family,par,misspecified=misspecified)
            p=clamp.(sigmoid.(intercept .+ slope.*h .- meta.dmat),1e-9,1-1e-9)
            push!(values,-sum(meta.ymat.*log.(p)+(1 .- meta.ymat).*log.(1 .- p)))
        end
        best[family]=minimum(values)
    end
    best["RL"]<=best["Bayesian"] ? "RL" : "Bayesian"
end

function write_namedtuples(path,rows)
    isempty(rows) && return
    names=propertynames(rows[1])
    open(path,"w") do io
        println(io,join(string.(names),","))
        for row in rows
            println(io,csv_row((getproperty(row,n) for n in names)...))
        end
    end
end

function freeze_dataset(path,data)
    keys0=sort(collect(keys(data)),by=string)
    open(path,"w") do io
        println(io,join(string.(keys0),","))
        for i in eachindex(data[keys0[1]])
            println(io,csv_row((data[k][i] for k in keys0)...))
        end
    end
end

function run_observational(; architectures=("spaced","tutor"), suffix="")
    rows=NamedTuple[]
    for architecture in architectures, a in (0.0,0.40,0.70,1.00)
        reps=a==0.70 ? 1000 : 400
        for r in 0:reps-1
            seed=MASTER_SEED+(architecture=="spaced" ? 100000 : 1000000)+Int(round(1000a))*2000+r
            data=architecture=="spaced" ? simulate_spaced(seed,30,a) :
                 simulate_tutor(seed,a,n_people=30,n_skills=6,n_opps=8)[1]
            fits,meta=observational_fits(data,architecture)
            for fit in fits
                push!(rows,(architecture=architecture,adaptivity=a,replication=r,
                            prediction=meta.prediction,outcome_rate=meta.outcome_rate,
                            plan_actual_corr=meta.plan_actual_corr,action_state_corr=meta.action_state_corr,
                            estimator=fit.estimator,truth=fit.truth,estimate=fit.estimate,se=fit.se,
                            df=fit.df,covered=fit.covered,sign_error=fit.sign_error,converged=fit.converged))
            end
            if a==0.70 && r<5
                freeze_dataset(joinpath(VALID,"$(architecture)_primary_rep$(r).csv"),data)
            end
        end
    end
    write_namedtuples(joinpath(OUT,"observational_replications$(suffix).csv"),rows)
end

function run_oracle_and_context()
    oracle=NamedTuple[]; context=NamedTuple[]
    for r in 0:499
        data=simulate_spaced(MASTER_SEED+5000000+r,30,0.70,oracle_policy=true)
        fits,meta=observational_fits(data,"spaced")
        for fit in fits
            push!(oracle,(replication=r,estimator=fit.estimator,truth=fit.truth,
                          estimate=fit.estimate,se=fit.se,df=fit.df,covered=fit.covered,
                          sign_error=fit.sign_error,converged=fit.converged))
        end
    end
    for strength in (0.0,0.5,1.0), r in 0:399
        data=simulate_spaced(MASTER_SEED+6000000+Int(1000strength)*1000+r,30,0.70,context_strength=strength)
        fit=context_fit(data)
        push!(context,(context_strength=strength,proxy_r2=NaN,replication=r,
                       truth=fit.truth,estimate=fit.estimate,se=fit.se,
                       covered=fit.covered,converged=fit.converged))
    end
    for r2 in (0.0,0.25,0.50,0.75,1.0), r in 0:399
        data=simulate_spaced(MASTER_SEED+7000000+Int(1000r2)*1000+r,30,0.70,
                             context_strength=1.0,context_proxy_r2=r2)
        fit=context_fit(data,include_proxy=true)
        push!(context,(context_strength=1.0,proxy_r2=r2,replication=r,
                       truth=fit.truth,estimate=fit.estimate,se=fit.se,
                       covered=fit.covered,converged=fit.converged))
    end
    write_namedtuples(joinpath(OUT,"oracle_policy_replications.csv"),oracle)
    write_namedtuples(joinpath(OUT,"context_sensitivity_replications.csv"),context)
end

function run_mrt(; architectures=("spaced","tutor"), suffix="")
    effect=NamedTuple[]; nulls=NamedTuple[]; compliance_rows=NamedTuple[]
    for arch in architectures, n in (12,30,60,100), r in 0:999
        seed=MASTER_SEED+(arch=="spaced" ? 3000000 : 4000000)+n*10000+r
        if arch=="spaced"
            data=simulate_spaced(seed,n,0.70,context_strength=1.0,randomized_delta=0.30)
            fit=wcls_fit(data,data[:failure])
        else
            data=simulate_tutor(seed,0.70,n_people=n,randomized_delta=0.25)[1]
            fit=wcls_fit(data,1 .- data[:success])
        end
        push!(effect,(architecture=arch,n_people=n,replication=r,truth=fit.truth,
                      estimate=fit.estimate,se=fit.se,df=fit.df,covered=fit.covered,
                      detected=fit.detected,rejected=fit.rejected))
        if n==30 && r<3 freeze_dataset(joinpath(VALID,"$(arch)_mrt_rep$(r).csv"),data) end
    end
    for arch in architectures, n in (12,30), r in 0:999
        seed=MASTER_SEED+(arch=="spaced" ? 4500000 : 4700000)+n*10000+r
        if arch=="spaced"
            data=simulate_spaced(seed,n,0.70,context_strength=1.0,force_randomize=true)
            fit=wcls_fit(data,data[:failure])
        else
            data=simulate_tutor(seed,0.70,n_people=n,force_randomize=true)[1]
            fit=wcls_fit(data,1 .- data[:success])
        end
        push!(nulls,(architecture=arch,n_people=n,replication=r,estimate=fit.estimate,
                     se=fit.se,df=fit.df,rejected=fit.rejected))
    end
    for arch in architectures, r in 0:999
        seed=MASTER_SEED+(arch=="spaced" ? 8000000 : 8500000)+r
        if arch=="spaced"
            data=simulate_spaced(seed,30,0.70,context_strength=1.0,randomized_delta=0.30,compliance=0.70)
            fit=wcls_fit(data,data[:failure])
        else
            data=simulate_tutor(seed,0.70,n_people=30,randomized_delta=0.25,compliance=0.70)[1]
            fit=wcls_fit(data,1 .- data[:success])
        end
        push!(compliance_rows,(architecture=arch,n_people=30,compliance=0.70,
                               replication=r,truth=fit.truth,estimate=fit.estimate,se=fit.se,
                               df=fit.df,covered=fit.covered,detected=fit.detected))
    end
    write_namedtuples(joinpath(OUT,"mrt_effect_replications$(suffix).csv"),effect)
    write_namedtuples(joinpath(OUT,"mrt_null_replications$(suffix).csv"),nulls)
    write_namedtuples(joinpath(OUT,"mrt_compliance_replications$(suffix).csv"),compliance_rows)
end

function run_dynamic_recovery(; adaptivities=(0.0,0.70), suffix="")
    rows=NamedTuple[]
    for architecture in ("spaced","tutor"), a in adaptivities, r in 0:399
        seed=MASTER_SEED+9000000+(architecture=="spaced" ? 0 : 2000000)+Int(1000a)*1000+r
        fits=architecture=="spaced" ? spaced_dynamic_parameter_recovery(seed,a) :
                                      tutor_dynamic_parameter_recovery(seed,a)
        for fit in fits
            push!(rows,(architecture=architecture,mechanism="clean",adaptivity=a,
                        context_strength=0.0,observation_strength=0.0,
                        replication=r,parameter=fit.parameter,truth=fit.truth,
                        estimate=fit.estimate,se=fit.se,covered=fit.covered,
                        observed_rate=fit.observed_rate))
        end
    end
    write_namedtuples(joinpath(OUT,"dynamic_parameter_recovery_replications$(suffix).csv"),rows)
end

function run_observation_selection_stress()
    rows=NamedTuple[]
    conditions=((0.0,0.0,"clean"),(0.5,0.0,"outcome_omission"),
                (1.0,0.0,"outcome_omission"),
                (1.0,0.5,"informative_observation"),(1.0,1.0,"informative_observation"))
    for (context_strength,observation_strength,mechanism) in conditions, r in 0:399
        seed=MASTER_SEED+13000000+Int(1000context_strength)*10000+
             Int(1000observation_strength)*1000+r
        fits=spaced_dynamic_parameter_recovery(seed,0.70,
                                               context_strength=context_strength,
                                               observation_strength=observation_strength)
        for fit in fits
            push!(rows,(architecture="spaced",mechanism=mechanism,adaptivity=0.70,
                        context_strength=context_strength,
                        observation_strength=observation_strength,
                        replication=r,parameter=fit.parameter,truth=fit.truth,
                        estimate=fit.estimate,se=fit.se,covered=fit.covered,
                        observed_rate=fit.observed_rate))
        end
    end
    write_namedtuples(joinpath(OUT,"observation_selection_stress_replications.csv"),rows)
end

function run_update_class()
    rows=NamedTuple[]
    for a in (0.0,0.70), family in ("RL","Bayesian"), miss in (false,true), r in 0:399
        seed=MASTER_SEED+10000000+Int(1000a)*10000+(family=="RL" ? 0 : 4000)+(miss ? 2000 : 0)+r
        _,meta=simulate_tutor(seed,a,n_people=30,n_skills=8,n_opps=4,family=family)
        selected=recovery_choice(meta,misspecified=miss)
        push!(rows,(adaptivity=a,true_family=family,
                    observation_model=(miss ? "misspecified" : "correct"),replication=r,
                    selected_family=selected,recovered=(selected==family)))
    end
    write_namedtuples(joinpath(OUT,"update_class_replications.csv"),rows)
end

function main()
    run_observational(architectures=("spaced",),suffix="_spaced")
    run_observational(architectures=("tutor",),suffix="_tutor")
    run_oracle_and_context()
    run_mrt(architectures=("spaced",),suffix="_spaced")
    run_mrt(architectures=("tutor",),suffix="_tutor")
    run_dynamic_recovery()
    run_observation_selection_stress()
    run_update_class()
    open(joinpath(OUT,"run_manifest.txt"),"w") do io
        println(io,"master_seed=$(MASTER_SEED)")
        println(io,"julia_version=$(VERSION)")
        println(io,"primary_observational_replications=1000")
        println(io,"adaptivity_grid_nonprimary_replications=400")
        println(io,"mrt_replications_per_cell=1000")
        println(io,"dynamic_recovery_replications_per_cell=400")
        println(io,"observation_selection_stress_replications_per_cell=400")
        println(io,"update_class_replications_per_cell=400")
    end
    println("Julia simulations completed: $(OUT)")
end

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
