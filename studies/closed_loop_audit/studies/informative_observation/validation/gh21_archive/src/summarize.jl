include("InformativeObservation.jl")
using .InformativeObservation, TOML, Statistics, Printf
c=config(); n=c["replications_per_condition"]; base=joinpath(studyroot,"outputs","confirmatory")
truth=[c["response"]["delay_coefficient"],c["response"]["state_coefficient"],c["response"]["previous_failure_coefficient"],c["learner"]["retention"]]
names=["delay_coefficient","state_coefficient","previous_failure_coefficient","retention"]
all=Dict(lambda=>[TOML.parsefile(joinpath(base,"lambda$lambda",@sprintf("rep%04d.toml",rep))) for rep in 1:n] for lambda in 0:2)
summary=NamedTuple[]; long=NamedTuple[]; contrasts=NamedTuple[]
for lambda in 0:2
    fs=all[lambda]; @assert length(fs)==n && length(unique(f["replication"] for f in fs))==n
    valid=[f["valid"] for f in fs]; nv=sum(valid); good=fs[valid]
    nv>1||error("Insufficient valid fits")
    for j in 1:4
        e=[f["estimate"][j] for f in good]; se=[f["se"][j] for f in good]
        covered=[f["valid"] && f["lower"][j]<=truth[j]<=f["upper"][j] for f in fs]
        pc=sum(covered)/nv; pu=sum(covered)/n
        push!(summary,(lambda=lambda,parameter=names[j],truth=truth[j],n_total=n,n_valid=nv,n_failed=n-nv,
          mean_estimate=mean(e),bias=mean(e)-truth[j],bias_mcse=std(e)/sqrt(nv),empirical_sd=std(e),mean_se=mean(se),
          se_sd_ratio=mean(se)/std(e),coverage_valid=pc,coverage_valid_mcse=sqrt(pc*(1-pc)/nv),
          coverage_unconditional=pu,coverage_unconditional_mcse=sqrt(pu*(1-pu)/n),
          mean_observed_rate=mean(f["observed_rate"] for f in fs)))
        for (rep,f) in enumerate(fs)
            push!(long,(lambda=lambda,replication=rep,parameter=names[j],truth=truth[j],estimate=f["estimate"][j],se=f["se"][j],lower=f["lower"][j],upper=f["upper"][j],valid=f["valid"],status=f["status"],covered=covered[rep],observed_rate=f["observed_rate"],seed=f["seed"],gradient_max=f["gradient_max"],hessian_condition=f["hessian_condition"]))
        end
    end
end
for lambda in 1:2, j in 1:4
    control=all[0]; treated=all[lambda]
    valid=[a["valid"]&&b["valid"] for (a,b) in zip(control,treated)]
    estimates=[treated[i]["estimate"][j]-control[i]["estimate"][j] for i in 1:n if valid[i]]
    coverage(f)=f["valid"]&&f["lower"][j]<=truth[j]<=f["upper"][j]
    diffcoverage=[Int(coverage(treated[i]))-Int(coverage(control[i])) for i in 1:n]
    paired=diffcoverage[valid]
    push!(contrasts,(lambda=lambda,parameter=names[j],paired_valid=sum(valid),n_total=n,bias_difference=mean(estimates),bias_difference_mcse=std(estimates)/sqrt(length(estimates)),coverage_difference_paired_valid=mean(paired),coverage_difference_paired_valid_mcse=std(paired)/sqrt(length(paired)),coverage_difference_unconditional=mean(diffcoverage),coverage_difference_unconditional_mcse=std(diffcoverage)/sqrt(n)))
end
write_csv(joinpath(base,"summary.csv"),summary)
write_csv(joinpath(base,"replication_results.csv"),long)
write_csv(joinpath(base,"paired_contrasts.csv"),contrasts)
foreach(println,summary)
