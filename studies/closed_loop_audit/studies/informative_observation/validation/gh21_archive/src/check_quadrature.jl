include("InformativeObservation.jl")
using .InformativeObservation, TOML, Printf, Test
c=config(); s=settings(); base=joinpath(studyroot,"outputs","confirmatory"); rows=NamedTuple[]
for lambda in 0:2
    fs=[TOML.parsefile(joinpath(base,"lambda$lambda",@sprintf("rep%04d.toml",rep))) for rep in 1:c["replications_per_condition"]]
    worst=argmax([f["hessian_condition"] for f in fs])
    ids=sort(unique([collect(1:c["estimation"]["preselected_node_check_replicates_per_condition"]);worst]))
    for rep in ids
        f=fs[rep]; d=simulate(c,f["seed"],lambda,f["alpha"])
        check=fit(d,c;nodes=c["estimation"]["gh_validation_nodes"])
        write_toml(joinpath(studyroot,"validation","quadrature",@sprintf("lambda%d_rep%04d.toml",lambda,rep)),check)
        for j in 1:4
            de=abs(f["estimate"][j]-check["estimate"][j]); ds=abs(f["se"][j]-check["se"][j]); dn=abs(f["nll"]-check["nll"])
            pass=f["valid"]==check["valid"] && f["valid"] && de<=s["quadrature_parameter_atol"]+s["quadrature_parameter_rtol"]*abs(f["estimate"][j]) && ds<=s["quadrature_se_atol"]+s["quadrature_se_rtol"]*abs(f["se"][j]) && dn<=s["quadrature_nll_atol"]
            push!(rows,(lambda=lambda,replication=rep,parameter=j,nodes_main=f["nodes"],nodes_check=check["nodes"],estimate_difference=de,se_difference=ds,nll_difference=dn,main_valid=f["valid"],check_valid=check["valid"],pass=pass))
        end
    end
end
write_csv(joinpath(studyroot,"validation","quadrature_comparison.csv"),rows)
println("Quadrature checks: ",count(x->x.pass,rows),"/",length(rows))
@test all(x->x.pass,rows)
