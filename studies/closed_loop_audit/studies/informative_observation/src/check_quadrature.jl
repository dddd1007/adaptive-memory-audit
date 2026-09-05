include("InformativeObservation.jl")
using .InformativeObservation, TOML, Printf, Test
c=config(); s=settings(); base=joinpath(studyroot,"outputs","confirmatory"); rows=NamedTuple[]
for lambda in 0:2
 fs=[TOML.parsefile(joinpath(base,"lambda$lambda",@sprintf("rep%04d.toml",rep))) for rep in 1:c["replications_per_condition"]]
 validids=findall(f->f["valid"],fs); invalidids=findall(f->!f["valid"],fs)
 worst=argmax([f["hessian_condition"] for f in fs]);worstvalid=validids[argmax([fs[j]["hessian_condition"] for j in validids])]
 ids=sort(unique([collect(1:c["estimation"]["preselected_node_check_replicates_per_condition"]);worst;worstvalid;invalidids]))
 for rep in ids
  f=fs[rep]; d=simulate(c,f["seed"],lambda,f["alpha"])
  check=fit(d,c;nodes=s["integration_nodes_check"])
  write_toml(joinpath(studyroot,"validation","quadrature",@sprintf("lambda%d_rep%04d.toml",lambda,rep)),check)
  H1=exact_hessian(f["theta"],d,c,gl(s["integration_nodes_main"]);step=s["hessian_step"])
  H2=exact_hessian(f["theta"],d,c,gl(s["integration_nodes_main"]);step=s["hessian_step_check"])
  hd=maximum(abs.(H1-H2))
  for j in 1:4
   de=abs(f["estimate"][j]-check["estimate"][j]); ds=abs(f["se"][j]-check["se"][j]); dn=abs(f["nll"]-check["nll"])
   equivalent=f["valid"]==check["valid"] && f["status"]==check["status"] && dn<=s["quadrature_nll_atol"] && de<=s["quadrature_parameter_atol"]+s["quadrature_parameter_rtol"]*abs(f["estimate"][j])
   se_pass= !f["valid"] || (ds<=s["quadrature_se_atol"]+s["quadrature_se_rtol"]*abs(f["se"][j]))
   pass=equivalent && se_pass && hd<s["hessian_check_atol"]
   push!(rows,(lambda=lambda,replication=rep,parameter=j,nodes_main=f["nodes"],nodes_check=check["nodes"],estimate_difference=de,se_difference=ds,nll_difference=dn,hessian_step_difference=hd,main_valid=f["valid"],check_valid=check["valid"],main_status=f["status"],check_status=check["status"],se_check_applicable=f["valid"],pass=pass))
  end
 end
end
write_csv(joinpath(studyroot,"validation","quadrature_comparison.csv"),rows)
println("Numerical stability checks: ",count(x->x.pass,rows),"/",length(rows),"; SE checks applicable: ",count(x->x.se_check_applicable,rows))
@test all(x->x.pass,rows)
