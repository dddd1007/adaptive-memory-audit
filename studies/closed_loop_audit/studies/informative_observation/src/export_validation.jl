include("InformativeObservation.jl")
using .InformativeObservation, TOML, Printf
c=config(); stage=ARGS[1]; rows=NamedTuple[]
ids=stage=="smoke" ? [10001] : [1,2,3]
for lambda in (0,1,2), rep in ids
    dataset=joinpath("outputs",stage,"frozen_data",@sprintf("lambda%d_rep%04d.csv",lambda,rep))
    d=read_dataset(joinpath(studyroot,dataset))
    f=TOML.parsefile(joinpath(studyroot,"outputs",stage,"lambda$lambda",@sprintf("rep%04d.toml",rep)))
    for aware in (stage=="smoke" ? (false,true) : (false,))
        fitresult=aware ? fit(d,c;aware=true) : f
        if aware; write_toml(joinpath(studyroot,"validation",@sprintf("aware_lambda%d_rep%04d.toml",lambda,rep)),fitresult) end
        x=fitresult["theta"]; e=fitresult["estimate"]; se=fitresult["se"]
        push!(rows,(dataset=dataset,lambda=lambda,replication=rep,aware=aware,nodes=fitresult["nodes"],nll=fitresult["nll"],
          theta1=x[1],theta2=x[2],theta3=x[3],theta4=x[4],estimate1=e[1],estimate2=e[2],estimate3=e[3],estimate4=e[4],
          se1=se[1],se2=se[2],se3=se[3],se4=se[4],valid=fitresult["valid"]))
    end
end
write_csv(joinpath(studyroot,"validation","$(stage)_julia_refits.csv"),rows)
