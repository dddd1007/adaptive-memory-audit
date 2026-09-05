include("InformativeObservation.jl")
using .InformativeObservation, TOML, Printf
c=config(); rows=NamedTuple[]
for (lambda,rep) in ((2,2),(0,33),(1,30),(2,13))
    f=TOML.parsefile(joinpath(studyroot,"outputs","confirmatory","lambda$lambda",@sprintf("rep%04d.toml",rep)))
    d=simulate(c,f["seed"],lambda,f["alpha"])
    for nodes in (41,81,161)
        z=fit(d,c;nodes=nodes)
        write_toml(joinpath(studyroot,"validation","node_upgrade",@sprintf("lambda%d_rep%04d_q%d.toml",lambda,rep,nodes)),z)
        e=z["estimate"]
        push!(rows,(lambda=lambda,replication=rep,nodes=nodes,status=z["status"],nll=z["nll"],gradient=z["gradient_max"],e1=e[1],e2=e[2],e3=e[3],e4=e[4]))
        println(last(rows));flush(stdout)
    end
end
write_csv(joinpath(studyroot,"validation","node_upgrade_diagnostic.csv"),rows)
