include("InformativeObservation.jl")
using .InformativeObservation, Statistics, TOML
c=config(); s=settings(); dest=joinpath(studyroot,"outputs","calibration"); mkpath(dest)
seed=c["master_seed"]+c["calibration_seed_offset"]+s["calibration_replication_id"]
holdoutseed=c["master_seed"]+c["calibration_seed_offset"]+s["holdout_replication_id"]
target=c["target_observation_rate"]; rows=NamedTuple[]; conditions=Dict[]
for lambda in c["lambda"]
    rate(a)=mean(simulate(c,seed,lambda,a;nseq=s["calibration_sequences"]).observed)
    grid=collect(-6.0:0.5:6.0); rates=rate.(grid)
    append!(rows,[(lambda=lambda,alpha=a,observed_rate=r,stage="scan") for (a,r) in zip(grid,rates)])
    if any(diff(rates).< -0.002); error("Calibration curve substantially nonmonotone; inspect scan") end
    left=findlast(<(target),rates); isnothing(left)&&error("Target not bracketed")
    left==length(grid)&&error("Target not bracketed")
    lo=grid[left]; hi=grid[left+1]
    if lambda==0
        alpha=log(target/(1-target))
    else
        for iter in 1:24
            mid=(lo+hi)/2
            if rate(mid)<target; lo=mid else hi=mid end
        end
        alpha=(lo+hi)/2
    end
    training=rate(alpha)
    holdout=mean(simulate(c,holdoutseed,lambda,alpha;nseq=s["holdout_sequences"]).observed)
    push!(conditions,Dict("lambda"=>lambda,"alpha"=>alpha,"training_rate"=>training,"holdout_rate"=>holdout))
    println("lambda=",lambda," alpha=",alpha," training=",training," holdout=",holdout); flush(stdout)
    abs(holdout-target)<=c["calibration_holdout_tolerance"]||error("Calibration holdout failed")
end
write_csv(joinpath(dest,"calibration_scan.csv"),rows)
write_toml(joinpath(studyroot,"config","calibrated_design.toml"),Dict("status"=>"calibrated_on_independent_data",
 "calibration_seed"=>seed,"holdout_seed"=>holdoutseed,"calibration_sequences"=>s["calibration_sequences"],
 "holdout_sequences"=>s["holdout_sequences"],"conditions"=>conditions,"design_sha256"=>filehash(joinpath(studyroot,"config","design.toml"))))
