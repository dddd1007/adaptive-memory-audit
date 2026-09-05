using Distributed, TOML, SHA, Dates, LinearAlgebra, Printf
include("InformativeObservation.jl")
using .InformativeObservation
BLAS.set_num_threads(1)
const C=config(); const S=settings()
const CAL=TOML.parsefile(joinpath(studyroot,"config","calibrated_design.toml"))
function code_manifest()
    files=[joinpath(studyroot,"src",f) for f in ("InformativeObservation.jl","run_study.jl","calibrate.jl","test_likelihood.jl","test_piecewise.jl")]
    append!(files,[joinpath(studyroot,"config",f) for f in ("design.toml","execution.toml","calibrated_design.toml")])
    push!(files,joinpath(studyroot,"validation","independent_likelihood.R"))
    append!(files,[joinpath(studyroot,"src",f) for f in ("summarize.jl","check_quadrature.jl","export_validation.jl")])
    push!(files,joinpath(studyroot,"validation","summarize_and_plot.R"))
    Dict(relpath(f,studyroot)=>filehash(f) for f in files)
end
function run_one(job,mode,destination,manifest)
    condition,rep=job; lambda=condition["lambda"]; alpha=condition["alpha"]
    offset=mode=="smoke" ? C["validation_seed_offset"] : C["confirmatory_seed_offset"]
    seed=C["master_seed"]+offset+rep
    output=joinpath(destination,"lambda$(Int(lambda))",@sprintf("rep%04d.toml",rep))
    if isfile(output)
        old=TOML.parsefile(output)
        old["seed"]==seed&&old["code_manifest"]==manifest||error("Checkpoint mismatch: $output")
        return (lambda=lambda,replication=rep,status="retained",elapsed_seconds=0.0)
    end
    data=simulate(C,seed,lambda,alpha)
    result=try
        fit(data,C,S)
    catch err
        Dict("valid"=>false,"status"=>"exception","error"=>sprint(showerror,err),
          "theta"=>fill(NaN,4),"estimate"=>fill(NaN,4),"se"=>fill(NaN,4),"lower"=>fill(NaN,4),"upper"=>fill(NaN,4),
          "nll"=>NaN,"gradient_max"=>NaN,"hessian_condition"=>Inf,"hessian_min_eigenvalue"=>NaN,
          "nodes"=>get(S,"integration_nodes_main",C["estimation"]["gh_nodes"]),"aware"=>false,"observed_rate"=>sum(data.observed)/length(data.observed),"elapsed_seconds"=>NaN)
    end
    result["lambda"]=lambda; result["alpha"]=alpha; result["replication"]=rep; result["seed"]=seed
    result["stage"]=mode; result["code_manifest"]=manifest
    write_toml(output,result)
    if mode=="smoke" || rep<=3
        write_dataset(joinpath(destination,"frozen_data",@sprintf("lambda%d_rep%04d.csv",Int(lambda),rep)),data)
    end
    println("lambda=",lambda," rep=",rep," status=",result["status"]," seconds=",round(result["elapsed_seconds"],digits=2));flush(stdout)
    (lambda=lambda,replication=rep,status=result["status"],elapsed_seconds=result["elapsed_seconds"])
end
function main()
    mode=isempty(ARGS) ? "smoke" : ARGS[1]
    mode in ("smoke","benchmark","confirmatory")||error("Unknown mode")
    nw=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 1
    reps=mode=="smoke" ? (10001:10005) : mode=="benchmark" ? (1:12) : (1:C["replications_per_condition"])
    destination=joinpath(studyroot,"outputs",mode)
    manifest=code_manifest(); freeze=joinpath(studyroot,"config","frozen_run_manifest.toml")
    if mode!="smoke"
        isfile(freeze)||error("Freeze the design and code before formal data")
        TOML.parsefile(freeze)["files"]==manifest||error("Frozen code/design changed")
    end
    jobs=[(condition,rep) for rep in reps for condition in CAL["conditions"]]
    if nw>1
        addprocs(nw;exeflags=Cmd(["--startup-file=no","--threads=1"]))
        source=abspath(@__FILE__)
        @everywhere workers() begin
            include($source)
        end
    end
    started=now(); before=time()
    result=nw>1 ? pmap(job->run_one(job,mode,destination,manifest),jobs) : map(job->run_one(job,mode,destination,manifest),jobs)
    write_toml(joinpath(destination,"run_manifest.toml"),Dict("stage"=>mode,"started"=>string(started),"finished"=>string(now()),"wall_seconds"=>time()-before,"workers"=>nw,"blas_threads"=>1,"julia"=>string(VERSION),"code_manifest"=>manifest,"jobs"=>length(jobs)))
    println("Completed ",mode," jobs=",length(result)," seconds=",time()-before)
end
if abspath(PROGRAM_FILE)==@__FILE__; main() end
