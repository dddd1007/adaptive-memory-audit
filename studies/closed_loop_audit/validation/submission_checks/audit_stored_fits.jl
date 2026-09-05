# Read-only audit of stored confirmatory fits. No simulation or refitting.
using TOML, SHA, Printf
length(ARGS)==2 || error("Usage: audit_stored_fits.jl STUDY_ROOT OUTPUT_DIRECTORY")
const study = abspath(ARGS[1])
const dest = abspath(ARGS[2])
mkpath(dest)
const settings = TOML.parsefile(joinpath(study, "config/execution.toml"))
const resultroot = joinpath(study, "outputs/confirmatory")
const frozen = TOML.parsefile(joinpath(study,"config/frozen_run_manifest.toml"))["files"]
@assert all(bytes2hex(sha256(read(joinpath(study,n))))==h for (n,h) in frozen)
rows = []
hashes = Dict{String,String}()
for lambda in 0:2, rep in 1:1000
    path = joinpath(resultroot, "lambda$lambda", @sprintf("rep%04d.toml",rep))
    hashes[relpath(path, study)] = bytes2hex(sha256(read(path)))
    f = TOML.parsefile(path)
    @assert f["replication"] == rep && f["lambda"] == lambda
    @assert f["seed"] == 220260905 + rep
    @assert f["code_manifest"] == frozen
    @assert f["aware"] == false && f["integration_method"] == "piecewise_legendre"
    flags = (gamma_boundary=abs(f["theta"][4]) >= settings["gamma_boundary"],
             coefficient_boundary=maximum(abs.(f["theta"][1:3])) >= settings["parameter_boundary"],
             gradient_failure=!(f["gradient_max"] < settings["gradient_tolerance"]),
             information_failure=!(f["hessian_min_eigenvalue"] > 0),
             conditioning_failure=!(f["hessian_condition"] < settings["maximum_hessian_condition"]),
             nonfinite=!(all(isfinite,f["theta"]) && isfinite(f["nll"])))
    @assert f["valid"] == !any(values(flags))
    push!(rows, (lambda=lambda,replication=rep,valid=f["valid"],status=f["status"],
        gamma=f["theta"][4],retention=f["estimate"][4], flags...))
end
open(joinpath(dest,"invalid_fit_diagnostics.csv"),"w") do io
    println(io, join(string.(keys(first(rows))),","))
    for r in rows
        r.valid || println(io, join(values(r),","))
    end
end
open(joinpath(dest,"fit_diagnostic_counts.csv"),"w") do io
    println(io,"lambda,total,valid,invalid,gamma_boundary,coefficient_boundary,gradient_failure,information_failure,conditioning_failure,nonfinite")
    for l in 0:2
        z=filter(r->r.lambda==l,rows)
        counts = [count(r->getproperty(r,k),z) for k in (:gamma_boundary,:coefficient_boundary,:gradient_failure,:information_failure,:conditioning_failure,:nonfinite)]
        vals = [l,length(z),count(r->r.valid,z),count(r->!r.valid,z),counts...]
        println(io,join(vals,",")); println(join(vals,","))
    end
end
open(joinpath(dest,"audited_result_hashes.toml"),"w") do io
    TOML.print(io,Dict("inputs"=>hashes);sorted=true)
end
println("Audited 3000 stored fits; no DGP, estimator, seed, or fit rule changed.")
