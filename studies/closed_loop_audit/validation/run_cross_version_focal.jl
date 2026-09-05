#!/usr/bin/env julia

# Focal cross-version validation for the archived Julia implementation.
# Usage: julia --startup-file=no run_cross_version_focal.jl output.csv

include(joinpath(@__DIR__, "..", "src", "julia", "run_revision_simulations.jl"))

destination = length(ARGS) > 0 ? abspath(ARGS[1]) :
              joinpath(@__DIR__, "cross_version_focal.csv")
rows = NamedTuple[]

function addrow(case_name, metric, value)
    push!(rows, (julia_version=string(VERSION), case=case_name,
                 metric=metric, value=Float64(value)))
end

for architecture in ("spaced", "tutor")
    seed = MASTER_SEED + (architecture == "spaced" ? 0 : 1000000) + 700000 + 17
    if architecture == "spaced"
        data = simulate_spaced(seed, 30, 0.70)
    else
        data = simulate_tutor(seed, 0.70, n_people=30)[1]
    end
    fits, meta = observational_fits(data, architecture)
    for fit in fits
        addrow("$(architecture)_observational", "$(fit.estimator)_estimate", fit.estimate)
        addrow("$(architecture)_observational", "$(fit.estimator)_se", fit.se)
        addrow("$(architecture)_observational", "$(fit.estimator)_df", fit.df)
    end
    addrow("$(architecture)_observational", "conditional_truth", meta.truth)
end

for architecture in ("spaced", "tutor")
    seed = MASTER_SEED + 9000000 + (architecture == "spaced" ? 0 : 2000000) + 700000 + 17
    fits = architecture == "spaced" ? spaced_dynamic_parameter_recovery(seed, 0.70) :
                                       tutor_dynamic_parameter_recovery(seed, 0.70)
    for fit in fits
        addrow("$(architecture)_dynamic", "$(fit.parameter)_estimate", fit.estimate)
        addrow("$(architecture)_dynamic", "$(fit.parameter)_se", fit.se)
    end
end

for architecture in ("spaced", "tutor")
    seed = MASTER_SEED + (architecture == "spaced" ? 3000000 : 4000000) + 300000 + 17
    if architecture == "spaced"
        data = simulate_spaced(seed, 30, 0.70, context_strength=1.0, randomized_delta=0.30)
        fit = wcls_fit(data, data[:failure])
    else
        data = simulate_tutor(seed, 0.70, n_people=30, randomized_delta=0.25)[1]
        fit = wcls_fit(data, 1 .- data[:success])
    end
    addrow("$(architecture)_mrt", "estimate", fit.estimate)
    addrow("$(architecture)_mrt", "se", fit.se)
    addrow("$(architecture)_mrt", "df", fit.df)
    addrow("$(architecture)_mrt", "truth", fit.truth)
end

write_namedtuples(destination, rows)
println("Wrote focal cross-version results to $(destination)")
