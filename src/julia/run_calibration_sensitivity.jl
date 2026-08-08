include(joinpath(@__DIR__, "run_closed_loop_simulation.jl"))

# Empirical signatures are reconstructed from the deidentified Anki transition
# table.  The loss scales were fixed before the search so that a practically
# meaningful mismatch in any one signature contributes comparably.
const CALIBRATION_TARGET_FAILURE = 0.09731051344743276
const CALIBRATION_TARGET_SPEARMAN = 0.8146452768576546
const CALIBRATION_TARGET_NAIVE_LOG_OR = -0.2158255188692115
const CALIBRATION_SCALE_FAILURE = 0.015
const CALIBRATION_SCALE_SPEARMAN = 0.05
const CALIBRATION_SCALE_NAIVE_LOG_OR = 0.10

const CALIBRATION_ADAPTIVITY_GRID = (0.50, 0.75, 1.00)
const CALIBRATION_DEVIATION_SD_GRID = (0.35, 0.50, 0.65)
const CALIBRATION_INTERCEPT_GRID = (-2.65, -2.40, -2.15)

const CALIBRATION_SEED = 20260731
const CALIBRATION_N_PARTICIPANTS = 12
const CALIBRATION_STATE_RELIABILITY = 0.80

function calibration_loss(failure_rate::Float64, scheduled_actual_spearman::Float64,
                          naive_log_odds::Float64)
    ((failure_rate - CALIBRATION_TARGET_FAILURE) /
     CALIBRATION_SCALE_FAILURE)^2 +
    ((scheduled_actual_spearman - CALIBRATION_TARGET_SPEARMAN) /
     CALIBRATION_SCALE_SPEARMAN)^2 +
    ((naive_log_odds - CALIBRATION_TARGET_NAIVE_LOG_OR) /
     CALIBRATION_SCALE_NAIVE_LOG_OR)^2
end

function calibration_dataset(seed::Int, adaptivity::Float64,
                             deviation_sd::Float64, outcome_intercept::Float64)
    condition = SimulationCondition(
        CALIBRATION_N_PARTICIPANTS,
        adaptivity=adaptivity,
        state_reliability=CALIBRATION_STATE_RELIABILITY,
    )
    simulate_closed_loop(
        seed, condition,
        natural_deviation_sd=deviation_sd,
        outcome_intercept=outcome_intercept,
    )
end

const SEARCH_COLUMNS = [
    "rank", "adaptivity", "natural_deviation_sd", "outcome_intercept",
    "n_replications", "mean_failure_rate", "target_failure_rate",
    "mean_scheduled_actual_spearman", "target_scheduled_actual_spearman",
    "mean_naive_log_odds", "target_naive_log_odds", "loss",
]

function run_calibration_search(output_dir::String; n_replications::Int=40,
                                seed_offset::Int=160_000_000)
    candidates = Any[]
    for adaptivity in CALIBRATION_ADAPTIVITY_GRID
        for deviation_sd in CALIBRATION_DEVIATION_SD_GRID
            for outcome_intercept in CALIBRATION_INTERCEPT_GRID
                failure_total = 0.0
                spearman_total = 0.0
                naive_total = 0.0
                println("[Julia calibration search] a=$(adaptivity), sd=$(deviation_sd), b0=$(outcome_intercept)")
                for replication in 0:n_replications - 1
                    # Common random numbers make candidate differences more
                    # precise and keep the ranking reproducible.
                    data, metrics = calibration_dataset(
                        CALIBRATION_SEED + seed_offset + replication,
                        adaptivity, deviation_sd, outcome_intercept,
                    )
                    naive = logistic_cluster_fit(
                        design_matrix(data, [:actual_gap_z]),
                        data[:failure], data[:participant],
                    )
                    failure_total += metrics.failure_rate
                    spearman_total += metrics.planned_actual_spearman
                    naive_total += naive.estimate
                end
                mean_failure = failure_total / n_replications
                mean_spearman = spearman_total / n_replications
                mean_naive = naive_total / n_replications
                push!(candidates, (
                    adaptivity=adaptivity,
                    natural_deviation_sd=deviation_sd,
                    outcome_intercept=outcome_intercept,
                    n_replications=n_replications,
                    mean_failure_rate=mean_failure,
                    mean_scheduled_actual_spearman=mean_spearman,
                    mean_naive_log_odds=mean_naive,
                    loss=calibration_loss(mean_failure, mean_spearman, mean_naive),
                ))
            end
        end
    end
    order = sortperm(collect(eachindex(candidates)), by=index -> candidates[index].loss)
    path = joinpath(output_dir, "calibration_search_candidates.csv")
    open(path, "w") do io
        write_csv_header(io, SEARCH_COLUMNS)
        for (rank, index) in enumerate(order)
            candidate = candidates[index]
            println(io, csv_row(
                rank, candidate.adaptivity, candidate.natural_deviation_sd,
                candidate.outcome_intercept, candidate.n_replications,
                candidate.mean_failure_rate, CALIBRATION_TARGET_FAILURE,
                candidate.mean_scheduled_actual_spearman,
                CALIBRATION_TARGET_SPEARMAN, candidate.mean_naive_log_odds,
                CALIBRATION_TARGET_NAIVE_LOG_OR, candidate.loss,
            ))
        end
    end
    [candidates[index] for index in order]
end

const ROBUSTNESS_COLUMNS = [
    "calibration_rank", "replication", "n_participants", "n_items", "n_stages",
    "adaptivity", "natural_deviation_sd", "outcome_intercept",
    "state_reliability", "failure_rate", "planned_actual_spearman",
    "estimator", "estimate", "standard_error", "ci_low", "ci_high",
    "converged", "truth", "bias", "squared_error", "covered", "sign_error",
]

function run_calibration_robustness(output_dir::String, ranked_candidates;
                                    n_top::Int=5, n_replications::Int=100,
                                    seed_offset::Int=170_000_000)
    path = joinpath(output_dir, "calibration_robustness_replications.csv")
    open(path, "w") do io
        write_csv_header(io, ROBUSTNESS_COLUMNS)
        for rank in 1:min(n_top, length(ranked_candidates))
            candidate = ranked_candidates[rank]
            println("[Julia calibration robustness] rank=$(rank), a=$(candidate.adaptivity), sd=$(candidate.natural_deviation_sd), b0=$(candidate.outcome_intercept)")
            for replication in 0:n_replications - 1
                # An independent stream is used for post-search robustness.
                data, metrics = calibration_dataset(
                    CALIBRATION_SEED + seed_offset + replication,
                    candidate.adaptivity, candidate.natural_deviation_sd,
                    candidate.outcome_intercept,
                )
                for result in fit_recovery_estimators(
                    data, metrics.true_delay_coefficient_per_sd,
                )
                    println(io, csv_row(
                        rank, replication, CALIBRATION_N_PARTICIPANTS,
                        DEFAULT_N_ITEMS, DEFAULT_N_STAGES,
                        candidate.adaptivity, candidate.natural_deviation_sd,
                        candidate.outcome_intercept,
                        CALIBRATION_STATE_RELIABILITY, metrics.failure_rate,
                        metrics.planned_actual_spearman, result.estimator,
                        result.estimate, result.standard_error,
                        result.ci_low, result.ci_high, result.converged,
                        result.truth, result.bias, result.squared_error,
                        result.covered, result.sign_error,
                    ))
                end
            end
        end
    end
    path
end

function main_calibration()
    project_dir = normpath(joinpath(@__DIR__, "..", ".."))
    output_dir = ensure_directory(get(
        ENV, "BRM_OUTPUT_DIR", joinpath(project_dir, "outputs"),
    ))
    search_replications = parse(Int, get(ENV, "BRM_CALIBRATION_SEARCH_REPS", "40"))
    robustness_replications = parse(Int, get(ENV, "BRM_CALIBRATION_ROBUST_REPS", "100"))
    search_seed_offset = parse(Int, get(
        ENV, "BRM_CALIBRATION_SEARCH_SEED_OFFSET", "160000000",
    ))
    robustness_seed_offset = parse(Int, get(
        ENV, "BRM_CALIBRATION_ROBUST_SEED_OFFSET", "170000000",
    ))
    ranked = run_calibration_search(
        output_dir, n_replications=search_replications,
        seed_offset=search_seed_offset,
    )
    run_calibration_robustness(
        output_dir, ranked, n_replications=robustness_replications,
        seed_offset=robustness_seed_offset,
    )
    println("Julia calibration-sensitivity simulations completed in $(output_dir)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_calibration()
end
