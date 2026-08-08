include(joinpath(@__DIR__, "BRMCore.jl"))
using .BRMCore

const SEED = 20260731
const BETA_DELAY_RAW = 0.25
const DEFAULT_N_ITEMS = 100
const DEFAULT_N_STAGES = 2

struct SimulationCondition
    n_participants::Int
    adaptivity::Float64
    state_reliability::Float64
    context_strength::Float64
    randomized_perturbation::Float64
    perturbation_compliance::Float64
end

SimulationCondition(n_participants::Int;
    adaptivity::Float64=0.75,
    state_reliability::Float64=0.80,
    context_strength::Float64=0.0,
    randomized_perturbation::Float64=0.0,
    perturbation_compliance::Float64=1.0,
) = SimulationCondition(
    n_participants,
    adaptivity,
    state_reliability,
    context_strength,
    randomized_perturbation,
    perturbation_compliance,
)

function design_matrix(data::Dict{Symbol, Any}, columns::Vector{Symbol})
    n = length(data[columns[1]])
    matrix = zeros(Float64, n, length(columns))
    for column_index in eachindex(columns)
        values = data[columns[column_index]]
        for row in 1:n
            matrix[row, column_index] = Float64(values[row])
        end
    end
    matrix
end

function simulate_closed_loop(
    seed::Int,
    condition::SimulationCondition;
    n_items::Int=DEFAULT_N_ITEMS,
    n_stages::Int=DEFAULT_N_STAGES,
    natural_deviation_sd::Float64=0.50,
    outcome_intercept::Float64=-2.40,
)
    rng = DeterministicRNG(seed)
    n_participants = condition.n_participants
    n_pairs = n_participants * n_items
    n_rows = n_pairs * n_stages
    participant = Vector{Int}(undef, n_rows)
    item = Vector{Int}(undef, n_rows)
    stage = Vector{Float64}(undef, n_rows)
    memory_state_all = Vector{Float64}(undef, n_rows)
    previous_failure_all = Vector{Float64}(undef, n_rows)
    context_all = Vector{Float64}(undef, n_rows)
    planned_log_gap = Vector{Float64}(undef, n_rows)
    actual_log_gap = Vector{Float64}(undef, n_rows)
    assignment_all = Vector{Float64}(undef, n_rows)
    would_comply_all = Vector{Float64}(undef, n_rows)
    failure = Vector{Float64}(undef, n_rows)
    failure_probability = Vector{Float64}(undef, n_rows)

    learner_ability = [0.65 * rand_normal!(rng) for _ in 1:n_participants]
    item_difficulty = [0.65 * rand_normal!(rng) for _ in 1:n_items]
    memory_state = Vector{Float64}(undef, n_pairs)
    previous_failure = zeros(Float64, n_pairs)
    for pair in 1:n_pairs
        person = div(pair - 1, n_items) + 1
        card = mod(pair - 1, n_items) + 1
        memory_state[pair] = learner_ability[person] - item_difficulty[card] +
                             0.35 * rand_normal!(rng)
    end

    potential_effect_total = 0.0
    potential_effect_count = 0
    for current_stage in 0:n_stages - 1
        offset = current_stage * n_pairs
        next_memory = Vector{Float64}(undef, n_pairs)
        for pair in 1:n_pairs
            row = offset + pair
            person = div(pair - 1, n_items) + 1
            card = mod(pair - 1, n_items) + 1
            context = rand_normal!(rng)
            planned = 0.75 + condition.adaptivity * memory_state[pair] -
                      0.55 * previous_failure[pair] + 0.20 * rand_normal!(rng)
            natural_deviation = natural_deviation_sd * rand_normal!(rng) +
                                0.40 * condition.context_strength * context
            base_actual = planned + natural_deviation
            assignment = condition.randomized_perturbation > 0.0 ?
                         rand_sign!(rng) : 0.0
            would_comply = condition.randomized_perturbation > 0.0 ?
                           rand_bernoulli!(rng, condition.perturbation_compliance) : 0.0
            randomized_shift = condition.randomized_perturbation * assignment * would_comply
            actual = base_actual + randomized_shift
            base_predictor = outcome_intercept + BETA_DELAY_RAW * base_actual -
                             0.55 * memory_state[pair] +
                             0.35 * previous_failure[pair] +
                             0.50 * condition.context_strength * context
            probability = sigmoid(base_predictor + BETA_DELAY_RAW * randomized_shift)
            response = rand_bernoulli!(rng, probability)

            if condition.randomized_perturbation > 0.0
                probability_plus = sigmoid(
                    base_predictor + BETA_DELAY_RAW *
                    condition.randomized_perturbation * would_comply,
                )
                probability_minus = sigmoid(
                    base_predictor - BETA_DELAY_RAW *
                    condition.randomized_perturbation * would_comply,
                )
                potential_effect_total += probability_plus - probability_minus
                potential_effect_count += 1
            end

            participant[row] = person - 1
            item[row] = card - 1
            stage[row] = current_stage
            memory_state_all[row] = memory_state[pair]
            previous_failure_all[row] = previous_failure[pair]
            context_all[row] = context
            planned_log_gap[row] = planned
            actual_log_gap[row] = actual
            assignment_all[row] = assignment
            would_comply_all[row] = would_comply
            failure[row] = response
            failure_probability[row] = probability
            next_memory[pair] = 0.75 * memory_state[pair] +
                                0.45 * (1.0 - response) -
                                0.65 * response + 0.20 * rand_normal!(rng)
        end
        memory_state = next_memory
        for pair in 1:n_pairs
            previous_failure[pair] = failure[offset + pair]
        end
    end

    memory_state_z, _ = standardize(memory_state_all)
    observed_state = Vector{Float64}(undef, n_rows)
    if condition.state_reliability >= 1.0
        observed_state .= memory_state_z
    else
        error_scale = sqrt(
            (1.0 - condition.state_reliability) /
            max(condition.state_reliability, 1.0e-6),
        )
        for row in 1:n_rows
            observed_state[row] = memory_state_z[row] + error_scale * rand_normal!(rng)
        end
    end
    observed_state_z, _ = standardize(observed_state)
    actual_gap_z, actual_scale = standardize(actual_log_gap)
    planned_gap_z, _ = standardize(planned_log_gap)
    context_z, _ = standardize(context_all)
    stage_z, _ = standardize(stage)
    data = Dict{Symbol, Any}(
        :participant => participant,
        :item => item,
        :stage => stage,
        :memory_state => memory_state_all,
        :previous_failure => previous_failure_all,
        :context => context_all,
        :planned_log_gap => planned_log_gap,
        :actual_log_gap => actual_log_gap,
        :assignment => assignment_all,
        :would_comply => would_comply_all,
        :failure => failure,
        :failure_probability => failure_probability,
        :actual_gap_z => actual_gap_z,
        :planned_gap_z => planned_gap_z,
        :memory_state_z => memory_state_z,
        :observed_state_z => observed_state_z,
        :context_z => context_z,
        :stage_z => stage_z,
    )
    metrics = (
        true_delay_coefficient_per_sd=BETA_DELAY_RAW * actual_scale,
        failure_rate=mean1(failure),
        planned_actual_spearman=spearman_correlation(planned_log_gap, actual_log_gap),
        true_randomized_risk_difference=potential_effect_count > 0 ?
            potential_effect_total / potential_effect_count : NaN,
    )
    data, metrics
end

function fit_recovery_estimators(data::Dict{Symbol, Any}, truth::Float64)
    designs = [
        ("naive", [:actual_gap_z]),
        ("history", [:actual_gap_z, :observed_state_z, :previous_failure, :stage_z]),
        ("schedule_conditioned", [
            :actual_gap_z, :planned_gap_z, :observed_state_z,
            :previous_failure, :stage_z,
        ]),
        ("oracle", [
            :actual_gap_z, :memory_state_z, :context_z,
            :previous_failure, :stage_z,
        ]),
    ]
    rows = Any[]
    for (estimator, columns) in designs
        fit = logistic_cluster_fit(
            design_matrix(data, columns),
            data[:failure],
            data[:participant],
        )
        estimate = fit.estimate
        push!(rows, (
            estimator=estimator,
            estimate=estimate,
            standard_error=fit.standard_error,
            ci_low=fit.ci_low,
            ci_high=fit.ci_high,
            converged=fit.converged,
            truth=truth,
            bias=estimate - truth,
            squared_error=(estimate - truth)^2,
            covered=fit.ci_low <= truth <= fit.ci_high,
            sign_error=estimate < 0.0,
        ))
    end
    rows
end

const RECOVERY_COLUMNS = [
    "study", "replication", "n_participants", "n_items", "n_stages",
    "adaptivity", "state_reliability", "context_strength",
    "true_delay_coefficient_per_sd", "failure_rate",
    "planned_actual_spearman", "true_randomized_risk_difference",
    "estimator", "estimate", "standard_error", "ci_low", "ci_high",
    "converged", "truth", "bias", "squared_error", "covered", "sign_error",
]

function write_recovery_row(io, study, replication, condition, n_items, n_stages,
                            metrics, result)
    println(io, csv_row(
        study, replication, condition.n_participants, n_items, n_stages,
        condition.adaptivity, condition.state_reliability,
        condition.context_strength, metrics.true_delay_coefficient_per_sd,
        metrics.failure_rate, metrics.planned_actual_spearman,
        metrics.true_randomized_risk_difference, result.estimator,
        result.estimate, result.standard_error, result.ci_low, result.ci_high,
        result.converged, result.truth, result.bias, result.squared_error,
        result.covered, result.sign_error,
    ))
end

function run_parameter_recovery(output_dir::String; n_replications::Int=200)
    path = joinpath(output_dir, "parameter_recovery_replications.csv")
    open(path, "w") do io
        write_csv_header(io, RECOVERY_COLUMNS)
        for n_participants in (12, 30, 100)
            for adaptivity in (0.0, 0.40, 0.75, 1.10)
                for reliability in (0.50, 0.80)
                    condition = SimulationCondition(
                        n_participants,
                        adaptivity=adaptivity,
                        state_reliability=reliability,
                    )
                    println("[Julia recovery] N=$(n_participants), adaptivity=$(adaptivity), reliability=$(reliability)")
                    for replication in 0:n_replications - 1
                        seed = SEED + n_participants * 100_000 +
                               round(Int, adaptivity * 1_000) * 100 +
                               round(Int, reliability * 100) * 10 + replication
                        data, metrics = simulate_closed_loop(seed, condition)
                        results = fit_recovery_estimators(
                            data,
                            metrics.true_delay_coefficient_per_sd,
                        )
                        for result in results
                            write_recovery_row(
                                io, "parameter_recovery", replication,
                                condition, DEFAULT_N_ITEMS, DEFAULT_N_STAGES,
                                metrics, result,
                            )
                        end
                    end
                end
            end
        end
    end
    path
end

function run_context_violation(output_dir::String; n_replications::Int=160)
    path = joinpath(output_dir, "context_violation_replications.csv")
    open(path, "w") do io
        write_csv_header(io, RECOVERY_COLUMNS)
        for n_participants in (12, 30, 100)
            for context_strength in (0.0, 0.50, 1.00)
                condition = SimulationCondition(
                    n_participants,
                    adaptivity=0.75,
                    state_reliability=0.80,
                    context_strength=context_strength,
                )
                println("[Julia context] N=$(n_participants), context=$(context_strength)")
                for replication in 0:n_replications - 1
                    seed = SEED + 50_000_000 + n_participants * 100_000 +
                           round(Int, context_strength * 100) * 1_000 + replication
                    data, metrics = simulate_closed_loop(seed, condition)
                    results = fit_recovery_estimators(
                        data,
                        metrics.true_delay_coefficient_per_sd,
                    )
                    for result in results
                        write_recovery_row(
                            io, "unlogged_context", replication,
                            condition, DEFAULT_N_ITEMS, DEFAULT_N_STAGES,
                            metrics, result,
                        )
                    end
                end
            end
        end
    end
    path
end

const MRT_COLUMNS = [
    "study", "replication", "n_participants", "n_items", "n_stages",
    "perturbation_log_days", "perturbation_compliance", "context_strength",
    "estimator", "truth", "estimate", "standard_error", "ci_low", "ci_high",
    "converged", "bias", "squared_error", "covered", "detected_positive",
    "failure_rate", "planned_actual_spearman",
]

function run_micro_randomized(output_dir::String; n_replications::Int=250)
    path = joinpath(output_dir, "micro_randomized_replications.csv")
    open(path, "w") do io
        write_csv_header(io, MRT_COLUMNS)
        for n_participants in (12, 30, 60, 100)
            for perturbation in (0.15, 0.30, 0.45)
                for compliance in (0.70, 1.00)
                    condition = SimulationCondition(
                        n_participants,
                        adaptivity=0.75,
                        state_reliability=0.80,
                        context_strength=1.00,
                        randomized_perturbation=perturbation,
                        perturbation_compliance=compliance,
                    )
                    println("[Julia MRT] N=$(n_participants), delta=$(perturbation), compliance=$(compliance)")
                    for replication in 0:n_replications - 1
                        seed = SEED + 90_000_000 + n_participants * 100_000 +
                               round(Int, perturbation * 1_000) * 100 +
                               round(Int, compliance * 100) * 10 + replication
                        data, metrics = simulate_closed_loop(seed, condition)
                        assignment_01 = (data[:assignment] .+ 1.0) ./ 2.0
                        unadjusted_data = Dict{Symbol, Any}(:assignment_01 => assignment_01)
                        unadjusted = ols_cluster_fit(
                            design_matrix(unadjusted_data, [:assignment_01]),
                            data[:failure], data[:participant],
                        )
                        adjusted_data = copy(data)
                        adjusted_data[:assignment_01] = assignment_01
                        adjusted = ols_cluster_fit(
                            design_matrix(adjusted_data, [
                                :assignment_01, :planned_gap_z, :observed_state_z,
                                :previous_failure, :stage_z,
                            ]),
                            data[:failure], data[:participant],
                        )
                        truth = metrics.true_randomized_risk_difference
                        for (label, fit) in (
                            ("unadjusted_itt", unadjusted),
                            ("covariate_adjusted_itt", adjusted),
                        )
                            estimate = fit.estimate
                            println(io, csv_row(
                                "micro_randomized_experiment", replication,
                                n_participants, DEFAULT_N_ITEMS, DEFAULT_N_STAGES,
                                perturbation, compliance, 1.0, label, truth,
                                estimate, fit.standard_error, fit.ci_low, fit.ci_high,
                                fit.converged, estimate - truth, (estimate - truth)^2,
                                fit.ci_low <= truth <= fit.ci_high,
                                fit.ci_low > 0.0, metrics.failure_rate,
                                metrics.planned_actual_spearman,
                            ))
                        end
                    end
                end
            end
        end
    end
    path
end

const WCLS_COLUMNS = [
    "study", "replication", "n_participants", "n_items", "n_stages",
    "perturbation_log_days", "perturbation_compliance", "context_strength",
    "estimator", "truth", "estimate", "standard_error",
    "normal_ci_low", "normal_ci_high", "participant_t_ci_low",
    "participant_t_ci_high", "converged", "bias", "squared_error",
    "normal_covered", "participant_t_covered", "normal_detected_positive",
    "participant_t_detected_positive", "failure_rate",
    "planned_actual_spearman", "lpm_estimate", "lpm_standard_error",
    "absolute_estimate_difference", "absolute_standard_error_difference",
]

"""
Run the fixed-probability WCLS audit on the same declared MRT grid used by the
additive linear-probability analysis.  Each simulated data set contributes an
unadjusted and a pre-treatment-covariate-adjusted result.  The corresponding
LPM fit is retained in every row so numerical equivalence can be audited rather
than asserted.  Participant-*t* intervals use `N - 1` degrees of freedom.
"""
function run_micro_randomized_wcls(output_dir::String; n_replications::Int=250)
    path = joinpath(output_dir, "micro_randomized_wcls_replications.csv")
    open(path, "w") do io
        write_csv_header(io, WCLS_COLUMNS)
        for n_participants in (12, 30, 60, 100)
            t_critical = student_t_critical_975(n_participants - 1)
            for perturbation in (0.15, 0.30, 0.45)
                for compliance in (0.70, 1.00)
                    condition = SimulationCondition(
                        n_participants,
                        adaptivity=0.75,
                        state_reliability=0.80,
                        context_strength=1.00,
                        randomized_perturbation=perturbation,
                        perturbation_compliance=compliance,
                    )
                    println("[Julia WCLS] N=$(n_participants), delta=$(perturbation), compliance=$(compliance)")
                    for replication in 0:n_replications - 1
                        # Match the declared primary MRT stream: equivalence is
                        # checked on identical simulated decision points.
                        seed = SEED + 90_000_000 + n_participants * 100_000 +
                               round(Int, perturbation * 1_000) * 100 +
                               round(Int, compliance * 100) * 10 + replication
                        data, metrics = simulate_closed_loop(seed, condition)
                        assignment_01 = (data[:assignment] .+ 1.0) ./ 2.0
                        nuisance_unadjusted = zeros(Float64, length(data[:failure]), 0)
                        nuisance_adjusted = design_matrix(data, [
                            :planned_gap_z, :observed_state_z,
                            :previous_failure, :stage_z,
                        ])
                        for (label, nuisance) in (
                            ("unadjusted_itt", nuisance_unadjusted),
                            ("covariate_adjusted_itt", nuisance_adjusted),
                        )
                            wcls = wcls_cluster_fit(
                                assignment_01, nuisance, data[:failure],
                                data[:participant], randomization_probability=0.5,
                            )
                            lpm_design = zeros(Float64, length(data[:failure]),
                                               size(nuisance, 2) + 1)
                            lpm_design[:, 1] .= assignment_01
                            if size(nuisance, 2) > 0
                                lpm_design[:, 2:end] .= nuisance
                            end
                            lpm = ols_cluster_fit(
                                lpm_design, data[:failure], data[:participant],
                                focal_index=1,
                            )
                            truth = metrics.true_randomized_risk_difference
                            estimate = wcls.estimate
                            standard_error = wcls.standard_error
                            normal_low = estimate - 1.96 * standard_error
                            normal_high = estimate + 1.96 * standard_error
                            t_low = estimate - t_critical * standard_error
                            t_high = estimate + t_critical * standard_error
                            println(io, csv_row(
                                "micro_randomized_wcls", replication,
                                n_participants, DEFAULT_N_ITEMS, DEFAULT_N_STAGES,
                                perturbation, compliance, 1.0, label, truth,
                                estimate, standard_error, normal_low, normal_high,
                                t_low, t_high, wcls.converged, estimate - truth,
                                (estimate - truth)^2,
                                normal_low <= truth <= normal_high,
                                t_low <= truth <= t_high,
                                normal_low > 0.0, t_low > 0.0,
                                metrics.failure_rate,
                                metrics.planned_actual_spearman,
                                lpm.estimate, lpm.standard_error,
                                abs(estimate - lpm.estimate),
                                abs(standard_error - lpm.standard_error),
                            ))
                        end
                    end
                end
            end
        end
    end
    path
end

function orthogonal_context_proxy(
    context_z::Vector{Float64},
    noise::Vector{Float64},
    target_r2::Float64,
)
    0.0 <= target_r2 <= 1.0 || error("target_r2 must lie in [0, 1]")
    noise_z, _ = standardize(noise)
    projection_denominator = sum(value * value for value in context_z)
    projection_denominator > 0.0 || error("context has zero variance")
    projection = sum(context_z[index] * noise_z[index]
                     for index in eachindex(context_z)) /
                 projection_denominator
    residual = [noise_z[index] - projection * context_z[index]
                for index in eachindex(context_z)]
    residual_z, _ = standardize(residual)
    proxy = sqrt(target_r2) .* context_z .+
            sqrt(1.0 - target_r2) .* residual_z
    proxy_z, _ = standardize(proxy)
    proxy_z
end

const CONTEXT_PROXY_COLUMNS = [
    "study", "replication", "n_participants", "n_items", "n_stages",
    "context_strength", "proxy_r2_target", "proxy_r2_realized",
    "estimator", "truth", "estimate", "standard_error", "normal_ci_low",
    "normal_ci_high", "participant_t_ci_low", "participant_t_ci_high",
    "converged", "bias", "squared_error", "normal_covered",
    "participant_t_covered", "sign_error", "failure_rate",
    "planned_actual_spearman",
]

function fit_context_proxy_estimators(
    data::Dict{Symbol, Any}, truth::Float64, proxy::Vector{Float64},
)
    proxy_data = copy(data)
    proxy_data[:context_proxy_z] = proxy
    designs = [
        ("plan_conditioned", [
            :actual_gap_z, :planned_gap_z, :observed_state_z,
            :previous_failure, :stage_z,
        ]),
        ("plan_plus_proxy", [
            :actual_gap_z, :planned_gap_z, :observed_state_z,
            :previous_failure, :stage_z, :context_proxy_z,
        ]),
        ("oracle", [
            :actual_gap_z, :memory_state_z, :context_z,
            :previous_failure, :stage_z,
        ]),
    ]
    rows = Any[]
    for (estimator, columns) in designs
        fit = logistic_cluster_fit(
            design_matrix(proxy_data, columns), data[:failure],
            data[:participant],
        )
        estimate = fit.estimate
        push!(rows, (
            estimator=estimator,
            estimate=estimate,
            standard_error=fit.standard_error,
            ci_low=fit.ci_low,
            ci_high=fit.ci_high,
            converged=fit.converged,
            bias=estimate - truth,
            squared_error=(estimate - truth)^2,
            sign_error=estimate < 0.0,
        ))
    end
    rows
end

"""
Quantify how progressively better measurement of an otherwise unlogged context
repairs schedule-conditioned delay-effect recovery.  A sample-orthogonalized
noise component makes the realized squared correlation equal the declared
target (0, .25, .50, .75, or 1) to floating-point tolerance in every data set.
The data and orthogonal noise are held fixed across proxy-quality levels within
replication, making the quality gradient a paired Monte Carlo comparison.
"""
function run_context_proxy(output_dir::String; n_replications::Int=160)
    path = joinpath(output_dir, "context_proxy_replications.csv")
    n_participants = 30
    t_critical = student_t_critical_975(n_participants - 1)
    condition = SimulationCondition(
        n_participants,
        adaptivity=0.75,
        state_reliability=0.80,
        context_strength=1.00,
    )
    open(path, "w") do io
        write_csv_header(io, CONTEXT_PROXY_COLUMNS)
        for replication in 0:n_replications - 1
            data_seed = SEED + 70_000_000 + n_participants * 100_000 + replication
            data, metrics = simulate_closed_loop(data_seed, condition)
            noise_rng = DeterministicRNG(data_seed + 700_000_000)
            noise = [rand_normal!(noise_rng) for _ in eachindex(data[:context_z])]
            for target_r2 in (0.0, 0.25, 0.50, 0.75, 1.00)
                proxy = orthogonal_context_proxy(data[:context_z], noise, target_r2)
                realized_r2 = correlation(proxy, data[:context_z])^2
                results = fit_context_proxy_estimators(
                    data, metrics.true_delay_coefficient_per_sd, proxy,
                )
                println("[Julia context proxy] replication=$(replication), R2=$(target_r2)")
                for result in results
                    normal_low = result.estimate - 1.96 * result.standard_error
                    normal_high = result.estimate + 1.96 * result.standard_error
                    t_low = result.estimate - t_critical * result.standard_error
                    t_high = result.estimate + t_critical * result.standard_error
                    println(io, csv_row(
                        "context_proxy", replication, n_participants,
                        DEFAULT_N_ITEMS, DEFAULT_N_STAGES, 1.0,
                        target_r2, realized_r2, result.estimator,
                        metrics.true_delay_coefficient_per_sd,
                        result.estimate, result.standard_error,
                        normal_low, normal_high, t_low, t_high,
                        result.converged, result.bias, result.squared_error,
                        normal_low <= metrics.true_delay_coefficient_per_sd <= normal_high,
                        t_low <= metrics.true_delay_coefficient_per_sd <= t_high,
                        result.sign_error, metrics.failure_rate,
                        metrics.planned_actual_spearman,
                    ))
                end
            end
        end
    end
    path
end

function main()
    project_dir = normpath(joinpath(@__DIR__, "..", ".."))
    output_dir = ensure_directory(get(ENV, "BRM_OUTPUT_DIR", joinpath(project_dir, "outputs")))
    recovery_reps = parse(Int, get(ENV, "BRM_RECOVERY_REPS", "200"))
    context_reps = parse(Int, get(ENV, "BRM_CONTEXT_REPS", "160"))
    mrt_reps = parse(Int, get(ENV, "BRM_MRT_REPS", "250"))
    wcls_reps = parse(Int, get(ENV, "BRM_WCLS_REPS", "250"))
    proxy_reps = parse(Int, get(ENV, "BRM_PROXY_REPS", "160"))
    extensions_only = lowercase(get(ENV, "BRM_EXTENSIONS_ONLY", "false")) in
                      ("1", "true", "yes")
    if !extensions_only
        run_parameter_recovery(output_dir, n_replications=recovery_reps)
        run_context_violation(output_dir, n_replications=context_reps)
        run_micro_randomized(output_dir, n_replications=mrt_reps)
    end
    run_micro_randomized_wcls(output_dir, n_replications=wcls_reps)
    run_context_proxy(output_dir, n_replications=proxy_reps)
    println("Julia closed-loop simulations completed in $(output_dir)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
