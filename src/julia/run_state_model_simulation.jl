include(joinpath(@__DIR__, "BRMCore.jl"))
using .BRMCore

const STATE_SEED = 20260731
const STATE_N_PARTICIPANTS = 30
const STATE_N_ITEMS = 40
const BETA_DELAY = 0.25
const OUTCOME_INTERCEPT = -1.85
const STATE_COEFFICIENT = -0.65
const PREVIOUS_FAILURE_COEFFICIENT = 0.35
const SCHEDULER_ADAPTIVITY = 0.75
const RL_LEARNING_RATE = 0.25
const BAYES_VOLATILITY = 0.20
const RL_SCHEDULER_RATE = 0.20
const BAYES_SCHEDULER_VOLATILITY = 0.12
const INITIAL_BAYES_VARIANCE = 0.80
const RL_PARAMETER_GRID = (0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
const BAYES_PARAMETER_GRID = (0.02, 0.06, 0.12, 0.20, 0.30, 0.45, 0.65)

struct StateModelCondition
    learner_family::String
    scheduler_family::String
    n_stages::Int
end

struct ObservationCondition
    label::String
    intercept_shift::Float64
    state_slope_multiplier::Float64
end

const OBSERVATION_CONDITIONS = (
    ObservationCondition("correct", 0.00, 1.00),
    ObservationCondition("intercept_minus_0.25", -0.25, 1.00),
    ObservationCondition("intercept_plus_0.25", 0.25, 1.00),
    ObservationCondition("state_slope_0.80", 0.00, 0.80),
    ObservationCondition("combined_plus_0.25_slope_0.80", 0.25, 0.80),
)

function family_update(family::String, probability::Float64,
                       mean_log_odds::Float64, variance::Float64,
                       success::Float64, parameter::Float64)
    if family == "RL"
        updated_probability = min(max(
            probability + parameter * (success - probability),
            0.01,
        ), 0.99)
        return updated_probability, clipped_logit(updated_probability), variance
    elseif family == "Bayesian"
        predicted_variance = variance + parameter
        information = max(probability * (1.0 - probability), 1.0e-5)
        posterior_variance = 1.0 / (1.0 / predicted_variance + information)
        posterior_mean = min(max(
            mean_log_odds + posterior_variance * (success - probability),
            -4.60,
        ), 4.60)
        return sigmoid(posterior_mean), posterior_mean, posterior_variance
    end
    error("Unknown model family: $(family)")
end

function state_design(data::Dict{Symbol, Any}, columns::Vector{Symbol})
    n = length(data[columns[1]])
    matrix = zeros(Float64, n, length(columns))
    for column in eachindex(columns), row in 1:n
        matrix[row, column] = Float64(data[columns[column]][row])
    end
    matrix
end

function simulate_state_model_dataset(seed::Int, condition::StateModelCondition;
                                      n_participants::Int=STATE_N_PARTICIPANTS,
                                      n_items::Int=STATE_N_ITEMS)
    rng = DeterministicRNG(seed)
    n_pairs = n_participants * n_items
    n_rows = n_pairs * condition.n_stages
    ability = [0.55 * rand_normal!(rng) for _ in 1:n_participants]
    difficulty = [0.55 * rand_normal!(rng) for _ in 1:n_items]
    initial_learner = Vector{Float64}(undef, n_pairs)
    observed_initial = Vector{Float64}(undef, n_pairs)
    initial_scheduler = Vector{Float64}(undef, n_pairs)
    for pair in 1:n_pairs
        person = div(pair - 1, n_items) + 1
        card = mod(pair - 1, n_items) + 1
        initial_learner[pair] = 0.80 + ability[person] - difficulty[card] +
                                0.25 * rand_normal!(rng)
        observed_initial[pair] = initial_learner[pair] + 0.35 * rand_normal!(rng)
        initial_scheduler[pair] = 0.75 * initial_learner[pair] +
                                  0.45 * rand_normal!(rng)
    end

    learner_probability = sigmoid.(initial_learner)
    learner_mean = copy(initial_learner)
    learner_variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    scheduler_probability = sigmoid.(initial_scheduler)
    scheduler_mean = copy(initial_scheduler)
    scheduler_variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    learner_rl_probability = sigmoid.(observed_initial)
    learner_rl_mean = copy(observed_initial)
    learner_rl_variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    learner_bayes_probability = sigmoid.(observed_initial)
    learner_bayes_mean = copy(observed_initial)
    learner_bayes_variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    policy_rl_probability = sigmoid.(initial_scheduler)
    policy_rl_mean = copy(initial_scheduler)
    policy_rl_variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    policy_bayes_probability = sigmoid.(initial_scheduler)
    policy_bayes_mean = copy(initial_scheduler)
    policy_bayes_variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    previous_failure = zeros(Float64, n_pairs)

    participant = Vector{Int}(undef, n_rows)
    item = Vector{Int}(undef, n_rows)
    stage = Vector{Float64}(undef, n_rows)
    failure = Vector{Float64}(undef, n_rows)
    success = Vector{Float64}(undef, n_rows)
    previous_failure_all = Vector{Float64}(undef, n_rows)
    planned_log_gap = Vector{Float64}(undef, n_rows)
    actual_log_gap = Vector{Float64}(undef, n_rows)
    learner_true_log_odds = Vector{Float64}(undef, n_rows)
    scheduler_true_log_odds = Vector{Float64}(undef, n_rows)
    learner_rl_log_odds = Vector{Float64}(undef, n_rows)
    learner_bayes_log_odds = Vector{Float64}(undef, n_rows)
    policy_rl_log_odds = Vector{Float64}(undef, n_rows)
    policy_bayes_log_odds = Vector{Float64}(undef, n_rows)
    true_unsigned_pe = Vector{Float64}(undef, n_rows)
    rl_unsigned_pe = Vector{Float64}(undef, n_rows)
    bayes_unsigned_pe = Vector{Float64}(undef, n_rows)
    success_matrix = zeros(Float64, n_pairs, condition.n_stages)

    learner_parameter = condition.learner_family == "RL" ?
                        RL_LEARNING_RATE : BAYES_VOLATILITY
    scheduler_parameter = condition.scheduler_family == "RL" ?
                          RL_SCHEDULER_RATE : BAYES_SCHEDULER_VOLATILITY

    for current_stage in 0:condition.n_stages - 1
        offset = current_stage * n_pairs
        new_previous_failure = zeros(Float64, n_pairs)
        for pair in 1:n_pairs
            row = offset + pair
            learner_log_odds = clipped_logit(learner_probability[pair])
            scheduler_log_odds = clipped_logit(scheduler_probability[pair])
            learner_rl_log = clipped_logit(learner_rl_probability[pair])
            learner_bayes_log = clipped_logit(learner_bayes_probability[pair])
            policy_rl_log = clipped_logit(policy_rl_probability[pair])
            policy_bayes_log = clipped_logit(policy_bayes_probability[pair])
            planned = 0.75 + SCHEDULER_ADAPTIVITY * scheduler_log_odds -
                      0.55 * previous_failure[pair] + 0.20 * rand_normal!(rng)
            actual = planned + 0.50 * rand_normal!(rng)
            probability = sigmoid(
                OUTCOME_INTERCEPT + BETA_DELAY * actual +
                STATE_COEFFICIENT * learner_log_odds +
                PREVIOUS_FAILURE_COEFFICIENT * previous_failure[pair],
            )
            response = rand_bernoulli!(rng, probability)
            current_success = 1.0 - response
            success_matrix[pair, current_stage + 1] = current_success
            person = div(pair - 1, n_items) + 1
            card = mod(pair - 1, n_items) + 1
            participant[row] = person - 1
            item[row] = card - 1
            stage[row] = current_stage
            failure[row] = response
            success[row] = current_success
            previous_failure_all[row] = previous_failure[pair]
            planned_log_gap[row] = planned
            actual_log_gap[row] = actual
            learner_true_log_odds[row] = learner_log_odds
            scheduler_true_log_odds[row] = scheduler_log_odds
            learner_rl_log_odds[row] = learner_rl_log
            learner_bayes_log_odds[row] = learner_bayes_log
            policy_rl_log_odds[row] = policy_rl_log
            policy_bayes_log_odds[row] = policy_bayes_log
            true_unsigned_pe[row] = abs(current_success - learner_probability[pair])
            rl_unsigned_pe[row] = abs(current_success - learner_rl_probability[pair])
            bayes_unsigned_pe[row] = abs(current_success - learner_bayes_probability[pair])

            learner_probability[pair], learner_mean[pair], learner_variance[pair] =
                family_update(
                    condition.learner_family, learner_probability[pair],
                    learner_mean[pair], learner_variance[pair], current_success,
                    learner_parameter,
                )
            scheduler_probability[pair], scheduler_mean[pair], scheduler_variance[pair] =
                family_update(
                    condition.scheduler_family, scheduler_probability[pair],
                    scheduler_mean[pair], scheduler_variance[pair], current_success,
                    scheduler_parameter,
                )
            learner_rl_probability[pair], learner_rl_mean[pair], learner_rl_variance[pair] =
                family_update(
                    "RL", learner_rl_probability[pair], learner_rl_mean[pair],
                    learner_rl_variance[pair], current_success, RL_LEARNING_RATE,
                )
            learner_bayes_probability[pair], learner_bayes_mean[pair], learner_bayes_variance[pair] =
                family_update(
                    "Bayesian", learner_bayes_probability[pair], learner_bayes_mean[pair],
                    learner_bayes_variance[pair], current_success, BAYES_VOLATILITY,
                )
            policy_rl_probability[pair], policy_rl_mean[pair], policy_rl_variance[pair] =
                family_update(
                    "RL", policy_rl_probability[pair], policy_rl_mean[pair],
                    policy_rl_variance[pair], current_success, RL_SCHEDULER_RATE,
                )
            policy_bayes_probability[pair], policy_bayes_mean[pair], policy_bayes_variance[pair] =
                family_update(
                    "Bayesian", policy_bayes_probability[pair], policy_bayes_mean[pair],
                    policy_bayes_variance[pair], current_success,
                    BAYES_SCHEDULER_VOLATILITY,
                )
            new_previous_failure[pair] = response
        end
        previous_failure = new_previous_failure
    end

    actual_z, actual_scale = standardize(actual_log_gap)
    planned_z, _ = standardize(planned_log_gap)
    learner_true_z, _ = standardize(learner_true_log_odds)
    learner_rl_z, _ = standardize(learner_rl_log_odds)
    learner_bayes_z, _ = standardize(learner_bayes_log_odds)
    policy_rl_z, _ = standardize(policy_rl_log_odds)
    policy_bayes_z, _ = standardize(policy_bayes_log_odds)
    stage_z, _ = standardize(stage)
    inferred_plan_rl = 0.75 .+ SCHEDULER_ADAPTIVITY .* policy_rl_log_odds .-
                       0.55 .* previous_failure_all
    inferred_plan_bayes = 0.75 .+ SCHEDULER_ADAPTIVITY .* policy_bayes_log_odds .-
                          0.55 .* previous_failure_all
    inferred_plan_rl_z, _ = standardize(inferred_plan_rl)
    inferred_plan_bayes_z, _ = standardize(inferred_plan_bayes)

    data = Dict{Symbol, Any}(
        :participant => participant, :item => item, :stage => stage,
        :failure => failure, :success => success,
        :previous_failure => previous_failure_all,
        :planned_log_gap => planned_log_gap, :actual_log_gap => actual_log_gap,
        :learner_true_log_odds => learner_true_log_odds,
        :scheduler_true_log_odds => scheduler_true_log_odds,
        :learner_rl_log_odds => learner_rl_log_odds,
        :learner_bayes_log_odds => learner_bayes_log_odds,
        :policy_rl_log_odds => policy_rl_log_odds,
        :policy_bayes_log_odds => policy_bayes_log_odds,
        :true_unsigned_pe => true_unsigned_pe,
        :rl_unsigned_pe => rl_unsigned_pe,
        :bayes_unsigned_pe => bayes_unsigned_pe,
        :actual_z => actual_z, :planned_z => planned_z,
        :learner_true_z => learner_true_z, :learner_rl_z => learner_rl_z,
        :learner_bayes_z => learner_bayes_z,
        :policy_rl_z => policy_rl_z, :policy_bayes_z => policy_bayes_z,
        :stage_z => stage_z, :inferred_plan_rl => inferred_plan_rl,
        :inferred_plan_bayes => inferred_plan_bayes,
        :inferred_plan_rl_z => inferred_plan_rl_z,
        :inferred_plan_bayes_z => inferred_plan_bayes_z,
    )
    metadata = (
        truth=BETA_DELAY * actual_scale,
        initial_learner_log_odds=observed_initial,
        success_matrix=success_matrix,
        failure_rate=mean1(failure),
    )
    data, metadata
end

function fit_delay_estimators(data::Dict{Symbol, Any}, truth::Float64)
    designs = [
        ("naive", [:actual_z]),
        ("state_RL", [:actual_z, :learner_rl_z, :previous_failure, :stage_z]),
        ("state_Bayesian", [:actual_z, :learner_bayes_z, :previous_failure, :stage_z]),
        ("dual_RL_RL", [:actual_z, :learner_rl_z, :inferred_plan_rl_z, :previous_failure, :stage_z]),
        ("dual_RL_Bayesian", [:actual_z, :learner_rl_z, :inferred_plan_bayes_z, :previous_failure, :stage_z]),
        ("dual_Bayesian_RL", [:actual_z, :learner_bayes_z, :inferred_plan_rl_z, :previous_failure, :stage_z]),
        ("dual_Bayesian_Bayesian", [:actual_z, :learner_bayes_z, :inferred_plan_bayes_z, :previous_failure, :stage_z]),
        ("logged_plan_RL", [:actual_z, :planned_z, :learner_rl_z, :previous_failure, :stage_z]),
        ("logged_plan_Bayesian", [:actual_z, :planned_z, :learner_bayes_z, :previous_failure, :stage_z]),
        ("oracle", [:actual_z, :learner_true_z, :previous_failure, :stage_z]),
    ]
    rows = Any[]
    for (estimator, columns) in designs
        fit = logistic_cluster_fit(
            state_design(data, columns), data[:failure], data[:participant],
        )
        estimate = fit.estimate
        push!(rows, (
            estimator=estimator, estimate=estimate,
            standard_error=fit.standard_error, ci_low=fit.ci_low,
            ci_high=fit.ci_high, converged=fit.converged, truth=truth,
            bias=estimate - truth, squared_error=(estimate - truth)^2,
            covered=fit.ci_low <= truth <= fit.ci_high,
            sign_error=estimate < 0.0,
        ))
    end
    rows
end

function reconstruct_state_history(initial_state::Vector{Float64},
                                   success_matrix::Matrix{Float64},
                                   family::String, parameter::Float64)
    n_pairs, n_stages = size(success_matrix)
    probability = sigmoid.(initial_state)
    mean_state = copy(initial_state)
    variance = fill(INITIAL_BAYES_VARIANCE, n_pairs)
    history = zeros(Float64, n_pairs, n_stages)
    for current_stage in 1:n_stages
        for pair in 1:n_pairs
            history[pair, current_stage] = clipped_logit(probability[pair])
            probability[pair], mean_state[pair], variance[pair] = family_update(
                family, probability[pair], mean_state[pair], variance[pair],
                success_matrix[pair, current_stage], parameter,
            )
        end
    end
    history
end

function score_state_matrix(data::Dict{Symbol, Any}, state_matrix::Matrix{Float64},
                            observation::ObservationCondition)
    outcome = data[:failure]
    negative_log_likelihood = 0.0
    n_pairs, n_stages = size(state_matrix)
    for current_stage in 1:n_stages, pair in 1:n_pairs
        row = (current_stage - 1) * n_pairs + pair
        predictor = OUTCOME_INTERCEPT + observation.intercept_shift +
                    BETA_DELAY * data[:actual_log_gap][row] +
                    STATE_COEFFICIENT * observation.state_slope_multiplier *
                    state_matrix[pair, current_stage] +
                    PREVIOUS_FAILURE_COEFFICIENT * data[:previous_failure][row]
        probability = min(max(sigmoid(predictor), 1.0e-8), 1.0 - 1.0e-8)
        negative_log_likelihood -= outcome[row] * log(probability) +
                                   (1.0 - outcome[row]) * log(1.0 - probability)
    end
    negative_log_likelihood
end

function score_candidate_models_detailed(data::Dict{Symbol, Any}, initial_state,
                                         success_matrix,
                                         observation::ObservationCondition)
    scores = Any[]
    for (family, grid) in (("RL", RL_PARAMETER_GRID),
                           ("Bayesian", BAYES_PARAMETER_GRID))
        for parameter in grid
            state_matrix = reconstruct_state_history(
                initial_state, success_matrix, family, parameter,
            )
            push!(scores, (family=family, parameter=parameter,
                           nll=score_state_matrix(data, state_matrix, observation)))
        end
    end
    rl_scores = [score for score in scores if score.family == "RL"]
    bayes_scores = [score for score in scores if score.family == "Bayesian"]
    rl_best = rl_scores[argmin([score.nll for score in rl_scores])]
    bayes_best = bayes_scores[argmin([score.nll for score in bayes_scores])]
    selected = rl_best.nll <= bayes_best.nll ? rl_best : bayes_best
    candidate_rows = Any[]
    later = [index for index in eachindex(data[:stage]) if data[:stage][index] > 0.0]
    true_state = subset(data[:learner_true_log_odds], later)
    true_pe = subset(data[:true_unsigned_pe], later)
    for best in (rl_best, bayes_best)
        state_matrix = reconstruct_state_history(
            initial_state, success_matrix, best.family, best.parameter,
        )
        state_vector = vec(state_matrix)
        candidate_state = subset(state_vector, later)
        candidate_pe_all = abs.(vec(success_matrix) .- sigmoid.(state_vector))
        candidate_pe = subset(candidate_pe_all, later)
        push!(candidate_rows, (
            candidate_family=best.family,
            best_parameter=best.parameter,
            nll=best.nll,
            selected_family=selected.family,
            state_rmse=rmse(state_vector, data[:learner_true_log_odds], later),
            state_correlation=correlation(candidate_state, true_state),
            unsigned_pe_correlation=correlation(candidate_pe, true_pe),
        ))
    end
    (
        selected_family=selected.family,
        selected_parameter=selected.parameter,
        rl_parameter_within_family=rl_best.parameter,
        bayesian_parameter_within_family=bayes_best.parameter,
        rl_nll=rl_best.nll,
        bayesian_nll=bayes_best.nll,
        candidates=candidate_rows,
    )
end

function score_candidate_models(data::Dict{Symbol, Any}, initial_state,
                                success_matrix)
    result = score_candidate_models_detailed(
        data, initial_state, success_matrix,
        ObservationCondition("correct", 0.0, 1.0),
    )
    (
        selected_family=result.selected_family,
        selected_parameter=result.selected_parameter,
        rl_parameter_within_family=result.rl_parameter_within_family,
        bayesian_parameter_within_family=result.bayesian_parameter_within_family,
        rl_nll=result.rl_nll,
        bayesian_nll=result.bayesian_nll,
    )
end

function rmse(left, right, indices)
    total = 0.0
    for index in indices
        delta = left[index] - right[index]
        total += delta * delta
    end
    sqrt(total / length(indices))
end

function subset(values, indices)
    [Float64(values[index]) for index in indices]
end

function recovery_diagnostics(data::Dict{Symbol, Any})
    later = [index for index in eachindex(data[:stage]) if data[:stage][index] > 0.0]
    all_indices = collect(eachindex(data[:stage]))
    true_state = subset(data[:learner_true_log_odds], later)
    true_pe = subset(data[:true_unsigned_pe], later)
    (
        rl_state_rmse=rmse(data[:learner_rl_log_odds], data[:learner_true_log_odds], later),
        bayesian_state_rmse=rmse(data[:learner_bayes_log_odds], data[:learner_true_log_odds], later),
        rl_state_correlation=correlation(subset(data[:learner_rl_log_odds], later), true_state),
        bayesian_state_correlation=correlation(subset(data[:learner_bayes_log_odds], later), true_state),
        rl_unsigned_pe_correlation=correlation(subset(data[:rl_unsigned_pe], later), true_pe),
        bayesian_unsigned_pe_correlation=correlation(subset(data[:bayes_unsigned_pe], later), true_pe),
        rl_policy_rmse=rmse(data[:inferred_plan_rl], data[:planned_log_gap], all_indices),
        bayesian_policy_rmse=rmse(data[:inferred_plan_bayes], data[:planned_log_gap], all_indices),
    )
end

const DELAY_COLUMNS = [
    "replication", "learner_family", "scheduler_family", "n_stages",
    "n_participants", "n_items", "failure_rate", "estimator", "estimate",
    "standard_error", "ci_low", "ci_high", "converged", "truth", "bias",
    "squared_error", "covered", "sign_error",
]

const DIAGNOSTIC_COLUMNS = [
    "replication", "learner_family", "scheduler_family", "n_stages",
    "failure_rate", "selected_family", "selected_parameter",
    "rl_parameter_within_family", "bayesian_parameter_within_family",
    "rl_nll", "bayesian_nll", "rl_state_rmse", "bayesian_state_rmse",
    "rl_state_correlation", "bayesian_state_correlation",
    "rl_unsigned_pe_correlation", "bayesian_unsigned_pe_correlation",
    "rl_policy_rmse", "bayesian_policy_rmse",
]

const OBSERVATION_MISSPEC_COLUMNS = [
    "replication", "learner_family", "scheduler_family", "n_stages",
    "observation_condition", "intercept_shift", "state_slope_multiplier",
    "failure_rate", "candidate_family", "best_parameter", "true_parameter",
    "parameter_error", "nll", "selected_family", "family_recovered",
    "state_rmse", "state_correlation", "unsigned_pe_correlation",
]

function run_observation_misspecification(output_dir::String;
                                          n_replications::Int=100,
                                          seed_offset::Int=130_000_000)
    path = joinpath(output_dir, "state_model_observation_misspec_replications.csv")
    open(path, "w") do io
        write_csv_header(io, OBSERVATION_MISSPEC_COLUMNS)
        for learner_family in ("RL", "Bayesian")
            true_parameter = learner_family == "RL" ? RL_LEARNING_RATE : BAYES_VOLATILITY
            for scheduler_family in ("RL", "Bayesian")
                condition = StateModelCondition(learner_family, scheduler_family, 8)
                println("[Julia observation] learner=$(learner_family), scheduler=$(scheduler_family)")
                for replication in 0:n_replications - 1
                    seed = STATE_SEED + seed_offset +
                           (learner_family == "RL" ? 0 : 10_000_000) +
                           (scheduler_family == "RL" ? 0 : 1_000_000) +
                           condition.n_stages * 10_000 +
                           replication
                    data, metadata = simulate_state_model_dataset(seed, condition)
                    for observation in OBSERVATION_CONDITIONS
                        result = score_candidate_models_detailed(
                            data, metadata.initial_learner_log_odds,
                            metadata.success_matrix, observation,
                        )
                        for candidate in result.candidates
                            println(io, csv_row(
                                replication, learner_family, scheduler_family, 8,
                                observation.label, observation.intercept_shift,
                                observation.state_slope_multiplier,
                                metadata.failure_rate, candidate.candidate_family,
                                candidate.best_parameter, true_parameter,
                                candidate.best_parameter - true_parameter,
                                candidate.nll, candidate.selected_family,
                                candidate.selected_family == learner_family,
                                candidate.state_rmse, candidate.state_correlation,
                                candidate.unsigned_pe_correlation,
                            ))
                        end
                    end
                end
            end
        end
    end
    path
end

function main()
    project_dir = normpath(joinpath(@__DIR__, "..", ".."))
    output_dir = ensure_directory(get(ENV, "BRM_OUTPUT_DIR", joinpath(project_dir, "outputs")))
    n_replications = parse(Int, get(ENV, "BRM_STATE_REPS", "100"))
    observation_replications = parse(
        Int, get(ENV, "BRM_STATE_OBS_REPS", string(n_replications)),
    )
    observation_seed_offset = parse(
        Int, get(ENV, "BRM_STATE_OBS_SEED_OFFSET", "130000000"),
    )
    observation_only = lowercase(get(ENV, "BRM_STATE_OBS_ONLY", "false")) in
                       ("true", "1", "yes")
    delay_path = joinpath(output_dir, "state_model_delay_recovery_replications.csv")
    diagnostic_path = joinpath(output_dir, "state_model_diagnostics_replications.csv")
    if !observation_only
        open(delay_path, "w") do delay_io
            open(diagnostic_path, "w") do diagnostic_io
                write_csv_header(delay_io, DELAY_COLUMNS)
                write_csv_header(diagnostic_io, DIAGNOSTIC_COLUMNS)
                for learner_family in ("RL", "Bayesian")
                    for scheduler_family in ("RL", "Bayesian")
                        for n_stages in (2, 8)
                            condition = StateModelCondition(
                                learner_family, scheduler_family, n_stages,
                            )
                            println("[Julia state] learner=$(learner_family), scheduler=$(scheduler_family), stages=$(n_stages)")
                            for replication in 0:n_replications - 1
                                seed = STATE_SEED + 80_000_000 +
                                       (learner_family == "RL" ? 0 : 10_000_000) +
                                       (scheduler_family == "RL" ? 0 : 1_000_000) +
                                       n_stages * 10_000 + replication
                                data, metadata = simulate_state_model_dataset(seed, condition)
                                for result in fit_delay_estimators(data, metadata.truth)
                                    println(delay_io, csv_row(
                                        replication, learner_family, scheduler_family,
                                        n_stages, STATE_N_PARTICIPANTS, STATE_N_ITEMS,
                                        metadata.failure_rate, result.estimator,
                                        result.estimate, result.standard_error,
                                        result.ci_low, result.ci_high, result.converged,
                                        result.truth, result.bias, result.squared_error,
                                        result.covered, result.sign_error,
                                    ))
                                end
                                model = score_candidate_models(
                                    data, metadata.initial_learner_log_odds,
                                    metadata.success_matrix,
                                )
                                diagnostic = recovery_diagnostics(data)
                                println(diagnostic_io, csv_row(
                                    replication, learner_family, scheduler_family,
                                    n_stages, metadata.failure_rate,
                                    model.selected_family, model.selected_parameter,
                                    model.rl_parameter_within_family,
                                    model.bayesian_parameter_within_family,
                                    model.rl_nll, model.bayesian_nll,
                                    diagnostic.rl_state_rmse,
                                    diagnostic.bayesian_state_rmse,
                                    diagnostic.rl_state_correlation,
                                    diagnostic.bayesian_state_correlation,
                                    diagnostic.rl_unsigned_pe_correlation,
                                    diagnostic.bayesian_unsigned_pe_correlation,
                                    diagnostic.rl_policy_rmse,
                                    diagnostic.bayesian_policy_rmse,
                                ))
                            end
                        end
                    end
                end
            end
        end
    end
    run_observation_misspecification(
        output_dir, n_replications=observation_replications,
        seed_offset=observation_seed_offset,
    )
    println("Julia state-model simulations completed in $(output_dir)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
