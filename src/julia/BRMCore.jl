module BRMCore

# Pure-Julia computational core.  Only Base is required so that the simulation
# layer remains portable across laboratory workstations, clusters, and CI jobs.

export DeterministicRNG, rand_uniform!, rand_normal!, rand_bernoulli!,
       rand_sign!, sigmoid, clipped_logit, standardize, mean1, sd1,
       correlation, spearman_correlation, logistic_cluster_fit,
       ols_cluster_fit, wcls_cluster_fit, student_t_critical_975,
       csv_row, write_csv_header, ensure_directory, isfinite_number

mutable struct DeterministicRNG
    state::UInt32
    has_spare::Bool
    spare::Float64
end

function DeterministicRNG(seed::Integer)
    state = UInt32(mod(seed, 4_294_967_295))
    state == UInt32(0) && (state = UInt32(0x6d2b79f5))
    DeterministicRNG(state, false, 0.0)
end

function next_u32!(rng::DeterministicRNG)
    x = rng.state
    x = xor(x, x << 13)
    x = xor(x, x >> 17)
    x = xor(x, x << 5)
    rng.state = x
    x
end

rand_uniform!(rng::DeterministicRNG) =
    (Float64(next_u32!(rng)) + 0.5) / 4_294_967_296.0

function rand_normal!(rng::DeterministicRNG)
    if rng.has_spare
        rng.has_spare = false
        return rng.spare
    end
    u1 = max(rand_uniform!(rng), 1.0e-15)
    u2 = rand_uniform!(rng)
    radius = sqrt(-2.0 * log(u1))
    angle = 2.0 * pi * u2
    rng.spare = radius * sin(angle)
    rng.has_spare = true
    radius * cos(angle)
end

rand_bernoulli!(rng::DeterministicRNG, probability::Real) =
    rand_uniform!(rng) < probability ? 1.0 : 0.0

rand_sign!(rng::DeterministicRNG) = rand_uniform!(rng) < 0.5 ? -1.0 : 1.0

clip(x::Real, lower::Real, upper::Real) = min(max(Float64(x), Float64(lower)), Float64(upper))

function sigmoid(x::Real)
    value = Float64(x)
    if value >= 0.0
        z = exp(-value)
        return 1.0 / (1.0 + z)
    end
    z = exp(value)
    z / (1.0 + z)
end

clipped_logit(probability::Real) = begin
    p = clip(probability, 0.01, 0.99)
    log(p / (1.0 - p))
end

function mean1(values)
    n = length(values)
    n == 0 && return NaN
    total = 0.0
    for value in values
        total += Float64(value)
    end
    total / n
end

function sd1(values; corrected::Bool=true)
    n = length(values)
    n <= (corrected ? 1 : 0) && return NaN
    center = mean1(values)
    total = 0.0
    for value in values
        delta = Float64(value) - center
        total += delta * delta
    end
    sqrt(total / (corrected ? n - 1 : n))
end

function standardize(values)
    center = mean1(values)
    scale = sd1(values)
    if !isfinite(scale) || scale <= 0.0
        return zeros(Float64, length(values)), 1.0
    end
    result = Vector{Float64}(undef, length(values))
    for index in eachindex(values)
        result[index] = (Float64(values[index]) - center) / scale
    end
    result, scale
end

function correlation(left, right)
    n = min(length(left), length(right))
    n < 3 && return NaN
    left_mean = mean1(left[1:n])
    right_mean = mean1(right[1:n])
    cross = 0.0
    left_ss = 0.0
    right_ss = 0.0
    for index in 1:n
        a = Float64(left[index]) - left_mean
        b = Float64(right[index]) - right_mean
        cross += a * b
        left_ss += a * a
        right_ss += b * b
    end
    denominator = sqrt(left_ss * right_ss)
    denominator <= 0.0 ? NaN : cross / denominator
end

function average_ranks(values)
    n = length(values)
    order = sortperm(collect(1:n), by=index -> Float64(values[index]))
    ranks = zeros(Float64, n)
    position = 1
    while position <= n
        last = position
        reference = Float64(values[order[position]])
        while last < n && Float64(values[order[last + 1]]) == reference
            last += 1
        end
        average_rank = (position + last) / 2.0
        for cursor in position:last
            ranks[order[cursor]] = average_rank
        end
        position = last + 1
    end
    ranks
end

spearman_correlation(left, right) = correlation(average_ranks(left), average_ranks(right))

function solve_linear(matrix::Matrix{Float64}, vector::Vector{Float64})
    n = size(matrix, 1)
    augmented = zeros(Float64, n, n + 1)
    for row in 1:n
        for column in 1:n
            augmented[row, column] = matrix[row, column]
        end
        augmented[row, n + 1] = vector[row]
    end
    for pivot in 1:n
        best_row = pivot
        best_value = abs(augmented[pivot, pivot])
        for row in pivot + 1:n
            candidate = abs(augmented[row, pivot])
            if candidate > best_value
                best_row = row
                best_value = candidate
            end
        end
        if best_row != pivot
            for column in pivot:n + 1
                temporary = augmented[pivot, column]
                augmented[pivot, column] = augmented[best_row, column]
                augmented[best_row, column] = temporary
            end
        end
        if abs(augmented[pivot, pivot]) < 1.0e-12
            augmented[pivot, pivot] += 1.0e-8
        end
        pivot_value = augmented[pivot, pivot]
        for column in pivot:n + 1
            augmented[pivot, column] /= pivot_value
        end
        for row in 1:n
            row == pivot && continue
            multiplier = augmented[row, pivot]
            multiplier == 0.0 && continue
            for column in pivot:n + 1
                augmented[row, column] -= multiplier * augmented[pivot, column]
            end
        end
    end
    augmented[:, n + 1]
end

function inverse_matrix(matrix::Matrix{Float64})
    n = size(matrix, 1)
    result = zeros(Float64, n, n)
    for column in 1:n
        target = zeros(Float64, n)
        target[column] = 1.0
        solution = solve_linear(matrix, target)
        for row in 1:n
            result[row, column] = solution[row]
        end
    end
    result
end

function matrix_product(left::Matrix{Float64}, right::Matrix{Float64})
    rows = size(left, 1)
    inner = size(left, 2)
    columns = size(right, 2)
    result = zeros(Float64, rows, columns)
    for row in 1:rows
        for column in 1:columns
            total = 0.0
            for cursor in 1:inner
                total += left[row, cursor] * right[cursor, column]
            end
            result[row, column] = total
        end
    end
    result
end

softplus(x::Real) = x > 30.0 ? Float64(x) : log1p(exp(Float64(x)))

function logistic_objective(design, outcome, beta, ridge)
    n = length(outcome)
    p0 = size(design, 2)
    value = 0.0
    for row in 1:n
        eta = beta[1]
        for column in 1:p0
            eta += design[row, column] * beta[column + 1]
        end
        eta = clip(eta, -30.0, 30.0)
        value += outcome[row] * eta - softplus(eta)
    end
    for column in 2:length(beta)
        value -= 0.5 * ridge * beta[column] * beta[column]
    end
    value
end

function design_value(design, row, column)
    column == 1 ? 1.0 : design[row, column - 1]
end

function logistic_cluster_fit(
    design::Matrix{Float64},
    outcome::Vector{Float64},
    clusters::Vector{Int};
    focal_index::Int=1,
    ridge::Float64=1.0e-7,
    max_iter::Int=80,
)
    n = length(outcome)
    p = size(design, 2) + 1
    beta = zeros(Float64, p)
    response_mean = clip(mean1(outcome), 1.0e-5, 1.0 - 1.0e-5)
    beta[1] = log(response_mean / (1.0 - response_mean))
    converged = false
    current_objective = logistic_objective(design, outcome, beta, ridge)

    for iteration in 1:max_iter
        information = zeros(Float64, p, p)
        score = zeros(Float64, p)
        for row in 1:n
            eta = beta[1]
            for column in 2:p
                eta += design[row, column - 1] * beta[column]
            end
            probability = sigmoid(clip(eta, -30.0, 30.0))
            weight = max(probability * (1.0 - probability), 1.0e-7)
            residual = outcome[row] - probability
            for column in 1:p
                x_column = design_value(design, row, column)
                score[column] += x_column * residual
                for other in 1:p
                    information[column, other] +=
                        weight * x_column * design_value(design, row, other)
                end
            end
        end
        for column in 2:p
            information[column, column] += ridge
            score[column] -= ridge * beta[column]
        end
        step = solve_linear(information, score)
        fraction = 1.0
        candidate = beta .+ step
        candidate_objective = logistic_objective(design, outcome, candidate, ridge)
        while candidate_objective < current_objective && fraction > 1.0e-4
            fraction *= 0.5
            candidate = beta .+ fraction .* step
            candidate_objective = logistic_objective(design, outcome, candidate, ridge)
        end
        beta = candidate
        current_objective = candidate_objective
        largest_step = 0.0
        for value in step
            largest_step = max(largest_step, abs(fraction * value))
        end
        if largest_step < 1.0e-7
            converged = true
            break
        end
    end

    information = zeros(Float64, p, p)
    probability = zeros(Float64, n)
    for row in 1:n
        eta = beta[1]
        for column in 2:p
            eta += design[row, column - 1] * beta[column]
        end
        probability[row] = sigmoid(clip(eta, -30.0, 30.0))
        weight = max(probability[row] * (1.0 - probability[row]), 1.0e-7)
        for column in 1:p
            x_column = design_value(design, row, column)
            for other in 1:p
                information[column, other] +=
                    weight * x_column * design_value(design, row, other)
            end
        end
    end
    for column in 2:p
        information[column, column] += ridge
    end
    bread = inverse_matrix(information)
    unique_clusters = sort(unique(clusters))
    score_by_cluster = Dict{Int, Vector{Float64}}()
    for cluster in unique_clusters
        score_by_cluster[cluster] = zeros(Float64, p)
    end
    for row in 1:n
        residual = outcome[row] - probability[row]
        target = score_by_cluster[clusters[row]]
        for column in 1:p
            target[column] += design_value(design, row, column) * residual
        end
    end
    meat = zeros(Float64, p, p)
    for cluster in unique_clusters
        cluster_score = score_by_cluster[cluster]
        for row in 1:p, column in 1:p
            meat[row, column] += cluster_score[row] * cluster_score[column]
        end
    end
    covariance = matrix_product(matrix_product(bread, meat), bread)
    group_count = length(unique_clusters)
    if group_count > 1 && n > p
        correction = (group_count / (group_count - 1.0)) * ((n - 1.0) / (n - p))
        covariance .*= correction
    end
    coefficient_index = focal_index + 1
    standard_error = sqrt(max(covariance[coefficient_index, coefficient_index], 0.0))
    estimate = beta[coefficient_index]
    (
        estimate=estimate,
        standard_error=standard_error,
        ci_low=estimate - 1.96 * standard_error,
        ci_high=estimate + 1.96 * standard_error,
        converged=converged,
    )
end

function ols_cluster_fit(
    design::Matrix{Float64},
    outcome::Vector{Float64},
    clusters::Vector{Int};
    focal_index::Int=1,
)
    n = length(outcome)
    p = size(design, 2) + 1
    cross_product = zeros(Float64, p, p)
    cross_outcome = zeros(Float64, p)
    for row in 1:n
        for column in 1:p
            x_column = design_value(design, row, column)
            cross_outcome[column] += x_column * outcome[row]
            for other in 1:p
                cross_product[column, other] +=
                    x_column * design_value(design, row, other)
            end
        end
    end
    for column in 1:p
        cross_product[column, column] += 1.0e-10
    end
    beta = solve_linear(cross_product, cross_outcome)
    residual = zeros(Float64, n)
    for row in 1:n
        fitted = beta[1]
        for column in 2:p
            fitted += design[row, column - 1] * beta[column]
        end
        residual[row] = outcome[row] - fitted
    end
    bread = inverse_matrix(cross_product)
    unique_clusters = sort(unique(clusters))
    score_by_cluster = Dict{Int, Vector{Float64}}()
    for cluster in unique_clusters
        score_by_cluster[cluster] = zeros(Float64, p)
    end
    for row in 1:n
        target = score_by_cluster[clusters[row]]
        for column in 1:p
            target[column] += design_value(design, row, column) * residual[row]
        end
    end
    meat = zeros(Float64, p, p)
    for cluster in unique_clusters
        cluster_score = score_by_cluster[cluster]
        for row in 1:p, column in 1:p
            meat[row, column] += cluster_score[row] * cluster_score[column]
        end
    end
    covariance = matrix_product(matrix_product(bread, meat), bread)
    group_count = length(unique_clusters)
    if group_count > 1 && n > p
        correction = (group_count / (group_count - 1.0)) * ((n - 1.0) / (n - p))
        covariance .*= correction
    end
    coefficient_index = focal_index + 1
    standard_error = sqrt(max(covariance[coefficient_index, coefficient_index], 0.0))
    estimate = beta[coefficient_index]
    (
        estimate=estimate,
        standard_error=standard_error,
        ci_low=estimate - 1.96 * standard_error,
        ci_high=estimate + 1.96 * standard_error,
        converged=true,
    )
end

"""
    wcls_cluster_fit(treatment, nuisance, outcome, clusters;
                     randomization_probability=0.5)

Fit the additive proximal-effect working model used in weighted and centered
least squares (WCLS).  `treatment` is coded 0/1 and the first design column is
`treatment - randomization_probability`; `nuisance` contains the pre-treatment
covariates used for precision.  In the micro-randomized experiment implemented
in this project the randomization probability is fixed at 0.5 and the reference
and actual randomization probabilities coincide.  Consequently all WCLS
weights equal one, and the treatment coefficient and cluster-robust standard
error must agree with the corresponding additive linear-probability model up
to floating-point arithmetic.  Keeping this as a separate estimator makes that
equivalence an explicit, testable implementation check.
"""
function wcls_cluster_fit(
    treatment::Vector{Float64},
    nuisance::Matrix{Float64},
    outcome::Vector{Float64},
    clusters::Vector{Int};
    randomization_probability::Float64=0.5,
)
    0.0 < randomization_probability < 1.0 ||
        error("randomization_probability must lie strictly between zero and one")
    n = length(outcome)
    length(treatment) == n || error("treatment and outcome lengths differ")
    length(clusters) == n || error("cluster and outcome lengths differ")
    size(nuisance, 1) == n || error("nuisance design and outcome lengths differ")

    design = zeros(Float64, n, size(nuisance, 2) + 1)
    for row in 1:n
        design[row, 1] = treatment[row] - randomization_probability
        for column in 1:size(nuisance, 2)
            design[row, column + 1] = nuisance[row, column]
        end
    end
    ols_cluster_fit(design, outcome, clusters, focal_index=1)
end

"""
    student_t_critical_975(degrees_of_freedom)

Two-sided 95% Student-*t* critical values used by the participant-cluster
Monte Carlo grid.  The supported degrees of freedom are fixed by the declared
sample-size grid, so tabulating the authoritative quantiles avoids adding a
non-Base Julia dependency to the reproducibility core.
"""
function student_t_critical_975(degrees_of_freedom::Int)
    critical = Dict(
        11 => 2.200985160091638,
        29 => 2.045229642132703,
        59 => 2.000995377048210,
        99 => 1.984216951508683,
    )
    haskey(critical, degrees_of_freedom) ||
        error("No predeclared t critical value for df=$(degrees_of_freedom)")
    critical[degrees_of_freedom]
end

function csv_escape(value)
    if value === missing || value === nothing
        return ""
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa AbstractFloat
        return isfinite(value) ? string(value) : ""
    end
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    text
end

csv_row(values...) = join((csv_escape(value) for value in values), ",")

write_csv_header(io, columns) = println(io, join(columns, ","))

ensure_directory(path::AbstractString) = isdir(path) ? path : (mkpath(path); path)

isfinite_number(value) = value isa Real && isfinite(Float64(value))

end
