using Distributions

Base.@kwdef mutable struct DyadicTruth

    distribution::String
    f::Vector{Float64} = fill(NaN, (n_evals))
    Var_f_given_A::Vector{Float64} = fill(NaN, (n_evals))
end



function DyadicTruth(params::Dict)

    n_evals = params["n_evals"]
    return DyadicTruth(
        params["data_distribution"],
        fill(NaN, (n_evals)),
        fill(NaN, (n_evals)),
    )
end



function DyadicTruth(params::Dict,
                     evals::DyadicEvaluationPoints)

    truth = DyadicTruth(params)

    if params["data_distribution"] == "triplet_gaussian_sum_product"
        generate_truth_triplet_gaussian_sum_product(truth, evals, params)
    else
        error("unknown distribution")
    end

    return truth
end



function generate_truth_triplet_gaussian_sum_product(
    truth::DyadicTruth, evals::DyadicEvaluationPoints, params::Dict)

    p = params["data_p"]
    generate_f_triplet_gaussian_sum_product(truth, evals, p)
    generate_Var_f_given_A_triplet_gaussian_sum_product(truth, evals, p)
end



function generate_f_triplet_gaussian_sum_product(
    truth::DyadicTruth, evals::DyadicEvaluationPoints, p::Vector{Float64})

    truth.f = (p[1]^2 + p[3]^2) * phi.(evals.w .- 1)
    truth.f += p[2]*(2 - p[2]) * phi.(evals.w)
    truth.f += 2*p[1]*p[3] * phi.(evals.w .+ 1)
end



function generate_Var_f_given_A_triplet_gaussian_sum_product(
    truth::DyadicTruth, evals::DyadicEvaluationPoints, p::Vector{Float64})

    truth.Var_f_given_A = (p[1] + p[3])^2 * p[2] * phi.(evals.w).^2
    truth.Var_f_given_A += p[3] * (p[3] * phi.(evals.w .- 1) + p[1] * phi.(evals.w .+ 1)).^2
    truth.Var_f_given_A += p[1] * (p[1] * phi.(evals.w .- 1) + p[3] * phi.(evals.w .+ 1)).^2
    truth.Var_f_given_A -= ((p[1]^2 + p[3]^2) * phi.(evals.w .- 1) +
        (p[1] + p[3]) * p[2] * phi.(evals.w) + 2 * p[1] * p[3] * phi.(evals.w .+ 1)).^2

end
