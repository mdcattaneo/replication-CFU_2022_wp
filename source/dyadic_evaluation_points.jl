using Printf

Base.@kwdef mutable struct DyadicEvaluationPoints

    n_evals::Int
    w::Vector{Float64} = fill(NaN, (n_evals))
    w_min::Float64 = NaN
    w_max::Float64 = NaN

end



function DyadicEvaluationPoints(params::Dict)

    @assert w_min < w_max

    evals = DyadicEvaluationPoints(n_evals = params["n_evals"])

    if params["evals_type"] == "deterministic"
        generate_evals_deterministic(evals, params["evals_w_min"], params["evals_w_max"])
    else
        error("unknown evals_type")
    end

    return evals

end



function generate_evals_deterministic(
    evals::DyadicEvaluationPoints,
    w_min::Float64, w_max::Float64)

    evals.w_min = w_min
    evals.w_max = w_max

    w_range = LinRange(w_min, w_max, evals.n_evals)

    for i in 1:evals.n_evals
        evals.w[i] = w_range[i]
    end

    return nothing
end
