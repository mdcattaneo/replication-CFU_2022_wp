Base.@kwdef mutable struct DyadicDiagnosis

    MSEhat::Vector{Float64} = Float64[]
    IMSEhat::Float64 = NaN
    MEhat::Float64 = NaN

    ucb_coverage::Bool = false
    ucb_max_width::Float64 = NaN
    ucb_mean_width::Float64 = NaN

    bcb_coverage::Bool = false
    bcb_max_width::Float64 = NaN
    bcb_mean_width::Float64 = NaN

    pcb_coverage::Bool = false
    pcb_max_width::Float64 = NaN
    pcb_mean_width::Float64 = NaN

end



function DyadicDiagnosis(
    truth::DyadicTruth,
    estimator::DyadicKernelDensityEstimator,
    evals::DyadicEvaluationPoints)

    ucb_diagnosis = diagnose_band(truth.f, estimator.ucb[1,:], estimator.ucb[2,:])
    bcb_diagnosis = diagnose_band(truth.f, estimator.bcb[1,:], estimator.bcb[2,:])
    pcb_diagnosis = diagnose_band(truth.f, estimator.pcb[1,:], estimator.pcb[2,:])

    diagnosis = DyadicDiagnosis(

        Float64[],
        NaN,
        NaN,

        ucb_diagnosis["coverage"],
        ucb_diagnosis["max_width"],
        ucb_diagnosis["mean_width"],

        bcb_diagnosis["coverage"],
        bcb_diagnosis["max_width"],
        bcb_diagnosis["mean_width"],

        pcb_diagnosis["coverage"],
        pcb_diagnosis["max_width"],
        pcb_diagnosis["mean_width"],
    )

    get_MSEhat(diagnosis, truth, estimator)
    get_IMSEhat(diagnosis, evals)
    get_MEhat(diagnosis, truth, estimator),

    return diagnosis

end



function diagnose_band(true_values::Vector{Float64},
                       lower_bounds::Vector{Float64},
                       upper_bounds::Vector{Float64})

    coverage = get_band_coverage(true_values, lower_bounds, upper_bounds)
    max_width = get_max_band_width(lower_bounds, upper_bounds)
    mean_width = get_mean_band_width(lower_bounds, upper_bounds)

    dict = Dict("coverage" => coverage,
                "max_width" => max_width,
                "mean_width" => mean_width)

    return dict

end



function get_max_band_width(lower_bounds::Vector{Float64},
                            upper_bounds::Vector{Float64})

    return maximum(upper_bounds .- lower_bounds)
end



function get_mean_band_width(lower_bounds::Vector{Float64},
                             upper_bounds::Vector{Float64})

    return sum(upper_bounds .- lower_bounds) / length(upper_bounds)
end



function get_MSEhat(diagnosis::DyadicDiagnosis, truth::DyadicTruth,
                    estimator::DyadicKernelDensityEstimator)

    diagnosis.MSEhat = (estimator.fhat .- truth.f).^2
end



function get_IMSEhat(diagnosis::DyadicDiagnosis,
                     evals::DyadicEvaluationPoints)

    diagnosis.IMSEhat = mean(diagnosis.MSEhat)
end



function get_MEhat(diagnosis::DyadicDiagnosis, truth::DyadicTruth,
                    estimator::DyadicKernelDensityEstimator)

    diagnosis.MEhat = maximum(abs.(estimator.fhat .- truth.f))
end



function get_band_coverage(truth_values::Vector{Float64},
                           lower_bounds::Vector{Float64},
                           upper_bounds::Vector{Float64})

    return all(lower_bounds .<= truth_values .<= upper_bounds)
end
