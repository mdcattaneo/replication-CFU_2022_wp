using Printf
using Dates


struct DyadicExperiment

    params::Dict
    data::DyadicData
    evals::DyadicEvaluationPoints
    estimator::DyadicKernelDensityEstimator
    truth::DyadicTruth
    diagnosis::DyadicDiagnosis
end


Base.@kwdef mutable struct DyadicLightweightExperiment

    # params
    params::Dict

    # estimator
    fhat::Vector{Float64}
    f::Vector{Float64}
    psd_error::Float64
    min_eig_Sigmahat::Float64
    min_eig_Sigmahatplus::Float64
    bandwidth_ROT::Float64

    # diagnosis
    IMSEhat::Float64
    MEhat::Float64
    ucb_coverage::Bool
    ucb_max_width::Float64
    ucb_mean_width::Float64
    bcb_coverage::Bool
    bcb_max_width::Float64
    bcb_mean_width::Float64
    pcb_coverage::Bool
    pcb_max_width::Float64
    pcb_mean_width::Float64
end


function DyadicLightweightExperiment(experiment::DyadicExperiment)

    DyadicLightweightExperiment(

        # params
        experiment.params,

        # estimator
        experiment.estimator.fhat,
        experiment.truth.f,
        experiment.estimator.psd_error,
        experiment.estimator.min_eig_Sigmahat,
        experiment.estimator.min_eig_Sigmahatplus,
        experiment.estimator.bandwidth_ROT,

        # diagnosis
        experiment.diagnosis.IMSEhat,
        experiment.diagnosis.MEhat,
        experiment.diagnosis.ucb_coverage,
        experiment.diagnosis.ucb_max_width,
        experiment.diagnosis.ucb_mean_width,
        experiment.diagnosis.bcb_coverage,
        experiment.diagnosis.bcb_max_width,
        experiment.diagnosis.bcb_mean_width,
        experiment.diagnosis.pcb_coverage,
        experiment.diagnosis.pcb_max_width,
        experiment.diagnosis.pcb_mean_width,
    )
end


function format_experiment(experiment::DyadicExperiment, params::Dict)

    if params["save_method"] == "full"
        return experiment
    elseif params["save_method"] == "lightweight"
        return DyadicLightweightExperiment(experiment)
    else
        error("unknown save_method")
    end
end


function run_dyadic_experiment(params::Dict)

    data = DyadicData(params)
    evals = DyadicEvaluationPoints(params)
    estimator = DyadicKernelDensityEstimator(params, data, evals)
    truth = DyadicTruth(params, evals)
    diagnosis = DyadicDiagnosis(truth, estimator, evals)
    experiment = DyadicExperiment(params, data, evals, estimator, truth, diagnosis)

    return format_experiment(experiment, params)
end



function run_dyadic_experiments(params_list::Vector{Dict{String, Any}})

    n_experiments = length(params_list)
    global counter = 0
    t00 = now()

    experiments = Vector{Any}(undef, n_experiments)
    Threads.@threads for rep in 1:n_experiments

        t0 = now()
        params = params_list[rep]
        experiments[rep] = run_dyadic_experiment(params)
        global counter += 1
        t1 = now()
        t_elapsed = canonicalize(t1 - t0)
        t_now = Time(t1)
        println("Experiment $counter / $n_experiments: finished at $t_now")
    end

    t1 = now()
    elapsed = canonicalize(t1 - t00)
    println("Experiments finished in ", elapsed)

    return experiments
end
