using Revise
#using UnicodePlots
#using ProfileView

println("Loading dyadic source code...")
include("../source/dyadic_helpers.jl")
include("../source/dyadic_data_generation.jl")
include("../source/dyadic_estimation.jl")
include("../source/dyadic_evaluation_points.jl")
include("../source/dyadic_truth.jl")
include("../source/dyadic_diagnosis.jl")
include("../source/dyadic_experiments.jl")
include("../source/dyadic_plotting.jl")

n_data = 500
n_evals = 100

params = Dict(
    "n_data" => n_data,
    "data_distribution" => "triplet_gaussian_sum_product",
    "data_p" => [0.2, 0.4, 0.4],
    "n_evals" => n_evals,
    "evals_type" => "deterministic",
    "evals_w_min" => -2.0,
    "evals_w_max" => 2.0,
    "kernel_name" => "epanechnikov_order_4",
    "bandwidth" => "ROT",
    "significance_level" => 0.05,
    "n_resample" => 10000,
    "verbose" => false,
    "psd_method" => "sdp",
    "fit_type" => "full",
    "save_method" => "lightweight",
)

n_repeats = 1
experiments = run_dyadic_experiments([copy(params) for _ in 1:n_repeats])
e1 = experiments[1]

println()
#average_fhat = zeros(n_evals)
#for i in 1:n_evals
    #average_fhat[i] = mean([e.fhat[i] for e in experiments])
#end

#bias = average_fhat .- e1.truth.f
#integrated_squared_bias = mean(bias.^2)
#println(params["kernel_name"])
#println(integrated_squared_bias)
#println(UnicodePlots.lineplot(bias))
#println(UnicodePlots.lineplot(experiments[1].truth.bias_fhat))

#println(UnicodePlots.lineplot(experiments[1].ucb))
#println(UnicodePlots.lineplot(experiments[1].fhat))
#println(UnicodePlots.lineplot(experiments[1].f))
#println(UnicodePlots.lineplot(experiments[1].truth.asymp_bias_fhat))


#println(UnicodePlots.lineplot(experiments[1].truth.Var_fhat))
#println(UnicodePlots.lineplot(experiments[1].truth.asymp_Var_fhat))

#display(experiments[1])


#display(experiments[1].truth.bias_fhat - experiments[1].truth.asymp_bias_fhat)

#println()
#print("ucb: ")
#println(sum([e.ucb_coverage for e in experiments]) / n_repeats)
#print("bcb: ")
#println(sum([e.bcb_coverage for e in experiments]) / n_repeats)
#print("pcb: ")
#println(sum([e.pcb_coverage for e in experiments]) / n_repeats)


#ylim = [0, 0.5]
#ordering = 1:5
#filepath = "../../plots/temp.pdf"
#for i in 1:5
    #plot_lightweight_experiment(experiments[i], ylim, ordering, filepath)
    #sleep(2)
#end


#println(experiments)

#println("ROT:   ", sum([e.estimator.bandwidth_ROT for e in experiments]) / n_repeats)
#println("IMSE:  ", experiments[1].truth.bandwidth_IMSE_optimal)
#println("AIMSE: ", experiments[1].truth.bandwidth_AIMSE_optimal)

#display(experiments[1])
#display(experiments[1].estimator.fhat)


#experiments = Vector{Any}(undef, n_repeats)
#for rep in 1:n_repeats
    #experiments[rep] = run_dyadic_experiment(params)
    #println("\rRunning experiment $rep / $n_repeats     ")
    #print("ucb: ")
    #println(sum([experiments[r].diagnosis.ucb_coverage for r in 1:rep]) / rep)
    #print("bcb: ")
    #println(sum([experiments[r].diagnosis.bcb_coverage for r in 1:rep]) / rep)
    #print("pcb: ")
    #println(sum([experiments[r].diagnosis.pcb_coverage for r in 1:rep]) / rep)
#end







#=
# get oracle variance estimator
experiments_oracle = run_dyadic_experiments([params for _ in 1:n_repeats])
Sigmahatplusses = [e.estimator.Sigmahatplus for e in experiments_oracle]
oracle_Sigmahatplus = sum(Sigmahatplusses) / n_repeats
min_eig_oracle_Sigmahatplus = minimum(eigen(oracle_Sigmahatplus).values)

# use oracle variance estimator for CBs
params["fit_type"] = "partial"
experiments = Vector{Any}(undef, n_repeats)
global counter = 0

for rep in 1:n_repeats
    data = DyadicData(params)
    evals = DyadicEvaluationPoints(params)
    estimator = DyadicKernelDensityEstimator(params)
    fit(params, estimator, data, evals)
    estimator.Sigmahatplus = oracle_Sigmahatplus
    estimate_uniform_confidence_band(estimator, data, evals)
    estimate_bonferroni_confidence_band(estimator, data, evals)
    estimate_pointwise_confidence_band(estimator, data, evals)
    truth = DyadicTruth(params, evals)
    diagnosis = DyadicDiagnosis(truth, estimator, evals)
    experiment = DyadicExperiment(params, data, evals, estimator, truth, diagnosis)
    experiments[rep] = experiment
    global counter += 1
    print("\rRunning experiment $counter / $n_repeats     ")
end

println()
print("ucb: ")
println(sum([e.diagnosis.ucb_coverage for e in experiments]) / n_repeats)
print("bcb: ")
println(sum([e.diagnosis.bcb_coverage for e in experiments]) / n_repeats)
print("pcb: ")
println(sum([e.diagnosis.pcb_coverage for e in experiments]) / n_repeats)
=#








#Sigmahatplusses = [e.estimator.Sigmahatplus for e in experiments]
#oracle_Sigmahatplus = sum(Sigmahatplusses) / n_repeats
#oracle_Sigmahatplus_diag = diag(oracle_Sigmahatplus)

#e = experiments[1]
#println()

#println("oracle:")
#display(n_data.^2 * bandwidth * sum(oracle_Sigmahatplus_diag) / n_evals)
#println()
#println("estimator:")
#display(n_data.^2 * bandwidth * sum(diag(e.estimator.Sigmahatplus)) / n_evals)
#println()
#println("truth:")
#display(n_data.^2 * bandwidth * sum(e.truth.Var_fhat) / n_evals)
#println()
#display(e.truth.Var_fhat ./ oracle_Sigmahatplus_diag)
#println("empirical:")
#empirical_Sigma_diag = sum([(e.estimator.fhat - e.truth.f).^2 for e in experiments]) / n_repeats
#display(n_data.^2 * bandwidth * sum(empirical_Sigma_diag) / n_evals)
#println()

#println("scaled oracle error:")
#println(n_data^1.5 * maximum(abs.(oracle_Sigmahatplus_diag - e.truth.Var_fhat)))
#println()

#println("scaled estimator error:")
#println(n_data^1.5 * maximum(abs.(diag(e.estimator.Sigmahatplus) - e.truth.Var_fhat)))
#println()







#ylim = [-0.1, 1]
#ordering = [1,2,3,4,5]
#filepath = "../../plots/temp.pdf"
#plot_lightweight_experiment(e, ylim, ordering, filepath)

#println(e.truth.asymp_Var_fhat .^ 0.5)
#println(e.Sigmahatplus_diag .^ 0.5)
#println(e.Sigmahat_diag .^ 0.5)

#println(sqrt(n_data) * maximum(abs.(e.Sigmahatplus_diag .^ 0.5)))
#println(sqrt(n_data) * maximum(abs.(e.truth.Var_fhat .^ 0.5)))

#println(sqrt(n_data) * maximum(abs.(e.Sigmahatplus_diag .^ 0.5 - e.truth.Var_fhat .^ 0.5)))

#for e in experiments
    #println(e.diagnosis.ucb_coverage)
#end
#@profview run_dyadic_experiments(params_list)

#truth = experiments[1].truth
#evals = experiments[1].evals

#generate_IMSE_AIMSE_optimal_bandwidth(truth, evals, params)
