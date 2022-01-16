# hypothesis testing
using Revise
using UnicodePlots

println("Loading source code...")
include("../source/dyadic_kde.jl")

# iid normal, unknown mean, variance 1
######################################

# data
n_data = 200
n_evals = 30
w_min = -2.0
w_max = 2.0
n_repeats = 100

# fit kde model
p = [0.5, 0.0, 0.5]
n_resample = 10000

params = Dict(
    "n_data" => n_data,
    "data_distribution" => "triplet_gaussian_sum_product",
    "data_p" => p,
    "n_evals" => n_evals,
    "evals_type" => "deterministic",
    "evals_w_min" => w_min,
    "evals_w_max" => w_max,
    "kernel_name" => "epanechnikov_order_4",
    "bandwidth" => "ROT",
    "bandwidth_type" => "ROT",
    "significance_level" => 0.05,
    "n_resample" => n_resample,
    "verbose" => false,
    "psd_method" => "sdp",
    "fit_type" => "full",
    "save_method" => "full",
)

accepts = zeros(n_repeats)

Threads.@threads for rep in 1:n_repeats
    println("Experiment ", rep, " / ", n_repeats)

    kde_experiment = run_dyadic_experiment(copy(params))

    # fit parametric model
    n_params = 51
    mu1s = range(-1.2, stop=-0.8, length=n_params)
    mu2s = range(0.8, stop=1.2, length=n_params)
    mixs = range(0.4, stop=0.6, length=n_params)
    abs_error = fill(Inf, (n_params, n_params, n_params))
    ws = kde_experiment.evals.w
    fhat_kde = kde_experiment.estimator.fhat
    sqrt_diag_Sigmahatplus = sqrt.(diag(kde_experiment.estimator.Sigmahatplus))

    for i in 1:n_params
        for j in 1:n_params
            if i <= j
                for r in 1:n_params
                    mix = mixs[r]
                    fi = mix .* phi.(ws .- mu1s[i])
                    fj = (1-mix) .* phi.(ws .- mu2s[j])
                    abs_error[i,j,r] = maximum(abs.(fi + fj - fhat_kde) ./ sqrt_diag_Sigmahatplus)
                end
            end
        end
    end

    best_idx = argmin(abs_error)
    println(best_idx)
    best_mu1 = mu1s[best_idx[1]]
    best_mu2 = mu2s[best_idx[2]]
    best_mix = mixs[best_idx[3]]
    f1 = best_mix .* phi.(ws .- best_mu1)
    f2 = (1-best_mix) .* phi.(ws .- best_mu2)
    println("best mu1: ", best_mu1)
    println("best mu2: ", best_mu2)
    println("best mix: ", best_mix)

    # accept or reject test
    ucb = kde_experiment.estimator.ucb
    accept = all(ucb[1,:] .<= f1 + f2 .<= ucb[2,:])
    accepts[rep] = accept
    println("accept? ", accept)
    println("acceptance rate: ", mean(accepts[1:rep]))

    display(UnicodePlots.lineplot(fhat_kde))
    display(UnicodePlots.lineplot(f1+f2))
end

println("acceptance rate: ", mean(accepts))
