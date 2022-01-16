# cluster script
using Revise
using Random
using Dates
using JLD2

println("Loading source code...")
include("../source/dyadic_kde.jl")
date_time = round(Dates.now(), Dates.Minute)


# specify parameters
n_data_plot = 100
n_evals_plot = 100
n_data_table = 500
n_evals_table = 50
n_repeats = 2000
degeneracies = ["total", "partial", "none"]
kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
w_min = -2.0
w_max = 2.0
psd_method = "sdp"
verbose = false
n_resample = 10000
ps = Dict(
    "total" => [0.5, 0.0, 0.5],
    "partial" => [0.25, 0.0, 0.75],
    "none" => [0.2, 0.2, 0.6],
)


# build tables params list
printstyled("Getting parameters for tables\n", color=:green)
params_list_table = Dict{String, Any}[]

for degen in degeneracies
    for kernel_name in kernel_names

        params = Dict(
            "n_data" => n_data_table,
            "data_distribution" => "triplet_gaussian_sum_product",
            "data_p" => ps[degen],
            "n_evals" => n_evals_table,
            "evals_type" => "deterministic",
            "evals_w_min" => w_min,
            "evals_w_max" => w_max,
            "kernel_name" => kernel_name,
            "bandwidth" => "ROT",
            "bandwidth_type" => "ROT",
            "significance_level" => 0.05,
            "n_resample" => n_resample,
            "verbose" => verbose,
            "psd_method" => psd_method,
            "fit_type" => "full",
            "save_method" => "lightweight",
            "purpose" => "table",
            "degen" => degen,
            "date_time" => date_time,
        )

        for rep in 1:n_repeats
            push!(params_list_table, copy(params))
        end

    end
end





# build plots params list
printstyled("Generating parameters for plots\n", color=:green)
params_list_plot = Dict{String, Any}[]

for degen in degeneracies

    params = Dict(
        "n_data" => n_data_plot,
        "data_distribution" => "triplet_gaussian_sum_product",
        "data_p" => ps[degen],
        "n_evals" => n_evals_plot,
        "evals_type" => "deterministic",
        "evals_w_min" => w_min,
        "evals_w_max" => w_max,
        "kernel_name" => "epanechnikov_order_2",
        "bandwidth" => "ROT",
        "bandwidth_type" => "ROT",
        "significance_level" => 0.05,
        "n_resample" => n_resample,
        "verbose" => verbose,
        "psd_method" => psd_method,
        "fit_type" => "full",
        "save_method" => "full",
        "choose_optimal_bandwidth" => false,
        "purpose" => "plot",
        "degen" => degen,
        "date_time" => date_time,
    )

    push!(params_list_plot, copy(params))

end

# run plot experiments
Random.seed!(3142)
printstyled("Running $(length(params_list_plot)) experiments for plots\n", color=:green)
experiments_plot = run_dyadic_experiments(params_list_plot)
#JLD2.save("../../data/data_plot.jld", "experiments_plot", experiments_plot)

# run table experiments
#Random.seed!(3142)
#printstyled("Running $(length(params_list_table)) experiments for table\n", color=:green)
#experiments_table = run_dyadic_experiments(params_list_table)
#JLD2.save("../../data/data_table.jld", "experiments_table", experiments_table)

#printstyled("Experiments saved to disk\n", color=:green)

# experiments plots
println("Making experiment plots")
degeneracies = unique([e.params["degen"] for e in experiments_plot])
ylim = (-0.01, 0.45)
legend_ordering = [3, 4, 1, 2]

for degen in degeneracies
    for experiment in experiments_plot
        if experiment.params["purpose"] == "plot"
            if experiment.params["degen"] == degen
                plot_experiment(
                    experiment, ylim,
                    "../../plots/result_$(degen).pdf")
            end
        end
    end
end
