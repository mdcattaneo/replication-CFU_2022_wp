# local script
using Revise
using JLD
using Dates

# load source code
include("../source/dyadic_kde.jl")


# load data
printstyled("Loading data\n", color=:green)
#experiments_table = load("../../data/data_table.jld")["experiments_table"]
#println("Found $(length(experiments_table)) experiments for tables")
experiments_plot = load("../../data/data_plot.jld")["experiments_plot"]
println("Found $(length(experiments_plot)) experiments for plots")


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


# illustration plots
#println("Making illustration plots")
#include("illustration_plots.jl")


# tables
#println("Making tables")
#include("tables.jl")
#make_table(experiments_table, "../../tables/")
