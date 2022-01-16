using Revise
using PyPlot

include("../source/dyadic_kde.jl")

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["mathtext.fontset"]  = "cm";
rcParams["font.family"]  = "serif";
rcParams["font.serif"]  = "cm";
plt.ioff()




function make_illustration_plot(p, w_min, w_max, n_evals, filepath::String)

    p1 = p[1]; p2 = p[2]; p3 = p[3]
    linewidth = 1.0

    @assert 0 <= p1 <= 1
    @assert 0 <= p2 <= 1
    @assert 0 <= p3 <= 1
    @assert p1 + p2 + p3 == 1.0

    params = Dict(
        "data_distribution" => "triplet_gaussian_sum_product",
        "data_p" => p,
        "n_evals" => n_evals,
        "evals_type" => "deterministic",
        "evals_w_min" => w_min,
        "evals_w_max" => w_max,
    )

    evals = DyadicEvaluationPoints(params)
    truth = DyadicTruth(params, evals)

    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(evals.w, truth.f, color="black", label="\$f_W(w)\$",
            lw=linewidth)
    ax.plot(evals.w, truth.Var_f_given_A.^0.5,
            color="black", linestyle=(0, (1, 1)),
            lw=linewidth,
            label="\$\\mathrm{Var}[f_{W \\mid A}(w \\mid A_i)]^{1/2}\$")
    plt.ylim((-0.01, 0.45))
    plt.xlabel("Evaluation point, \$w\$")
    plt.yticks(range(0.0, stop=0.4, step=0.1))
    legend(loc="upper left")
    plt.ylabel("Density", labelpad=4.0)
    plt.tight_layout()
    PyPlot.savefig(filepath)
    close("all")

end



function phi(x)

    return (2 * pi)^(-1/2) * exp(-(x^2) / 2)

end



# run plots
# ---------

w_min = -2.0
w_max = 2.0
n_evals = 1000

p_total = [0.5, 0.0, 0.5]
p_partial = [0.25, 0.0, 0.75]
p_none = [0.2, 0.2, 0.6]

make_illustration_plot(p_total, w_min, w_max, n_evals,
                       "../../plots/illustration_total.pdf")

make_illustration_plot(p_partial, w_min, w_max, n_evals,
                       "../../plots/illustration_partial.pdf")

make_illustration_plot(p_none, w_min, w_max, n_evals,
                       "../../plots/illustration_none.pdf")
