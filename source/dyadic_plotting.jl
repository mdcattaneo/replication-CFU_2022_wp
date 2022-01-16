using PyPlot
using Distributions

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["mathtext.fontset"]  = "cm";
rcParams["font.family"]  = "serif";
rcParams["font.serif"]  = "cm";
plt.ioff()

linewidth = 1.0
num_intervals_PCB = 10



function plot_experiment(
    experiment::DyadicExperiment,
    ylim, filepath::String)

    fig, ax = plt.subplots(figsize=(4,4))

    e = 0.22
    pcb_marker = PyPlot.matplotlib.path.Path([[-1,0], [-1,e], [-1,-e], [-1,0], [1,0], [1,e], [1,-e], [1,0]])
    handle_f = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth, label="\$f_W(w)\$")
    handle_fhat = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth, linestyle=(0, (1,1)), label="\$\\widehat f_W(w)\$")
    handle_ucb = PyPlot.matplotlib.patches.Patch(facecolor="lightgray", edgecolor="lightgray", label="UCB")
    handle_pcb = PyPlot.matplotlib.lines.Line2D([0], [0], color="black", lw=0, markersize=17, marker=pcb_marker, label="PCI")
    handles = [handle_f, handle_fhat, handle_ucb, handle_pcb]

    plot_ucb(ax, experiment)
    plot_pcb(ax, experiment)
    plot_f(ax, experiment)
    plot_fhat(ax, experiment)

    # format and save
    PyPlot.xlabel("Evaluation point, \$w\$")
    plt.ylim(ylim)
    plt.yticks(range(0.0, stop=0.4, step=0.1))
    x_min = minimum(experiment.evals.w)
    x_max = maximum(experiment.evals.w)
    plt.xlim((x_min, x_max))
    legend(handles=handles, loc="upper left")
    plt.ylabel("Density", labelpad=4.0)
    plt.tight_layout()
    PyPlot.savefig(filepath)
    close("all")

end



function plot_f(ax, experiment)

    ax.plot(
        experiment.evals.w,
        experiment.truth.f,
        color = "black",
        linewidth=linewidth,
    )

end



function plot_fhat(ax, experiment)

    ax.plot(
        experiment.evals.w,
        experiment.estimator.fhat,
        color = "black",
        linewidth=linewidth,
        linestyle=(0, (1,1)),
    )

end



function plot_ucb(ax, experiment)

    plot_confidence_band(
        ax,
        experiment.evals.w,
        experiment.estimator.ucb[1,:],
        experiment.estimator.ucb[2,:],
        "lightgray",
    )

end



function plot_pcb(ax, experiment)

    plot_confidence_intervals(
        ax,
        experiment.evals.w,
        experiment.estimator.pcb[1,:],
        experiment.estimator.pcb[2,:],
        num_intervals_PCB,
        "black",
        linewidth,
        "-"
    )

end



function plot_confidence_band(ax, x, y_lower, y_upper, fill_color)

    ax.fill_between(
        x, y_lower, y_upper,
        color=fill_color,
        linewidth=0.0,
    )

end



function plot_confidence_intervals(ax, x, y_lower, y_upper, num_intervals,
                                   line_color, line_width, line_style)

    inds = range(1, stop=length(x), length=2*num_intervals+1)
    inds = [inds[i] for i in 1:length(inds) if i%2 == 0]
    inds = Int.(round.(inds))
    @assert size(inds, 1) == num_intervals
    xr = x[inds]
    ylr = y_lower[inds]
    yur = y_upper[inds]

    ax.vlines(
        xr[1], ylr[1], yur[1],
        linewidth=line_width,
        linestyle=line_style,
        color=line_color,
    )

    for i in 1:num_intervals

        ax.vlines(
            xr[i], ylr[i], yur[i],
            linewidth=line_width,
            linestyle=line_style,
            color=line_color,
        )

        marker = PyPlot.matplotlib.path.Path([[-1,0], [1,0], [-1,0], [1,0]])

        ax.plot(
            xr[i], ylr[i],
            linewidth=0.0,
            color=line_color,
            marker=marker,
            markersize=4.0,
        )

        ax.plot(
            xr[i], yur[i],
            linewidth=0.0,
            color=line_color,
            marker=marker,
            markersize=4.0,
        )

    end

end
