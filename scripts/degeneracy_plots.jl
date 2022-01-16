using PyPlot

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["mathtext.fontset"]  = "cm";
rcParams["font.family"]  = "serif";
rcParams["font.serif"]  = "cm";
plt.ioff()



function make_illustration_plot(p, w_min, w_max, n_points, filepath::String)

    q = 1-p
    ws = range(w_min, stop=w_max, length=n_points)
    f = Vector{Float64}(undef, n_points)
    Var_f_given_A = Vector{Float64}(undef, n_points)

    for i in 1:n_points
        w = ws[i]
        f[i] = (p^2 + q^2) * phi(w-1) + 2 * p * q * phi(w+1)
        Var_f_given_A[i] = p * q * (2*p-1)^2 * (phi(w-1) - phi(w+1))^2
    end

    fig, ax = plt.subplots(figsize=(5,4))
    ax.plot(ws, f, color="black", label="Density, \$f(w)\$")
    ax.plot(ws, Var_f_given_A.^0.5, color="black", linestyle=(0, (1, 1)),
            label="\$\\mathrm{Var}[f_{W \\mid A}(w \\mid A_i)]^{1/2}\$")
    xlabel("Evaluation point, \$w\$")
    ylabel("Density")
    plt.ylim((-0.01, 0.32))
    legend(loc="upper left")
    savefig(filepath)
    close("all")

end



function phi(x)

    return (2 * pi)^(-1/2) * exp(-(x^2) / 2)

end



# run plots
# ---------

w_min = -3
w_max = 3
n_points = 1000

make_illustration_plot(0.25, w_min, w_max, n_points,
                       "../../plots/illustration_p_0.25.pdf")

make_illustration_plot(0.5, w_min, w_max, n_points,
                       "../../plots/illustration_p_0.5.pdf")
