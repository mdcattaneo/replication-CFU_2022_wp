using Printf






function make_table(experiments, out_path)

    n_data = experiments[1].params["n_data"]
    n_evals = experiments[1].params["n_evals"]
    psd_method = experiments[1].params["psd_method"]
    n_resample = experiments[1].params["n_resample"]
    date_time = experiments[1].params["date_time"]
    println("n_data_table: ", n_data)
    println("n_evals_table: ", n_evals)
    println("psd_method: ", psd_method)
    println("n_resample: ", n_resample)
    println("date_time: ", date_time)

    s = """
        \\begin{center}
            \\begin{tabular}{|c|c|c|c|c|cc|cc|}
                \\hline
                \\multirow{2}{*}{\$ \\pi \$}
                & \\multirow{2}{*}{Degeneracy type}
                & \\multirow{2}{*}{\$ h^*_{\\ROT} \$}
                & \\multirow{2}{*}{\$ p \$}
                & \\multirow{2}{*}{RIMSE}
                & \\multicolumn{2}{|c|}{UCB}
                & \\multicolumn{2}{|c|}{PCB} \\\\
                \\cline{6-9}
                & & & &
                & CR & AW
                & CR & AW \\\\
                \\hline
        """

    for degen in ["total", "partial", "none"]

        es_degen = [e for e in experiments if e.params["degen"] == degen]
        p = es_degen[1].params["data_p"]
        p_rationals = format_as_latex_fraction.(p)
        formatted_p = "\\left(" * join(p_rationals, ", ") * "\\right)"
        degen_name = uppercasefirst(degen)
        degen_name = replace(degen_name, "_" => " ")
        average_bandwidth_ROT = mean([e.bandwidth_ROT for e in es_degen])
        average_bandwidth_ROT = @sprintf("%.3f", average_bandwidth_ROT)

        s = s * """
            \\multirow{2}{*}{\$ $formatted_p \$}
            & \\multirow{2}{*}{$degen_name}
            & \\multirow{2}{*}{$average_bandwidth_ROT}
            """

        bandwidth_type = "ROT"
        kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]

        for n in 1:length(kernel_names)

            if n != 1
                s = s * """ & & """
            end

            # get all repeats corresponding to a single row
            kernel_name = kernel_names[n]
            es = [e for e in experiments if e.params["n_data"] == n_data]
            es = [e for e in es if e.params["degen"] == degen]
            es = [e for e in es if e.params["kernel_name"] == kernel_name]
            global n_repeats = length(es)

            row = make_table_row(es)

            for r in row
                s = s * " & " * r
            end

            s = s * """ \\\\
                    """
        end

        s = s * """
            \\hline
            """
    end


    s = s * """
        \\end{tabular}
        \\end{center}
        % n_data_table = $n_data
        % n_evals_table = $n_evals
        % psd_method = $psd_method
        % n_repeats = $n_repeats
        % n_resample = $n_resample
        % date_time = $date_time
        """

    open("$out_path/table.tex", "w") do file
        write(file, s)
    end

    open("$out_path/table_$date_time.tex", "w") do file
        write(file, s)
    end

    return nothing

end







function make_table_row(experiments::Vector{DyadicLightweightExperiment})


    e1 = experiments[1]

    # check most parameters are equal
    for e in experiments
        for f in keys(e1.params)
            if !(string(f) == "bandwidth")
                if e.params[f] != e1.params[f]
                    error("table: parameter $f must be equal")
                end
            end
        end
    end

    n_evals = e1.params["n_evals"]
    average_fhat = zeros(n_evals)

    for i in 1:n_evals
        average_fhat[i] = mean([e.fhat[i] for e in experiments])
    end

    kernel_order = e1.params["kernel_name"][end]

    average_IMSEhat = mean([e.IMSEhat for e in experiments])
    average_RIMSEhat = sqrt(average_IMSEhat)

    average_ucb_coverage = 100 * mean([e.ucb_coverage for e in experiments])
    average_ucb_mean_width = mean([e.ucb_mean_width for e in experiments])

    average_bcb_coverage = 100 * mean([e.bcb_coverage for e in experiments])
    average_bcb_mean_width = mean([e.bcb_mean_width for e in experiments])

    average_pcb_coverage = 100 * mean([e.pcb_coverage for e in experiments])
    average_pcb_mean_width = mean([e.pcb_mean_width for e in experiments])

    return [
        kernel_order,
        @sprintf("%.4f", average_RIMSEhat),
        @sprintf("%.1f", average_ucb_coverage) * "\\%",
        @sprintf("%.3f", average_ucb_mean_width),
        @sprintf("%.1f", average_pcb_coverage) * "\\%",
        @sprintf("%.3f", average_pcb_mean_width),
    ]

end
