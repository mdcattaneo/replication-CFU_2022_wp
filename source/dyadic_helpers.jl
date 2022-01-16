using Printf
using LinearAlgebra



"Use the trapezium rule to estimate integral f dx"
function trapezium_integrate(f::Vector{Float64}, x::Vector{Float64})

    if length(f) != length(x)
        error("vectors must be the same length")
    end

    n_areas = length(x) - 1
    areas = fill(NaN, n_areas)

    for i in 1:n_areas
        areas[i] = 0.5 * (f[i+1] + f[i]) * (x[i+1] - x[i])
    end

    area = sum(areas)

    return area
end



function get_quantile(A::Vector, q::Real)

    sorted_A = sort(A)
    idx = round(Int, length(A) * q)
    quantile = sorted_A[idx]

    return quantile
end



function phi(t::Float64)

    return (2 * pi)^(-0.5) * exp(-(t^2) / 2)
end



function phi_derivative(t::Float64, p::Int)

    if p == 0
        return phi(t)
    elseif p == 1
        return (-t) * phi(t)
    elseif p == 2
        return (t^2 - 1) * phi(t)
    elseif p == 3
        return (-t) * (t^2 - 3) * phi(t)
    elseif p == 4
        return (t^4 - 6*t^2 + 3) * phi(t)
    else
        error("unknown derivative order")
    end
end



function my_print(a, b)

    Printf.@printf("\n    %-25s %-5s %-10s", a, "=>", b)
    return nothing
end



function format_as_latex_fraction(p)

    if p == 0
        return "0"
    else
        p_rational = string(rationalize(p))
        p_rational = replace(p_rational, "//" => "}{")
        p_rational = "\\frac{" * p_rational * "}"
    end
end



"l2-project a vector v onto the l1-ball of radius b"
function l1proj(v::Vector{Float64}, b::Float64)

    # Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML

    @assert b > 0
    n = length(v)
    u = sort(abs.(v), rev=true)
    sv = cumsum(u)
    a = u .- (sv .- b) ./ (1:n)
    rho = maximum([j for j in 1:n if a[j] > 0])
    theta = max(0, (sv[rho] - b) / rho)
    w = sign.(v) .* max.(abs.(v) .- theta, 0)
    return w
end



"ADMM algorithm to find the nearest epsilon-PD matrix in maximum norm."
function ADMM(mat::Symmetric{Float64}, epsilon::Float64=1e-4, mu::Float64=10.0,
              it_max::Int=1000, etol::Float64=1e-4, etol_distance::Float64=1e-4,
              verbose::Bool=false)

    # mu: ADMM penalty parameter
    # it_max: maximum number of iterations
    # etol: tolerance parameter for convergence of primal and dual residual
    # etol_distance: tolerance parameter for convergence of the distance
    # https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf
    # https://github.com/celiaescribe/BDcocolasso

    t00 = time()
    p = size(mat, 1)
    R = Diagonal(mat)
    S = zeros((p,p))
    L = zeros((p,p))
    itr = 0

    while itr < it_max

        Rp = copy(R)
        Sp = copy(S)
        t0 = time()

        # Subproblem I: R step
        W = mat + S + mu * L
        W_eigdec = eigen(W)
        W_V = real.(W_eigdec.vectors)
        W_D = real.(W_eigdec.values)
        R = W_V * Diagonal(max.(W_D, epsilon)) * W_V'

        # Subproblem II: S step
        M = R - mat - mu * L

        # Lower triangle of M
        M_tril = zeros(Int(p*(p+1)/2))
        counter = 1
        for i in 1:p
            for j in 1:p
                if j <= i
                    M_tril[counter] = M[i,j]
                    counter +=1
                end
            end
        end

        l1proj_M = l1proj(M_tril, mu/2)

        # Update S
        counter = 1
        for i in 1:p
            for j in 1:p
                if j <= i
                    S[i,j] = M_tril[counter] - l1proj_M[counter]
                    S[j,i] = S[i,j]
                    counter +=1
                end
            end
        end

        # L step: update the Lagrange parameter
        L = L - (R - S - mat) / mu
        t1 = time()

        if verbose
            println("Iteration: ", itr)
            println("eps_R: ", maximum(abs.(R - Rp)))
            println("eps_S: ", maximum(abs.(S - Sp)))
            println("eps_primal: ", maximum(abs.(R - S - mat)))
            println("Iteration time: ", t1 - t0)
            println("Max norm distance: ", maximum(abs.(R - mat)))
        end

        # Check convergence
        if (((maximum(abs.(R - Rp)) < etol) &&
            (maximum(abs.(S - Sp)) < etol) &&
            (maximum(abs.(R - S - mat)) < etol)) ||
            (abs(maximum(abs.(Rp - mat)) - maximum(abs.(R - mat))) < etol_distance))
          itr = it_max
        else
          itr += 1
        end

        # Update learning rate
        if itr % 20 == 0
          mu /= 2
        end

  end

  if verbose
      println("Total time: ", time() - t00)
  end

  return Symmetric(R)

end



function matrix_lipschitz_number(mat::Symmetric{Float64}, v::Vector{Float64})

    @assert size(mat, 1) == length(v)
    steps = diff(mat, dims=1) ./ diff(v)
    return maximum(steps)
end
