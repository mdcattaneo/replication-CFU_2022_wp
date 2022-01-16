using Printf

Base.@kwdef struct DyadicData

    n_data::Int
    N_data::Int = n_data * (n_data-1) // 2
    A::Vector{Float64} = fill(NaN, (n_data))
    V::Array{Float64, 2} = fill(NaN, (n_data, n_data))
    W::Array{Float64, 2} = fill(NaN, (n_data, n_data))
    W_vec::Vector{Float64} = Float64[]

end



function DyadicData(params::Dict)

    data = DyadicData(n_data = params["n_data"])

    if params["data_distribution"] == "triplet_gaussian_sum_product"
        generate_A_triplet(data, params["data_p"])
        generate_V_gaussian(data)
        generate_W_sum_product(data)
    else
        error("unknown data_distribution")
    end

    return data
end



function generate_A_triplet(data::DyadicData, p::Vector{Float64})

    @assert length(p) == 3
    @assert all(0 .<= p .<= 1)
    @assert sum(p) == 1

    c = cumsum(p)
    r = rand(data.n_data)
    data.A .= -1.0 * (r .<= c[1]) + 1.0 * (r .>= c[2])
end



function generate_V_gaussian(data::DyadicData)

    for i in 1:data.n_data
        for j in 1:data.n_data
            if i < j
                data.V[i,j] = randn()
            else
                data.V[i,j] = NaN
            end
        end
    end

    return nothing
end



function generate_W_sum_product(data::DyadicData)

    for i in 1:data.n_data
        for j in 1:data.n_data
            if i < j
                W = data.A[i] * data.A[j] + data.V[i,j]
                data.W[i,j] = W
                push!(data.W_vec, W)
            else
                data.W[i,j] = NaN
            end
        end
    end

    return nothing
end
