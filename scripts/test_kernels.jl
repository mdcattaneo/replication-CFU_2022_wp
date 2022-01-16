include("../source/dyadic_helpers.jl")
include("../source/dyadic_data_generation.jl")
include("../source/dyadic_estimation.jl")

println("Conducting kernel tests")

h = 0.1
w_min = 0.0
w_max = 1.0
tol = 1e-4
len = 10000

ws = collect(range(0, stop=1, length=len))

# order 2
for i in 1:len
    w = ws[i]
    ks = [kernel(s, w, h, w_min, w_max, "epanechnikov_order_2") for s in ws]
    cs = (ws .- w) ./ h
    @assert abs(trapezium_integrate(ks, ws) - 1) <= tol
    @assert abs(trapezium_integrate(cs .* ks, ws) - 0) <= tol
    println(i, " / ", len, " tests passed")
end


# order 4
for i in 1:len
    w = ws[i]
    ks = [kernel(s, w, h, w_min, w_max, "epanechnikov_order_4") for s in ws]
    cs = (ws .- w) ./ h
    @assert abs(trapezium_integrate(ks, ws) - 1) <= tol
    @assert abs(trapezium_integrate(cs .* ks, ws) - 0) <= tol
    @assert abs(trapezium_integrate(cs.^2 .* ks, ws) - 0) <= tol
    @assert abs(trapezium_integrate(cs.^3 .* ks, ws) - 0) <= tol
    println(i, " / ", len, " tests passed")
end

println("Kernel tests passed")
