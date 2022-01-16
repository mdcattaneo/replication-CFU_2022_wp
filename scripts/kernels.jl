using UnicodePlots

include("../source/dyadic_helpers.jl")
include("../source/dyadic_data_generation.jl")
include("../source/dyadic_estimation.jl")

a = -1.0
b = 2.0
h = 1.0
w = 0.0
kernel_name = "epanechnikov_order_4"

for w in range(a, stop=b, length=100)
    ss = range(a, stop=b, length=1000)
    ks = [kernel(s,w,h,a,b,kernel_name) for s in ss]
    println(lineplot(ss, ks))
    #i0 = trapezium_integrate(ks, collect(ss))
    #i1 = trapezium_integrate((ss .- w) .* ks, collect(ss))
    sleep(0.1)
end
