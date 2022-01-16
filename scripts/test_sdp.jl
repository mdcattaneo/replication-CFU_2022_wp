using Convex
using Mosek
using MosekTools
using SCS
using BenchmarkTools
using MathOptInterface
const MOI = MathOptInterface

n = 100

A = randn((n,n))
A = A + A'
H = rand(Float64, (n,n))
H = H + H'

for j in 1:5

    B = A + randn(n,n) / 20
    C = Semidefinite(size(B, 1))
    objective = maximum(abs((H .* (C - B))))
    problem = minimize(objective)

    println("Mosek")
    @btime solve!($problem, Mosek.Optimizer, silent_solver=true, warmstart=true)
    println("SCS")
    @btime solve!($problem,
                 #MOI.OptimizerWithAttributes(SCS.Optimizer, "warm_start" => true),
                 SCS.Optimizer,
                 silent_solver=true)

end
