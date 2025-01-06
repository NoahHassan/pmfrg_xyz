import Pkg
Pkg.add("JLD2")
using JLD2

sol_yannik = load_object("dimer_flow_yannik.jld2")
sol_noah   = load_object("dimer_flow_noah.jld2")

sol_yannik[15][1]
sol_noah[15][1]

total = 0
n_zeros = 0
n_ones = 0
n_negones = 0
for val in sol_noah[1][1]
    total += 1
    if(val == 0.0)
        n_zeros += 1
    elseif(val == 1.0)
        n_ones += 1
    elseif(val == -1.0)
        n_negones += 1
    end
end

println(total)
println(n_zeros)
println(n_ones)
println(n_negones)
println(n_zeros + n_ones + n_negones)