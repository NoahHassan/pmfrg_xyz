using Test
using PMFRG_xyz

include("regression/dimer_anisotropy/regression_tests_dimer.jl")

@time @testset verbose=true "PMFRG_xyz tests" begin
    run_regression_tests()
end

