include("benchmark_getXBubble.jl")
include("../test/regression/dimer_anisotropy/regression_tests_dimer.jl")

function getXBubble_test_and_benchmark()
    run_getXbubble_regression_tests()
    benchmark_synthetic_square(N=10,lattice_size=5)
end
