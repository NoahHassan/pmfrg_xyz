using JSON


include("benchmark_getXBubble.jl")
include("../test/regression/dimer_anisotropy/regression_tests_dimer.jl")



function getXBubble_test_and_benchmark()
    println("Testing...")
    run_getXbubble_regression_tests()
    println("Benchmarking...")
    bench_result = benchmark_synthetic_square(N=10,lattice_size=5)
    open(joinpath(@__DIR__,"benchmark_getXBubble.db"),"w") do f

        write(f,JSON.json(bench_result))
        end
    return bench_result
end
