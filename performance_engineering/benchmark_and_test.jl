using JSON


include("benchmark_getXBubble.jl")
include("../test/regression/dimer_anisotropy/regression_tests_dimer.jl")
include("git_utils.jl")
import .GitUtils

function add_data_to_db(data)
    fname = joinpath(@__DIR__, "benchmark_getXBubble.db")
    existing_data = open(fname, "r") do f
        JSON.parse(read(f))
    end
    push!(existing_data, data)
    open(fname, "w") do f
        JSON.json(f, existing_data; pretty=true)
    end

end

function getXBubble_test_and_benchmark()
    println("Testing...")
    run_getXbubble_regression_tests()
    println("Benchmarking...")
    bench_result = benchmark_synthetic_square(N=10, lattice_size=5)

    funcnames = ["mean","minimum","maximum"]
    quantities = ["times","gctimes"]


    data = Dict("commit" => GitUtils.get_git_commit_short(),
                      "nthreads" => Threads.nthreads(),
                      "benchmark_data" => Dict("$(q)_$(fn)"=> eval(Meta.parse("$fn($(getfield(bench_result, Symbol(q))))"))
                for q in quantities
                    for fn in funcnames))
    add_data_to_db(data)
    return bench_result
end
