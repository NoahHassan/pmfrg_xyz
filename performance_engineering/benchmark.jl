#!/usr/bin/env julia
# Usage: julia benchmark.jl [example_name]

using Pkg
Pkg.activate(@__DIR__)

using BenchmarkTools, PMFRG_xyz

include("example_configs.jl")
using .ExampleSetups: example_setups

include("git_utils.jl")
using .GitUtils: get_git_commit_short, check_git_status

# level 0
function main()
    setup_output_dir()
    example_name = parse_args()
    check_git_status()

    par, isotropy = setup_example(example_name)
    run_compilation(par, isotropy)
    benchmark_results = run_benchmark(par, isotropy)
    display_results(benchmark_results)
    save_results(example_name, benchmark_results)
end

# level 1
function setup_output_dir()
    mkpath(joinpath(@__DIR__, "benchmark_data"))
end

function parse_args()
    length(ARGS) >= 1 ? ARGS[1] : "dimer"
end

function setup_example(example_name)
    println("Setting up $example_name example...")
    example_setups[example_name]()
end

function run_compilation(par, isotropy)
    println("Running initial execution (for compilation)...")
    SolveFRG(par, isotropy)
end

function run_benchmark(par, isotropy)
    println("\nBenchmarking SolveFRG...")
    @benchmark SolveFRG($par, $isotropy)
end

function display_results(benchmark_results)
    println("\nBenchmark results:")
    display(benchmark_results)
    println("\n")
    print_timing_summary(benchmark_results)
end

function save_results(example_name, benchmark_results)
    git_commit = get_git_commit_short()
    benchmark_file = get_benchmark_filename(example_name, git_commit)

    write_benchmark_file(benchmark_file, example_name, git_commit, benchmark_results)
    println("\nBenchmark results saved to: $benchmark_file")
end

# level 2
function print_timing_summary(b)
    println("Minimum time: ", minimum(b.times) / 1e9, " seconds")
    println("Median time:  ", median(b.times) / 1e9, " seconds")
    println("Mean time:    ", mean(b.times) / 1e9, " seconds")
    println("Maximum time: ", maximum(b.times) / 1e9, " seconds")
end

function get_benchmark_filename(example_name, git_commit)
    joinpath(@__DIR__, "benchmark_data", "benchmark_$(example_name)_$(git_commit).txt")
end

function write_benchmark_file(filepath, example_name, git_commit, b)
    open(filepath, "w") do io
        println(io, "Benchmark for $example_name example")
        println(io, "Git commit: $git_commit")
        println(io, "\nBenchmark results:")
        println(io, b)
        println(io, "\nSummary:")
        println(io, "Minimum time: ", minimum(b.times) / 1e9, " seconds")
        println(io, "Median time:  ", median(b.times) / 1e9, " seconds")
        println(io, "Mean time:    ", mean(b.times) / 1e9, " seconds")
        println(io, "Maximum time: ", maximum(b.times) / 1e9, " seconds")
    end
end

main()
