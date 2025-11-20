#!/usr/bin/env julia
# Regression test runner for getXBubble! function
#
# This script loads and runs the generated regression tests
# with appropriate floating-point comparison tolerances.

using Test
using Pkg

function main()::Int
    setup_test_environment()
    run_regression_tests()
    return 0
end

# level 1
function setup_test_environment()
    println("Setting up test environment...")
    test_dir = @__DIR__
    pmfrg_xyz_root = joinpath(test_dir, "../../..")
    examples_dir = joinpath(pmfrg_xyz_root, "examples")

    Pkg.activate(examples_dir)

    @eval using SpinFRGLattices
    @eval import PMFRG_xyz: Params, SolveFRG
    @eval using Recorder
end

function run_regression_tests()
    dimer_test_dir = joinpath(@__DIR__, "dimer_anisotropy")
    test_file = joinpath(dimer_test_dir, "regression_tests_getXBubble_dimer.jl")

    if !isfile(test_file)
        error("""
        Regression test file not found: $test_file

        Please generate test data first by running:
          1. Apply patch: cd pmfrg_xyz && patch -p1 < test/patches/add_recorder_getxbubble.patch
          2. Record data: julia --project=examples test/regression/getXBubble/dimer_anisotropy/record_data.jl
          3. Unapply patch: cd pmfrg_xyz && patch -R -p1 < test/patches/add_recorder_getxbubble.patch
        """)
    end

    println("\nRunning regression tests for getXBubble!...")
    println("Test file: $test_file\n")

    @testset "getXBubble! Regression Tests" begin
        include(test_file)
    end
end

#######

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
