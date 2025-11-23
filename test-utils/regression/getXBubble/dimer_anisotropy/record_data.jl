#!/usr/bin/env julia
# Script to record regression test data for getXBubble! function
#
using Recorder
using SpinFRGLattices
import PMFRG_xyz: Params, SolveFRG

function print_usage()
    println( "1. Apply the patch: git apply test-utils/patches/add_recorder_getxbubble.patch")
    println( "2. Add Recorder.jl to the project")
    println( "2. Run this script: julia --project=test-utils test/regression/getXBubble/dimer_anisotropy/record_data.jl")
    println( "3. Unapply patch: git apply -R < test-utils/patches/add_recorder_getxbubble.patch")
end
    
function main()::Int
    record_test_data()
    script_fname = generate_regression_tests()
    print_summary(script_fname)
    return 0
end

# level 1
function record_test_data()
    println("\nRunning dimer example with recording enabled...")
    println("This will record every 5th call to getXBubble! (calls 1, 6, 11, ..., 46)")

    @time run_dimer_example()

    # Get the first key from return_values (should be getXBubble!)
    keys_list = collect(keys(Recorder.gs.return_values))
    if !isempty(keys_list)
        call_count = length(Recorder.gs.return_values[keys_list[1]])
        println("\nRecorded $call_count calls to $(keys_list[1])")
    else
        println("\nWarning: No calls were recorded!")
        println("\nHave you applied the patch correctly?")
    end
end


function generate_regression_tests()
    println("\nGenerating regression test files...")

    current_dir = pwd()


    try
        cd(@__DIR__)
        fname, _  = Recorder.create_regression_tests(tag="getXBubble_dimer")
        return joinpath(@__DIR__,fname)
    finally
        cd(current_dir)
    end
end

function print_summary(script_fname)
    data_file = joinpath(@__DIR__, "regression_tests_getXBubble_dimer.data")
    test_file = script_fname 

    println("\n" * "="^80)
    println("Regression test data recorded successfully!")
    println("="^80)
    println("\nGenerated files:")
    println("  Data: $data_file")
    println("  Test: $test_file")
    println("\nNext steps:")
    println("  1. Unapply the patch:")
    println("     cd /home/michele/PMFRG/pmfrg_xyz")
    println("     git apply -R  < test-utils/patches/add_recorder_getxbubble.patch")
    println("\n  2. Review and customize the generated test file")
    println("     (especially floating-point comparisons)")
    println("\n  3. Run tests:")
    println("     julia --project=test-utils $test_file")
    println("="^80)
end

# level 2

function run_dimer_example()
    @eval begin
        System = SpinFRGLattices.getPolymer(2)
        par = Params(System,
                    N=8,
                    accuracy = 1e-10,
                    temp_max = 10.0,
                    temp_min = 1.0
                    )
        isotropy = zeros(System.Npairs, 3)

        for n in 1:System.Npairs
            isotropy[n, :] = [1.0, 0.5, 0.2]
        end

        results = SolveFRG(par,isotropy)
    end
end

#######

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
