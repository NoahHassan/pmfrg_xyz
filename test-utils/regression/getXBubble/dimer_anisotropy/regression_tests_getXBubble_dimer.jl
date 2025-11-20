using PMFRG_xyz
import PMFRG_xyz: getXBubble!
using Test
using JLD2
@testset verbose = true "Tests for getXBubble!" begin
    function compare_arguments_post(args_post_exp, args_post)
        # getXBubble! only modifies Workspace.X (first argument)
        # Extract the Workspace from the argument lists
        workspace_exp = args_post_exp[1]
        workspace = args_post[1]

        # Compare only the X field
        X_exp = workspace_exp.X
        X = workspace.X

        result = true
        if X != X_exp
            abs_diff = abs.(X .- X_exp)
            total_diff = sum(abs_diff)
            max_diff = maximum(abs_diff)

            # Use tolerance of 1e-14 (same as PMFRGCore reference)
            if total_diff > 1e-14
                different_val_places = X .!= X_exp

                println("Test failed: Workspace.X differs significantly")
                println("  Absolute Difference (total): ", total_diff)
                println("  Maximum Difference: ", max_diff)
                println("  Number of different values: ", sum(different_val_places))

                println("  First computed values: ")
                println("    ", first(X[different_val_places], 5))
                println("  First expected values: ")
                println("    ", first(X_exp[different_val_places], 5))

                result = false
            end
        end
        result
    end

    data = load_object((@__DIR__) * "/regression_tests_getXBubble_dimer.data")
    @testset for i = 1:length(data["return_value"])
        return_value = (data["return_value"])[i]
        arguments = (data["arguments"])[i]
        arguments_post = (data["arguments_post"])[i]
        # Call the function to mutate the arguments
        getXBubble!(arguments...)
        # Check that Workspace.X matches the expected post-call state
        @test compare_arguments_post(arguments_post, arguments)
    end
end
