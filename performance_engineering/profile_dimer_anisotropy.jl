#!/usr/bin/env julia
# Profiling script for dimer anisotropy example using PProf.jl
# Based on pmfrg_xyz/example-dimer-anisotropy.sh
# Activate parent project to access dependencies
using Pkg
Pkg.activate(@__DIR__)


using PProf
using Profile

# Create output directory for profile data
profile_output_dir = joinpath(@__DIR__, "profile_data")
mkpath(profile_output_dir)

# Get current git commit hash (first 6 characters)
function get_git_commit_short()
    try
        commit_hash = strip(read(`git rev-parse HEAD`, String))
        return commit_hash[1:6]
    catch
        @warn "Could not get git commit hash, using 'nogit' as fallback"
        return "nogit"
    end
end

# Check for uncommitted changes
function check_git_status()
    try
        status_output = strip(read(`git status --porcelain`, String))
        if !isempty(status_output)
            @warn """
            ⚠️  WARNING: There are uncommitted changes in the repository!
            The profile will be named after commit $(get_git_commit_short()),
            but your current code may differ from that commit.

            Uncommitted changes detected:
            $(status_output)

            Consider committing your changes before profiling for accurate tracking.
            """
            return false
        end
        return true
    catch
        @warn "Could not check git status"
        return true  # Don't block profiling if git check fails
    end
end

using SpinFRGLattices
using PMFRG_xyz

# Check for uncommitted changes
check_git_status()

println("Setting up dimer anisotropy system...")
system = SpinFRGLattices.getPolymer(2)
par = Params(system)
isotropy = zeros(3, length(system.couplings))
isotropy .= [1.0, 0.5, 0.2]
isotropy_matrix = Matrix(transpose(isotropy))

println("Running initial execution (for compilation)...")
SolveFRG(par, isotropy_matrix)

println("\nStarting profiling run...")
Profile.clear()
@profile SolveFRG(par, isotropy_matrix)

# Generate filename with git commit hash
git_commit = get_git_commit_short()
profile_file = joinpath(profile_output_dir, "profile_dimer_anisotropy_$(git_commit).pb.gz")

println("\nSaving profile data to: $profile_file")
pprof(out=profile_file)

println("\nGenerating interactive profile visualization...")
pprof()  # Opens interactive pprof viewer in browser

println("\nProfile data saved to: $profile_file")
println("Git commit: $git_commit")
println("\nTo view saved profile later, run:")
println("  julia -e 'using PProf; pprof(\"$profile_file\")'")
println("Or use the pprof command-line tool:")
println("  pprof -http=: \"$profile_file\"")
