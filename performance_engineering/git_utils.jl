module GitUtils

export get_git_commit_short, check_git_status

function get_git_commit_short()
    try
        commit_hash = strip(read(`git rev-parse HEAD`, String))
        return commit_hash[1:6]
    catch
        @warn "Could not get git commit hash, using 'nogit' as fallback"
        return "nogit"
    end
end

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
        return true
    end
end

end
