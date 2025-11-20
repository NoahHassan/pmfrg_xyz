module PProfVersionUtils

import Pkg
using PProf

export get_pprof_version, get_pprof_path

function get_pprof_version()
    deps_dict = Pkg.dependencies()
    pprof_jll_uuid = find_pprof_jll_uuid(deps_dict)
    extract_version_string(deps_dict, pprof_jll_uuid)
end

# level 1
function find_pprof_jll_uuid(deps_dict)
    findfirst(p -> p.name == "pprof_jll", deps_dict)
end

function extract_version_string(deps_dict, uuid)
    version = deps_dict[uuid].version
    replace(string(version), "+" => "_")
end

function get_pprof_path()
    PProf.pprof_jll.pprof_path
end

end # module
