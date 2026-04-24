using Pkg
Pkg.activate(".")

include("../src/PMFRG_general.jl")
using .PMFRG_general

include("../src/PMFRG_xyz.jl")
using .PMFRG_xyz

include("../src/PMFRG_old.jl")
using .PMFRG_old

using SpinFRGLattices
using JLD2

using DiffEqCallbacks


function main()
    run_general()
    run_xyz()
    compare_general_xyz()
end


# Struct that contains all the relevant output 
# for the GENERAL implementation
struct Observables_Full{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
    iSigma::Array{T,3} # 9 Sigma-flavors
    Gamma::Array{T,5} # 81 Gamma-flavors
end

function run_general()

    function pack_all(State, t, Par)
        chix = PMFRG_general.getChi_x(State, exp(t), Par)
        chiy = PMFRG_general.getChi_y(State, exp(t), Par)
        chiz = PMFRG_general.getChi_z(State, exp(t), Par)
        println("t: $(round(exp(t);digits=3))")
        return Observables_Full(
            copy(chix),
            copy(chiy),
            copy(chiz),
            copy(State.x[2]),
            copy(State.x[3]),
        )
    end


    System = getPolymer(2)

    anisotropy = zeros(System.Npairs, 3, 3)
    for Rij = 1:System.Npairs
        anisotropy[Rij, :, :] = [
            -0.2 0.0 0.0
            0.0 3.0 0.0
            0.0 0.0 1.0
        ]
    end

    Par = PMFRG_general.Params(
        System,
        N = 4,
        temp_max = 100.0,
        temp_min = 1.0,
        accuracy = 1e-6,
    )

    sol, saved_values = PMFRG_general.SolveFRG(
        Par,
        anisotropy,
        SavedValues(
            PMFRG_general._getFloatType(Par),
            Observables_Full{PMFRG_general._getFloatType(Par)},
        ),
        (State, t, integrator) -> pack_all(State, t, Par),
    )



    save_object(
        "general_[-0.2 3.0 1.0]_polymer.jld2",
        [
            (saved_values.saveval[n], exp(saved_values.t[n])) for
            n = 1:length(saved_values.t)
        ],
    )
end


# Struct that contains all the relevant output 
# for the XYZ-specific implementation
struct Observables_Full_XYZ{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
    iSigma_x::Array{T,2}
    iSigma_y::Array{T,2}
    iSigma_z::Array{T,2}
    Gamma::Array{T,5}
end

function run_xyz()
    function pack_all_xyz(State, t, Par)
        chix = PMFRG_xyz.getChi_x(State, exp(t), Par)
        chiy = PMFRG_xyz.getChi_y(State, exp(t), Par)
        chiz = PMFRG_xyz.getChi_z(State, exp(t), Par)
        println("t: $(round(exp(t); digits=3))")
        return Observables_Full_XYZ(
            copy(chix),
            copy(chiy),
            copy(chiz),
            copy(State.x[2]),
            copy(State.x[3]),
            copy(State.x[4]),
            copy(State.x[5]),
        )
    end

    System = getPolymer(2)

    anisotropy = zeros(System.Npairs, 3)
    for Rij = 1:System.Npairs
        anisotropy[Rij, :] = [-0.2, 3.0, 1.0]
    end

    Par = PMFRG_xyz.Params(System, N = 4, temp_max = 100.0, temp_min = 1.0, accuracy = 1e-6)

    let
        sol, saved_values = PMFRG_xyz.SolveFRG(
            Par,
            anisotropy,
            SavedValues(
                PMFRG_xyz._getFloatType(Par),
                Observables_Full_XYZ{PMFRG_xyz._getFloatType(Par)},
            ),
            (State, t, integrator) -> pack_all_xyz(State, t, Par),
        )
        save_object(
            "xyz_[-0.2 3.0 1.0]_polymer.jld2",
            [
                (saved_values.saveval[n], exp(saved_values.t[n])) for
                n = 1:length(saved_values.t)
            ],
        )
    end

end



function compare_general_xyz()
    data_gen = load_object("general_[-0.2 3.0 1.0]_polymer.jld2")
    data_xyz = load_object("xyz_[-0.2 3.0 1.0]_polymer.jld2")

    sig_gen = getfield.(getindex.(data_gen, 1), :iSigma)
    sigx_xyz = getfield.(getindex.(data_xyz, 1), :iSigma_x)

    gam_gen = getfield.(getindex.(data_gen, 1), :Gamma)
    gam_xyz = getfield.(getindex.(data_xyz, 1), :Gamma)

    # Before switching sign in addTo1PartBubble (n=1):
    # [5.096792215915629e-9 1.7645491252554948e-9 1.1042344391208137e-9 7.171418595963982e-10]
    # [0.0, -1.3111058764277844e-8, 0.0, -1.2991740763648296e-8]
    # [0.0, -1.2210483735231037e-9, 0.0, -1.7117963713175754e-9]
    #
    # After switching sign in addTo1PartBubble (n=2):
    # [5.096792215915629e-9 1.7645491252554948e-9 1.1042344391208137e-9 7.171418595963982e-10]
    # [0.0, -1.3111058764277844e-8, 0.0, -1.2991740763648296e-8]
    # [0.0, -1.2210483735231037e-9, 0.0, -1.7117963713175754e-9]

    let
        n = 400
        println(sig_gen[n][1, :, :] .- sigx_xyz[n])

        # See debug.jl for comparing the whole struct 
        # for gamma
        println(gam_gen[n][1, 2, 1, 1, :] .- gam_xyz[n][1, 2, 1, 1, :]) # xxxx
        println(gam_gen[n][11, 2, 1, 1, :] .- gam_xyz[n][10, 2, 1, 1, :]) # xyxy
    end

end


main()
