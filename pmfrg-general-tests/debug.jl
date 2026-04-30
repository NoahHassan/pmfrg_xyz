using Pkg
using Test

include("../src/PMFRG_general.jl")
using .PMFRG_general

include("../src/PMFRG_old.jl")
using .PMFRG_old

_N = 4

function GenIndexFromFlavors(d1::Int, d2::Int, d3::Int, d4::Int)
    return d1 * 27 + d2 * 9 + d3 * 3 + d4 + 1
end

function FlavorsFromGenIndex(m::Int)
    n = m - 1
    d1 = div(n, 27)
    d2 = div((n - d1 * 27), 9)
    d3 = div((n - d1 * 27 - d2 * 9), 3)
    d4 = (n - d1 * 27 - d2 * 9 - d3 * 3)
    return (d1, d2, d3, d4)
end

function XYZIndex_from_GenIndex(n::Int)
    if ((n - 1) % 40 == 0)
        return div(n - 1, 40) + 1
    end

    d1, d2, d3, d4 = FlavorsFromGenIndex(n)
    if (d1 == d2)
        if (d3 != d4)
            return -1
        end
        if (d1 == 2 && (d3 + 2) % 3 == 0)
            return 9
        else
            return 3 + d1 + (d3 + 2) % 3 + 1
        end
    end
    if (d1 == d3)
        if (d2 != d4)
            return -1
        end
        if (d1 == 2 && (d2 + 2) % 3 == 0)
            return 15
        else
            return 9 + d1 + (d2 + 2) % 3 + 1
        end
    end
    if (d1 == d4)
        if (d2 != d3)
            return -1
        end
        if (d1 == 2 && (d2 + 2) % 3 == 0)
            return 21
        else
            return 15 + d1 + (d2 + 2) % 3 + 1
        end
    end

    return -1
end

function random_XYZ_Gamma_general(Nsites, N)
    G = zeros(81, Nsites, N, N, N)
    for p1 = 0:2, p2 = 0:2
        if p1 != p2
            let p1122 = GenIndexFromFlavors(p1, p1, p2, p2) # 4-9
                G[p1122, :, :, :, :] .= rand(Nsites, N, N, N)
            end

            let p1212 = GenIndexFromFlavors(p1, p2, p1, p2) # 10-15
                G[p1212, :, :, :, :] .= rand(Nsites, N, N, N)
            end

            let p1221 = GenIndexFromFlavors(p1, p2, p2, p1) # 16-21
                G[p1221, :, :, :, :] .= rand(Nsites, N, N, N)
            end
        else
            let p1111 = GenIndexFromFlavors(p1, p1, p1, p1)
                G[p1111, :, :, :, :] .= rand(Nsites, N, N, N)
            end
        end
    end

    return G
end

function Gamma_general_to_XYZ(Gamma::AbstractArray{T}, Nsites, N) where {T}
    # layout Gamma[flavor, site, s, t, u]
    G_xyz = zeros(21, Nsites, N, N, N)
    for d1 = 0:2, d2 = 0:2
        if d1 != d2
            let i1122_gen = GenIndexFromFlavors(d1, d1, d2, d2),
                i1122_xyz = XYZIndex_from_GenIndex(i1122_gen)

                G_xyz[i1122_xyz, :, :, :, :] = Gamma[i1122_gen, :, :, :, :]
            end

            let i1212_gen = GenIndexFromFlavors(d1, d2, d1, d2),
                i1212_xyz = XYZIndex_from_GenIndex(i1212_gen)

                G_xyz[i1212_xyz, :, :, :, :] = Gamma[i1212_gen, :, :, :, :]
            end

            let i1221_gen = GenIndexFromFlavors(d1, d2, d2, d1),
                i1221_xyz = XYZIndex_from_GenIndex(i1221_gen)

                G_xyz[i1221_xyz, :, :, :, :] = Gamma[i1221_gen, :, :, :, :]
            end
        else
            let i1111_gen = GenIndexFromFlavors(d1, d1, d1, d1),
                i1111_xyz = XYZIndex_from_GenIndex(i1111_gen)

                G_xyz[i1111_xyz, :, :, :, :] = Gamma[i1111_gen, :, :, :, :]
            end
        end
    end

    return G_xyz
end

using SpinFRGLattices

function test_deriv()
    global _N
    G_gen = random_XYZ_Gamma_general(2, _N)
    G_xyz = Gamma_general_to_XYZ(G_gen, 2, _N)

    System = getPolymer(2)

    Par_gen = PMFRG_general.Params(
        System,
        N = _N,
        temp_max = 1.0,
        temp_min = 10.0,
        accuracy = 1e-4,
    )
    Par_xyz =
        PMFRG_old.Params(System, N = 4, temp_max = 1.0, temp_min = 10.0, accuracy = 1e-4)

    # Gamma_(index, ki, 0, -wmw1, -wpw1, flavTransform)
    PMFRG_general.TestFRGAt(Par_gen, zeros(System.Npairs, 3, 3), G_gen)
    PMFRG_old.TestFRGAt(Par_xyz, zeros(System.Npairs, 3), G_xyz)
end



using JLD2

function comparison(x, y, rel_threshold)
    if x == 0 && y == 0
        return false
    else
        abs(x - y) / (abs(x) + abs(y)) > rel_threshold
    end
end


function main()

    test_deriv()
    dg_gen = load_object("Dgamma_gen.jld2")
    dg_xyz = load_object("Dgamma_xyz.jld2")

    PMFRG_general.iSigma_(dg_gen, 1, 1)
    PMFRG_old.iSigma_(dg_xyz.x, 1, 0)

    dg_gen[:, 1, 1]
    dg_xyz.x

    sk_gen = load_object("SKat_gen.jld2")
    sk_xyz = load_object("SKat_xyz.jld2")

    sk_gen[3]
    sk_xyz[3]

    g_gen = load_object("G_gen.jld2")
    g_xyz = load_object("G_xyz.jld2")

    g_gen[1]
    g_xyz[1]

    -g_xyz[1]^2 * dg_xyz.x[1]

    M = g_gen[1]
    DG = reshape(dg_gen[:, 1, 1], 3, 3)'
    -M * M * DG

    # doing sign(P) in the addX-Calculation gives sufficient aggreement
    # try disabling Katanin
    X_gen = load_object("X_gen_test.jld2")
    X_xyz = load_object("X_xyz_test.jld2")

    X_gen[5, :, :, :, :] .- X_xyz[XYZIndex_from_GenIndex(5), :, :, :, :]

    X_xyz_from_general = similar(X_xyz)
    X_xyz_from_general .= 0.0

    for n = 1:81
        index = XYZIndex_from_GenIndex(n)
        if (index != -1)
            for m = 1:2, ns = 1:4, nt = 1:4, nu = 1:4
                X_xyz_from_general[index, m, ns, nt, nu] = X_gen[n, m, ns, nt, nu]
                X_xyz_from_general[index+21, m, ns, nt, nu] = X_gen[n+81, m, ns, nt, nu]


                if comparison(X_gen[n, m, ns, nt, nu], X_xyz[index, m, ns, nt, nu], 1.0e-9)
                    println(
                        "diff (X): ",
                        X_gen[n, m, ns, nt, nu] - X_xyz[index, m, ns, nt, nu],
                    )
                end

                if comparison(
                    X_gen[n+81, m, ns, nt, nu],
                    X_xyz[index+21, m, ns, nt, nu],
                    1.0e-9,
                )
                    println(
                        "diff (Y): ",
                        X_gen[n+81, m, ns, nt, nu] - X_xyz[index+21, m, ns, nt, nu],
                    )
                end

            end
        end
    end


    @testset verbose = true "global" begin
        @testset verbose = true "X/Y equality" begin
            @testset "global equivalence" begin
                @test X_xyz_from_general ≈ X_xyz
            end

            @testset "X part only" begin
                @test X_xyz_from_general[1:21, :, :, :, :] ≈ X_xyz[1:21, :, :, :, :]
            end
            @testset "Y part only" begin
                @test X_xyz_from_general[22:42, :, :, :, :] ≈ X_xyz[22:42, :, :, :, :]
            end
        end



        function compare_spropY(N)
            function map_xyz_gen(i_xyz)
                d = Dict(1 => 1, 2 => 5, 3 => 9)
                d[i_xyz]
            end

            @testset verbose = true "spropY equality" begin
                for ns = 1:N, nt = 1:N, nw = -N:N-1
                    gen = load_object("spropY/gen-$ns-$nt-$nw.jld2")
                    xyz = load_object("spropY/xyz-$ns-$nt-$nw.jld2")
                    for d1 = 1:3, d2 = 1:3
                        @test xyz[:, :, d1, d2] ≈
                              gen[:, :, map_xyz_gen(d1), map_xyz_gen(d2)]
                    end
                    nonzero_indices = (1, 5, 9)

                    for d1 = 1:9, d2 = 1:9
                        if !((d1 in nonzero_indices) && (d2 in nonzero_indices))
                            @test gen[1, 1, d1, d2] == 0.0
                        end
                    end
                end
            end

        end


        function compare_Y_detailed(X_xyz, X_xyz_from_general)
            Y = X_xyz[22:end, :, :, :, :]
            Y_from_general = X_xyz_from_general[22:end, :, :, :, :]

            Y_zero_counts = 0
            Y_not_found = 0
            for n = 1:21, m = 1:2, ns = 1:4, nt = 1:4, nu = 1:4
                diff = +Inf
                map = nothing
                # Search for the minimum difference,
                # hoping thus to find the corresponding value
                # (assuming it's a "mapping" problem)
                for n_ = 1:21, m_ = 1:2, ns_ = 1:4, nt_ = 1:4, nu_ = 1:4
                    d = abs(Y[n, m, ns, nt, nu] - Y_from_general[n_, m_, ns_, nt_, nu_])
                    if d < diff
                        map = (n_, m_, ns_, nt_, nu_)
                        diff = d
                    end
                end
                v = abs(Y[n, m, ns, nt, nu])
                n_, m_, ns_, nt_, nu_ = map
                if Y[n, m, ns, nt, nu] == 0 && Y_from_general[n_, m_, ns_, nt_, nu_] == 0
                    Y_zero_counts += 1
                else
                    if (diff / v < 1.0e-6) && ((n, m, ns, nt, nu) != map)
                        print("$n $m $ns $nt $nu => $n_ $m_ $ns_ $nt_ $nu_  ")
                        print("$(Y[n,m,ns,nt,nu]) ≈ $(Y_from_general[n_,m_,ns_,nt_,nu_])")
                        if m != m_ || ns != ns_ || nt != nt_ || nu != nu_
                            print(" Special case")
                        end
                        print("\n")
                    else
                        # println("$n $m $ns $nt $nu => not found")
                        Y_not_found += 1
                    end
                end

            end

            println("zero counts: $Y_zero_counts")
            println("not found: $Y_not_found")
        end

        compare_Y_detailed(X_xyz, X_xyz_from_general)
        compare_spropY(_N)
    end
end

if PROGRAM_FILE == @__FILE__
    main()
end
