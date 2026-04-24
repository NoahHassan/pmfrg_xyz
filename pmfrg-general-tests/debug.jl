using Pkg

include("../src/PMFRG_general.jl")
using .PMFRG_general

include("../src/PMFRG_old.jl")
using .PMFRG_old

function IndexFromFlavors(d1::Int, d2::Int, d3::Int, d4::Int)
    return d1 * 27 + d2 * 9 + d3 * 3 + d4 + 1
end

function FlavorsFromIndex(m::Int)
    n = m - 1
    d1 = div(n, 27)
    d2 = div((n - d1 * 27), 9)
    d3 = div((n - d1 * 27 - d2 * 9), 3)
    d4 = (n - d1 * 27 - d2 * 9 - d3 * 3)
    return (d1, d2, d3, d4)
end

function ConvertToXYZIndex(n::Int)
    if ((n - 1) % 40 == 0)
        return div(n - 1, 40) + 1
    end

    d1, d2, d3, d4 = FlavorsFromIndex(n)
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

function Generate_random_G(Nsites, N)
    G = zeros(81, Nsites, N, N, N)
    for p1 = 0:2, p2 = 0:2
        p1122 = IndexFromFlavors(p1, p1, p2, p2)
        p1212 = IndexFromFlavors(p1, p2, p1, p2)
        p1221 = IndexFromFlavors(p1, p2, p2, p1)

        G[p1122, :, :, :, :] .= rand(Nsites, N, N, N)
        G[p1212, :, :, :, :] .= rand(Nsites, N, N, N)
        G[p1221, :, :, :, :] .= rand(Nsites, N, N, N)
    end

    return G
end

function G_gen_to_xyz(Gamma::AbstractArray{T}, Nsites, N) where {T}
    # layout Gamma[flavor, site, s, t, u]
    G_xyz = zeros(21, Nsites, N, N, N)
    for d1 = 0:2, d2 = 0:2
        i1122_gen = IndexFromFlavors(d1, d1, d2, d2)
        i1122_xyz = ConvertToXYZIndex(i1122_gen)
        # println(i1122_gen, " ", i1122_xyz)
        i1212_gen = IndexFromFlavors(d1, d2, d1, d2)
        i1212_xyz = ConvertToXYZIndex(i1212_gen)
        # println(i1212_gen, " ", i1212_xyz)
        i1221_gen = IndexFromFlavors(d1, d2, d2, d1)
        i1221_xyz = ConvertToXYZIndex(i1221_gen)
        # println(i1221_gen, " ", i1221_xyz)

        G_xyz[i1122_xyz, :, :, :, :] = Gamma[i1122_gen, :, :, :, :]
        G_xyz[i1212_xyz, :, :, :, :] = Gamma[i1212_gen, :, :, :, :]
        G_xyz[i1221_xyz, :, :, :, :] = Gamma[i1221_gen, :, :, :, :]
    end

    return G_xyz
end

using SpinFRGLattices

function test_deriv()
    G_gen = Generate_random_G(2, 4)
    G_xyz = G_gen_to_xyz(G_gen, 2, 4)

    System = getPolymer(2)

    Par_gen = PMFRG_general.Params(
        System,
        N = 4,
        temp_max = 100.0,
        temp_min = 10.0,
        accuracy = 1e-4,
    )
    Par_xyz =
        PMFRG_old.Params(System, N = 4, temp_max = 100.0, temp_min = 10.0, accuracy = 1e-4)

    # Gamma_(index, ki, 0, -wmw1, -wpw1, flavTransform)
    PMFRG_general.TestFRGAt(Par_gen, zeros(System.Npairs, 3, 3), G_gen)
    PMFRG_old.TestFRGAt(Par_xyz, zeros(System.Npairs, 3), G_xyz)
end

test_deriv()


using JLD2

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

X_gen[5, :, :, :, :] .- X_xyz[ConvertToXYZIndex(5), :, :, :, :]

for n = 1:81
    index = ConvertToXYZIndex(n)
    if (index != -1)
        for m = 1:2, ns = 1:4, nt = 1:4, nu = 1:4
            if (X_gen[n, m, ns, nt, nu] - X_xyz[index, m, ns, nt, nu] != 0)
                println(X_gen[n, m, ns, nt, nu] - X_xyz[index, m, ns, nt, nu])
            end
        end
    end
end

A = [1, 2, 3, 4, 5, 6, 7, 8, 9]
reshape(A, 3, 3)'
