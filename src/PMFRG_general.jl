module PMFRG_general

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

using RecursiveArrayTools
using SpinFRGLattices, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools, StructArrays
using SpinFRGLattices.StaticArrays
using LinearAlgebra

setZero!(a::AbstractArray{T,N}) where {T,N} = fill!(a, zero(T))

function setZero!(PartArr::ArrayPartition)
    for arr in PartArr.x
        fill!(arr, 0.0)
    end
end

"""Recursively sets structure to zero"""
function setZero!(a::T) where {T}
    for f in fieldnames(T)
        setZero!(getfield(a, f))
    end
    return a
end

struct StateType{T}
    f_int::Vector{T}
    iSigma::Array{T,3} # 9 Sigma-flavors
    Gamma::Array{T,5} # 81 Gamma-flavors
end

struct Observables{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
end

struct NumericalParams{T<:Real}
    N::Int

    accuracy::T
    temp_min::T
    temp_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

struct OptionParams
    use_symmetry::Bool
    minimal_output::Bool
end

struct OneLoopParams_1{T,SType}
    System::SType
    NumericalParams::NumericalParams{T}
    Options::OptionParams
end

struct OneLoopWorkspace{T,ParType}
    State::StateType{T}
    Deriv::StateType{T}
    X::Array{T,5} # 81 + 81 = 162 flavor combinations
    Par::ParType
end

getVDims(Par) = (
    81,
    Par.System.Npairs,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
)
getBubbleVDims(Par) = (
    162,
    Par.System.Npairs,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
)
_getFloatType(Par) = typeof(Par.NumericalParams.accuracy)

function StateType(NUnique::Int, N::Int, VDims::Tuple, type = Float64)
    return StateType(zeros(type, NUnique), zeros(type, 9, NUnique, N), zeros(type, VDims))
end

StateType(Par) =
    StateType(Par.System.NUnique, Par.NumericalParams.N, getVDims(Par), _getFloatType(Par))
RecursiveArrayTools.ArrayPartition(x) = ArrayPartition(x.f_int, x.iSigma, x.Gamma)
StateType(Arr::ArrayPartition) = StateType(Arr.x...)

function NumericalParams(;
    N::Integer = 24,
    accuracy = 1e-6,
    temp_min = exp(-10.0),
    temp_max = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, lenIntw)),
)

    return NumericalParams(N, accuracy, temp_min, temp_max, lenIntw, lenIntw_acc)
end

function OneLoopWorkspace(State, Deriv, X, Par)
    setZero!(Deriv)
    setZero!(X)

    return OneLoopWorkspace(StateType(State.x...), StateType(Deriv.x...), X, Par)
end

OptionParams(; use_symmetry::Bool = true, MinimalOutput::Bool = false, kwargs...) =
    OptionParams(use_symmetry, MinimalOutput)
Params(System; kwargs...) =
    OneLoopParams_1(System, NumericalParams(; kwargs...), OptionParams(; kwargs...))

#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

function get_w(nw, T)
    return pi * (2 * nw + 1)
end

function get_sign_iw(nw::Integer, N::Integer)
    nw_bounds = min(nw, N - 1)
    return nw_bounds + 1
end

### The input used to be one of iSigma_x, iSigma_y, iSigma_z
### Now it will be something like iSigma[flavorIndex, :, :]
function iSigma_(iSigma::AbstractArray, x::Integer, nw::Integer)
    N = size(iSigma, 3)
    s = 1
    if nw < 0
        nw = -nw - 1
        s = -1
    end
    iw = get_sign_iw(nw, N)
    return s * iSigma[:, x, iw]
end

# the inverse is given by ∼ (δ_αβ ω + iΣ_αβ)

#iSigma inputted as 1d-array
function Calculate_G_Inverse(w::Real, iSigma::AbstractArray)
    return [
        (w+iSigma[1]) iSigma[2] iSigma[3]
        iSigma[4] (w+iSigma[5]) iSigma[6]
        iSigma[7] iSigma[8] (w+iSigma[9])
    ]
    # w*I + iSigma
end

function iG_(iSigma::AbstractArray, x::Integer, nw::Integer, T::Real)
    w = get_w(nw, T)
    return inv(Calculate_G_Inverse(w * sqrt(T), iSigma_(iSigma, x, nw)))
end

### by differentiating the above inverse by T
function iS_(iSigma::AbstractArray, x::Integer, nw::Integer, T::Real)
    iG = iG_(iSigma, x, nw, T)
    w = get_w(nw, T)
    return -iG .* (w / (2.0 * sqrt(T))) .* iG
end

function iSKat_(
    iSigma::AbstractArray,
    DSigma::AbstractArray,
    x::Integer,
    nw::Integer,
    T::Real,
)
    iG = iG_(iSigma, x, nw, T)
    w = get_w(nw, T)

    # return iS_(iSigma, x, nw, T)

    # Attention here
    return -iG * (reshape(iSigma_(DSigma, x, nw), 3, 3)' + (w / (2.0 * sqrt(T))) * I) * iG
end

####################################################
######### VERTICES ## VERTICES ## VERTICES #########
####################################################

# The XYZ model used a transformation to map flavors into each other
# after sign change in one of the frequencies. We have to find a
# different solution now. There are three 1-element equivalence classes,
# 9 2-element classes and 15 4-element classes. The 1- and 2-element classes
# Transform together under XYZ symmetry. The additional 15 off-diagonal
# classes transform differnetly.
#
# A class is defined as a set of flavor-combinations that transform into
# each other upon change of frequency signs.
#
# 1-element classes:
# xxxx (1 to 3)
# These transform as identity under Ms, Mt and Mu
#
# 2-element classes:
# xxyy (4-9)
# xyxy (10-15)
# xyyx (16-21)
# These have only one identity-transformation (Ms, Mt, Mu respectively)
#
# 4-element classes:
# xxxy (22-45)
# xxyz (46-69)
# xyxz (70-81)
# These transform under all Ms, Mt and Mu transformations

# m indexed starting form 1
function FlavorsFromIndex(m::Int)
    n = m - 1
    d1 = div(n, 27)
    d2 = div((n - d1 * 27), 9)
    d3 = div((n - d1 * 27 - d2 * 9), 3)
    d4 = (n - d1 * 27 - d2 * 9 - d3 * 3)
    return (d1, d2, d3, d4)
end

# d1, d2, d3, d4 indexed from zero
function IndexFromFlavors(d1::Int, d2::Int, d3::Int, d4::Int)
    return d1 * 27 + d2 * 9 + d3 * 3 + d4 + 1
end

function Klein4_Permute(n::Int, ns::Int, nt::Int, nu::Int)
    d1, d2, d3, d4 = FlavorsFromIndex(n)

    function Permute_S(d1, d2, d3, d4)
        return (d2, d1, d4, d3)
    end
    function Permute_T(d1, d2, d3, d4)
        return (d3, d4, d1, d2)
    end
    function Permute_U(d1, d2, d3, d4)
        return (d4, d3, d2, d1)
    end

    if (ns < 0)
        d1, d2, d3, d4 = Permute_S(d1, d2, d3, d4)
    end
    if (nt < 0)
        d1, d2, d3, d4 = Permute_T(d1, d2, d3, d4)
    end
    if (nu < 0)
        d1, d2, d3, d4 = Permute_U(d1, d2, d3, d4)
    end

    return IndexFromFlavors(d1, d2, d3, d4)
end

function ConvertFreqArgs(ns, nt, nu, Nw)
    ns, nt, nu = abs.((ns, nt, nu))

    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2) ### weird cutoff, idk why
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)

    return ns, nt, nu
end

using LinearAlgebra
using SparseArrays

# I include isFlavorTransform for optimization purposes. The integer n
# labels the Vertex flavor.
function V_(
    Vertex::AbstractArray,
    n::Int,
    ns::Int,
    nt::Int,
    nu::Int,
    Rij::Integer,
    Rji::Integer,
    N::Integer,
)

    n_transf = Klein4_Permute(n, ns, nt, nu)

    ns, nt, nu = ConvertFreqArgs(ns, nt, nu, N)
    Rij = ifelse(nt * nu < 0, Rji, Rij)
    return Vertex[n_transf, Rij, ns+1, nt+1, nu+1]
end

function mixedFrequencies(ns, nt, nu, nwpr)
    nw1 = Int((ns + nt + nu - 1) / 2)
    nw2 = Int((ns - nt - nu - 1) / 2)
    nw3 = Int((-ns + nt - nu - 1) / 2)
    nw4 = Int((-ns - nt + nu - 1) / 2)

    wpw1 = nwpr + nw1 + 1
    wpw2 = nwpr + nw2 + 1
    wpw3 = nwpr + nw3 + 1
    wpw4 = nwpr + nw4 + 1
    wmw1 = nwpr - nw1
    wmw2 = nwpr - nw2
    wmw3 = nwpr - nw3
    wmw4 = nwpr - nw4

    return wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4
end

function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
    (; State, X, Par) = Workspace
    N = Par.NumericalParams.N
    (; Npairs, Nsum, siteSum, invpairs) = Par.System

    Vert(n, Rij, s, t, u) = V_(State.Gamma, n, s, t, u, Rij, invpairs[Rij], N)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    # get fields of siteSum struct as Matrices for better use of LoopVectorization
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    X_sum = @MVector zeros(81)
    for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        fill!(X_sum, 0.0)
        sumsum = 0
        for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m, xk =
                S_ki[k_spl, Rij], S_kj[k_spl, Rij], S_m[k_spl, Rij], S_xk[k_spl, Rij]
            Ptm(n_a, n_b) = Props[xk, xk, n_a, n_b] * m

            V12 = @SVector [Vert(n, ki, ns, wpw1, -wpw2) for n = 1:81]
            V34 = @SVector [Vert(n, kj, ns, -wmw3, -wmw4) for n = 1:81]

            for _n = 1:81
                for _m = 1:81
                    d1, d2, d3, d4 = FlavorsFromIndex(_n)
                    g1, g2, g3, g4 = FlavorsFromIndex(_m)
                    v12_flavor = IndexFromFlavors(d1, d2, g1, g4)
                    v34_flavor = IndexFromFlavors(g2, g3, d3, d4)

                    X_sum[_n] +=
                        -V12[v12_flavor] *
                        V34[v34_flavor] *
                        Ptm(g1 * 3 + g2 + 1, g3 * 3 + g4 + 1)
                end
            end

        end
        X[1:81, Rij, is, it, iu] .+= X_sum
    end
    return
end

function addY!(
    Workspace,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Props;
    _l = 1.0,
)
    (; State, X, Par) = Workspace
    N = Par.NumericalParams.N
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    Vert(n, Rij, s, t, u) = V_(State.Gamma, n, s, t, u, Rij, invpairs[Rij], N)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    X_sum = @MVector zeros(81)

    # Xtilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        # loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij) 
        (; xi, xj) = PairTypes[Rij]

        # For some reason promoting P_ and PT_ to SMatrix objects
        # reduces performance slightly
        function P_(n::Int, m::Int)
            return Props[xi, xj, n, m]
        end

        function PT_(n::Int, m::Int)
            return Props[xj, xi, m, n]
        end

        V13 = @SVector [Vert(n, Rij, -wmw1, nt, wmw3) for n = 1:81]
        V24 = @SVector [Vert(n, Rij, wpw2, -nt, -wpw4) for n = 1:81]
        V31 = @SVector [Vert(n, Rij, wmw3, nt, -wmw1) for n = 1:81]
        V42 = @SVector [Vert(n, Rij, -wpw4, -nt, wpw2) for n = 1:81]

        fill!(X_sum, 0.0)

        for n = 1:81
            for m = 1:81
                d1, d2, d3, d4 = FlavorsFromIndex(n)
                g1, g2, g3, g4 = FlavorsFromIndex(m)
                v13_flavor = IndexFromFlavors(d1, g2, d3, g3)
                v24_flavor = IndexFromFlavors(d2, g1, d4, g4)
                v31_flavor = IndexFromFlavors(d1, g3, d3, g2)
                v42_flavor = IndexFromFlavors(d2, g4, d4, g1)

                X_sum[n] += (
                    V13[v13_flavor] *
                    V24[v24_flavor] *
                    P_(g1 * 3 + g2 + 1, g3 * 3 + g4 + 1) +
                    V31[v31_flavor] *
                    V42[v42_flavor] *
                    PT_(g1 * 3 + g2 + 1, g3 * 3 + g4 + 1)
                )
            end
        end

        X[82:end, Rij, is, it, iu] .+= X_sum
    end
end

function getXBubble!(Workspace, T::Real)
    Par = Workspace.Par
    (; N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System

    iG_ab(x, nw) = iG_(Workspace.State.iSigma, x, nw, T)
    iSKat_ab(x, nw) = iSKat_(Workspace.State.iSigma, Workspace.Deriv.iSigma, x, nw, T)

    function getKataninProp!(BubbleProp, nw1, nw2)

        for i = 1:Par.System.NUnique, j = 1:Par.System.NUnique
            for n_a = 1:9, n_b = 1:9
                x_a = div(n_a - 1, 3) + 1
                y_a = (n_a - 1) % 3 + 1
                x_b = div(n_b - 1, 3) + 1
                y_b = (n_b - 1) % 3 + 1
                BubbleProp[i, j, n_a, n_b] =
                    iSKat_ab(i, nw1)[x_a, y_a] * iG_ab(j, nw2)[x_b, y_b]
            end
        end

        return -BubbleProp
    end

    for is = 1:N, it = 1:N
        BubbleProp = zeros(NUnique, NUnique, 9, 9)
        ns = is - 1
        nt = it - 1
        for nw = -lenIntw:lenIntw-1 # Matsubara sum
            spropX = getKataninProp!(BubbleProp, nw, nw + ns)
            spropY = getKataninProp!(BubbleProp, nw, nw - nt)
            for iu = 1:N
                nu = iu - 1
                if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                    continue
                end
                addY!(Workspace, is, it, iu, nw, spropY) # add to XTilde-type bubble functions

                ### If no u--t symmetry, then add all the bubbles
                ### If use u--t symmetry, then only add for nu smaller then nt (all other obtained by symmetry)
                # if(!Par.Options.use_symmetry || nu<=nt)

                addX!(Workspace, is, it, iu, nw, spropX)
                # end
            end
        end
    end
end

function symmetrizeBubble!(X::Array{T,5}, Par) where {T}
    N = Par.NumericalParams.N
    (; Npairs, OnsitePairs) = Par.System
    use_symmetry = Par.Options.use_symmetry
    # use the u <--> t symmetry
    if (use_symmetry)
        # for it in 1:N
        #     for iu in it+1:N, is in 1:N, Rij in 1:Npairs
        #         X.a[Rij,is,it,iu] = -X.a[Rij,is,iu,it]
        #         X.b[Rij,is,it,iu] = -X.b[Rij,is,iu,it]
        #         X.c[Rij,is,it,iu] = (
        #         + X.a[Rij,is,it,iu]+
        #         - X.b[Rij,is,it,iu]+
        #         + X.c[Rij,is,iu,it])
        #     end
        # end
    end
    #local definitions of X.Tilde vertices
    for iu = 1:N
        for it = 1:N, is = 1:N, R in OnsitePairs
            for n = 1:81
                d1, d2, d3, d4 = FlavorsFromIndex(n)
                n_transf = IndexFromFlavors(d1, d3, d2, d4)
                X[81+n, R, is, it, iu] = -X[n_transf, R, it, is, iu]
            end
        end
    end
end

function addToVertexFromBubble!(Gamma::Array{T,5}, X::Array{T,5}) where {T}
    for iu in axes(Gamma, 5)
        for it in axes(Gamma, 4), is in axes(Gamma, 3), Rij in axes(Gamma, 2)
            for n = 1:81
                d1, d2, d3, d4 = FlavorsFromIndex(n)
                n_transf = IndexFromFlavors(d1, d2, d4, d3)

                Gamma[n, Rij, is, it, iu] += (
                    X[n, Rij, is, it, iu] + X[81+n, Rij, is, it, iu] -
                    X[81+n_transf, Rij, is, iu, it]
                )
            end
        end
    end
    return Gamma
end

function symmetrizeVertex!(Gamma::Array{T,5}, Par) where {T}
    N = Par.NumericalParams.N
    for iu = 1:N
        for it = 1:N, is = 1:N, R in Par.System.OnsitePairs
            for n = 1:81
                d1, d2, d3, d4 = FlavorsFromIndex(n)
                n_transf = IndexFromFlavors(d1, d3, d2, d4)
                Gamma[n, R, is, it, iu] = -Gamma[n_transf, R, it, is, iu]
            end
        end
    end
end

######################################################################
######### FLOW EQUATIONS ## FLOW EQUATIONS ## FLOW EQUATIONS #########
######################################################################

# function getDFint!(Workspace, T::Real)
#     (; State, Deriv, Par) = Workspace
#     (; lenIntw_acc) = Par.NumericalParams
#     NUnique = Par.System.NUnique

# 	iSigma_ab(n, x, nw) = iSigma_(State.iSigma[n, :, :], x, nw)
# 	iG_ab(n, x, nw) = iG_(State.iSigma[n, :, :], x, nw, T)
# 	iS_ab(n, x, nw) = iS_(State.iSigma[n, :, :], x, nw, T)

# 	for x in 1:NUnique
# 		sumres = 0.
# 		for nw in -lenIntw_acc:lenIntw_acc-1
# 			w = get_w(nw,T)
# 			sumres += iSx(x, nw) / iGy(x, nw) * iSigmax(x, nw) / w
#             sumres += iSy(x, nw) / iGy(x, nw) * iSigmay(x, nw) / w
#             sumres += iSz(x, nw) / iGz(x, nw) * iSigmaz(x, nw) / w
#         end
# 		Deriv.f_int[x] = -0.5 * sumres
# 	end
# end

function get_Self_Energy!(Workspace, T::Real)
    println(T)
    Par = Workspace.Par
    @inline iS_ab(x, nw) = vec(iS_(Workspace.State.iSigma, x, nw, T)') ./ 2
    compute1PartBubble!(Workspace.Deriv.iSigma, Workspace.State.Gamma, iS_ab, Par)
end

function compute1PartBubble!(Dgamma::Array{T,3}, Gamma::Array{T,5}, Props, Par) where {T}
    invpairs = Par.System.invpairs

    setZero!(Dgamma)
    @inline Gamma_(n, Rij, s, t, u) =
        V_(Gamma, n, s, t, u, Rij, invpairs[Rij], Par.NumericalParams.N)
    addTo1PartBubble!(Dgamma, Gamma_, Props, Par)
end

function IndexToString(n::Int)
    d1, d2, d3, d4 = FlavorsFromIndex(n)
    letters = ["x", "y", "z"]
    return letters[d1+1] * letters[d2+1] * letters[d3+1] * letters[d4+1]
end

function addTo1PartBubble!(Dgamma::Array{T,3}, Gamma_::Function, Props, Par) where {T}

    (; N, lenIntw_acc) = Par.NumericalParams
    (; siteSum, Nsum, OnsitePairs) = Par.System

    Threads.@threads for iw1 = 1:N
        nw1 = iw1 - 1
        for (x, Rx) in enumerate(OnsitePairs)
            for nw = -lenIntw_acc:lenIntw_acc-1
                jsum = zeros(9)
                wpw1 = nw1 + nw + 1
                wmw1 = nw - nw1
                for k_spl = 1:Nsum[Rx]
                    (; m, ki, xk) = siteSum[k_spl, Rx]
                    gam = @SVector [Gamma_(n, ki, 0, -wmw1, -wpw1) for n = 1:81]
                    for d1 = 0:2, d2 = 0:2
                        sig_index = d1 * 3 + d2 + 1
                        for d3 = 0:2, d4 = 0:2
                            prop_index = d3 * 3 + d4 + 1
                            gam_index = IndexFromFlavors(d3, d4, d1, d2)
                            jsum[sig_index] +=
                                Props(xk, nw)[prop_index] * gam[gam_index] * m
                        end
                    end
                end
                Dgamma[:, x, iw1] .+= -jsum
            end
        end
    end
end

using JLD2
function getDeriv!(Deriv, State, setup, Lam; saveArgs = true)

    (X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(State, Deriv, X, Par)

    # getDFint!(Workspace, Lam)
    get_Self_Energy!(Workspace, Lam)
    getXBubble!(Workspace, Lam)
    symmetrizeBubble!(Workspace.X, Par)
    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.X)
    symmetrizeVertex!(Workspace.Deriv.Gamma, Par)

    return
end

####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################

t_to_Lam(t) = exp(t)
Lam_to_t(t) = log(t)

function AllocateSetup(Par::OneLoopParams_1)
    println("Allocate Setup")
    ## Allocate Memory:
    floattype = _getFloatType(Par)
    X = zeros(floattype, getBubbleVDims(Par))
    return (X, Par)
end

function InitializeState(Par, isotropy)

    N = Par.NumericalParams.N
    (; couplings, NUnique) = Par.System

    VDims = getVDims(Par)
    #floattype = _getFloatType(Par)

    State = ArrayPartition(
        zeros(NUnique),          ### f_int
        zeros(9, NUnique, N),    ### iSigma
        zeros(VDims),            ### Gamma
    )

    Gamma = State.x[3]

    setToBareVertex!(Gamma, couplings, isotropy)

    println(getChi_x(State, 100.0, Par))
    println(getChi_y(State, 100.0, Par))
    println(getChi_z(State, 100.0, Par))

    return State

end

function InitializeAt(Par, Gamma)
    N = Par.NumericalParams.N
    (; couplings, NUnique) = Par.System

    VDims = getVDims(Par)
    #floattype = _getFloatType(Par)

    State = ArrayPartition(
        zeros(NUnique),          ### f_int
        zeros(9, NUnique, N),    ### iSigma
        copy(Gamma),             ### Gamma
    )

    Gamma = State.x[3]

    println(getChi_x(State, 100.0, Par))
    println(getChi_y(State, 100.0, Par))
    println(getChi_z(State, 100.0, Par))

    return State
end

function save_static_chis(State, t, Par)
    chi_x = getChi_x(State, t_to_Lam(t), Par)
    chi_y = getChi_y(State, t_to_Lam(t), Par)
    chi_z = getChi_z(State, t_to_Lam(t), Par)

    println("t: $(round(exp(t);digits=3))")

    return Observables(copy(chi_x), copy(chi_y), copy(chi_z))
end

function gettMesh(T_min, T_max, npoints)
    t_min = get_t_min(T_min)
    t_max = Lam_to_t(T_max)
    return LinRange(t_min, t_max, npoints)
end

function launchPMFRG!(
    State,
    setup,
    Deriv!::Function,
    saved_values::SavedValues,
    save_func::Function;
    method = DP5(),
    npoints = 600,
)

    Par = setup[end]
    (; temp_max, temp_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(temp_max)
    tend = get_t_min(temp_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    ObsSaveat = gettMesh(temp_min, temp_max, npoints)
    saveCB = SavingCallback(
        save_func,
        saved_values,
        save_everystep = false,
        saveat = ObsSaveat,
        tdir = -1,
    )

    problem = ODEProblem(Deriv_subst!, State, (t0, tend), setup) # function, initial state, timespan, ??
    sol = solve(
        problem,
        method,
        reltol = accuracy,
        abstol = accuracy,
        save_everystep = false,
        callback = saveCB,
        dt = Lam_to_t(0.2 * temp_max),
    )

    return sol, saved_values
end

function testPMFRG!(State, setup, Deriv!::Function)
    Par = setup[end]
    (; temp_max, temp_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(temp_max)
    tend = get_t_min(temp_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    der = copy(State)
    setZero!(der)

    Deriv_subst!(der, State, setup, t0)
end

SolveFRG(Par, isotropy; kwargs...) = launchPMFRG!(
    InitializeState(Par, isotropy),
    AllocateSetup(Par),
    getDeriv!,
    SavedValues(_getFloatType(Par), Observables{_getFloatType(Par)}),
    (State, t, _) -> save_static_chis(State, t, Par);
    kwargs...,
)
SolveFRG(Par, isotropy, saved_values, save_func; kwargs...) = launchPMFRG!(
    InitializeState(Par, isotropy),
    AllocateSetup(Par),
    getDeriv!,
    saved_values,
    save_func;
    kwargs...,
)

TestFRG(Par, isotropy; kwargs...) =
    testPMFRG!(InitializeState(Par, isotropy), AllocateSetup(Par), getDeriv!; kwargs...)
TestFRGAt(Par, isotropy, Gamma; kwargs...) =
    testPMFRG!(InitializeAt(Par, Gamma), AllocateSetup(Par), getDeriv!; kwargs...)

function get_t_min(Lam)
    Lam < exp(-30) && @warn "temp_min too small! Set to exp(-30) instead."
    max(Lam_to_t(Lam), -30.0)
end

function generateSubstituteDeriv(getDeriv!::Function)

    function DerivSubs!(Deriv, State, par, t; s = true)
        Lam = t_to_Lam(t)
        a = getDeriv!(Deriv, State, par, Lam, saveArgs = s)
        Deriv .*= Lam
        a
    end

end

function MapIndexToWord(n::Int)
    letters = ['x', 'y', 'z']
    d1, d2, d3, d4 = FlavorsFromIndex(n)
    return string(letters[d1+1], letters[d2+1], letters[d3+1], letters[d4+1])
end

function setToBareVertex!(
    Gamma::AbstractArray{T,5},
    couplings::AbstractVector,
    anisotropy::Array{T,3},
) where {T}
    epsilon_mat = reshape(
        [
            0 0 0
            0 0 1
            0 -1 0
            0 0 -1
            0 0 0
            1 0 0
            0 1 0
            -1 0 0
            0 0 0
        ],
        3,
        3,
        3,
    )

    epsilon(i, j, k) = epsilon_mat[k, j, i]

    for Rj in axes(Gamma, 2)
        for n_flavor = 1:81
            d1, d2, d3, d4 = FlavorsFromIndex(n_flavor)
            d1 += 1
            d2 += 1
            d3 += 1
            d4 += 1

            for g1 = 1:3, g2 = 1:3
                Gamma[n_flavor, Rj, :, :, :] .+=
                    -1.0 *
                    epsilon(d1, d2, g1) *
                    couplings[Rj] *
                    anisotropy[Rj, g1, g2] *
                    epsilon(d3, d4, g2)
            end
        end
    end

    for n = 1:81
        if (abs(Gamma[n, 2, 1, 1, 1]) > 0.0)
            word = MapIndexToWord(n)
            println("Γ_$(word) = $(Gamma[n, 2, 1, 1, 1])")
        end
    end

    return Gamma
end

#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################

getChi_z(State::ArrayPartition, T::Real, Par) = getChi_z(State.x[2], State.x[3], T, Par)
getChi_x(State::ArrayPartition, T::Real, Par) = getChi_x(State.x[2], State.x[3], T, Par)
getChi_y(State::ArrayPartition, T::Real, Par) = getChi_y(State.x[2], State.x[3], T, Par)

function getChi_z(iSigma::AbstractArray, Gamma::AbstractArray, T::Real, Par)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iG_ab(x, w) = iG_(iSigma, x, w, T)
    V_n(n, Rij, s, t, u) = V_(Gamma, n, s, t, u, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iG_ab(xi, nK)[1, 1] * iG_ab(xi, nK)[2, 2]
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                for g1 = 0:2, g2 = 0:2, g3 = 0:2, g4 = 0:2
                    GGGG = (
                        iG_ab(xi, nK)[1, g1+1] *
                        iG_ab(xi, nK)[2, g2+1] *
                        iG_ab(xj, nK2)[1, g3+1] *
                        iG_ab(xj, nK2)[2, g4+1]
                    )
                    Chi[Rij] +=
                        GGGG * V_n(IndexFromFlavors(g1, g2, g3, g4), Rij, 0, npwpw2, -w2mw)
                end
            end
        end
    end
    return (Chi)
end

function getChi_x(iSigma::AbstractArray, Gamma::AbstractArray, T::Real, Par)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iG_ab(x, w) = iG_(iSigma, x, w, T)
    V_n(n, Rij, s, t, u) = V_(Gamma, n, s, t, u, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)
    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iG_ab(xi, nK)[2, 2] * iG_ab(xi, nK)[3, 3]
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                for g1 = 0:2, g2 = 0:2, g3 = 0:2, g4 = 0:2
                    GGGG = (
                        iG_ab(xi, nK)[2, g1+1] *
                        iG_ab(xi, nK)[3, g2+1] *
                        iG_ab(xj, nK2)[2, g3+1] *
                        iG_ab(xj, nK2)[3, g4+1]
                    )
                    Chi[Rij] +=
                        GGGG * V_n(IndexFromFlavors(g1, g2, g3, g4), Rij, 0, npwpw2, -w2mw)
                end
            end
        end
    end
    return (Chi)
end

function getChi_y(iSigma::AbstractArray, Gamma::AbstractArray, T::Real, Par)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iG_ab(x, w) = iG_(iSigma, x, w, T)
    V_n(n, Rij, s, t, u) = V_(Gamma, n, s, t, u, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iG_ab(xi, nK)[3, 3] * iG_ab(xi, nK)[1, 1]
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                for g1 = 0:2, g2 = 0:2, g3 = 0:2, g4 = 0:2
                    GGGG = (
                        iG_ab(xi, nK)[3, g1+1] *
                        iG_ab(xi, nK)[1, g2+1] *
                        iG_ab(xj, nK2)[3, g3+1] *
                        iG_ab(xj, nK2)[1, g4+1]
                    )
                    Chi[Rij] +=
                        GGGG * V_n(IndexFromFlavors(g1, g2, g3, g4), Rij, 0, npwpw2, -w2mw)
                end
            end
        end
    end
    return (Chi)
end

export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z, TestFRGAt

end
