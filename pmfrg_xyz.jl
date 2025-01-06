### pinned SpinFRGLattices at v0.5.2
import Pkg;
Pkg.add("SpinFRGLattices")

Pkg.add("OrdinaryDiffEq")
Pkg.add("RecursiveArrayTools")
Pkg.add("CairoMakie")
Pkg.add("DiffEqCallbacks")

using RecursiveArrayTools

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

setZero!(a::AbstractArray{T,N}) where {T,N} = fill!(a,zero(T))

function setZero!(PartArr::ArrayPartition)
    for arr in PartArr.x
        fill!(arr,0.)
    end
end
# for arrays like [[1,2,3], [4,5,6]]

"""Recursively sets structure to zero"""
function setZero!(a::T) where T
    for f in fieldnames(T)
        setZero!( getfield(a,f))
    end
    return a
end

### A vertex will now be an Array{T, 5} with the first index
### labelling the flavor

# struct VertexType{T}
#     a::Array{T, 4}
#     b::Array{T, 4}
#     c::Array{T, 4}
# end

### There will now be 24 Bubbles (yey!)
### Also stored as an Array{T, 5}

# struct BubbleType{T}
#     a::Array{T, 4}
#     b::Array{T, 4}
#     c::Array{T, 4}

#     Ta::Array{T, 4}
#     Tb::Array{T, 4}
#     Tc::Array{T, 4}
#     Td::Array{T, 4}
# end

struct SigmaType{T}
    x::Array{T, 2}
    y::Array{T, 2}
    z::Array{T, 2}
end

struct StateType{T}
    f_int::Array{T}         ### additional index in f, Sigma and Gamma for inequivalent sites (x in geometry package)
    iSigma::SigmaType{T}
    Gamma::Array{T, 5}
end

struct NumericalParams
    T::Real             ### temperature
    N::Integer          ### number of matsubara freqs

    accuracy::Real
    lambda_min::Real
    lambda_max::Real

    lenIntw::Int
    lenIntw_acc::Int
end

struct OptionParams
    use_symmetry::Bool
    minimal_output::Bool
end

struct OneLoopParams
    System
    NumericalParams::NumericalParams
    Options::OptionParams
end

struct OneLoopWorkspace{T}
    State::StateType{T} ### Stores the current state
    Deriv::StateType{T} ### Stores the derivative
    X::Array{T, 5}      ### Stores the bubble function X and XTilde
    Par                 ### Params
end

getVDims(Par) = (21, Par.System.Npairs, Par.NumericalParams.N, Par.NumericalParams.N, Par.NumericalParams.N)
getBubbleVDims(Par) = (42, Par.System.Npairs, Par.NumericalParams.N, Par.NumericalParams.N, Par.NumericalParams.N)
_getFloatType(Par) = typeof(Par.NumericalParams.T)

# function VertexType(VDims::Tuple)
#     return VertexType(
#         zeros(VDims),
#         zeros(VDims),
#         zeros(VDims)
#     )
# end
# VertexType(Par) = VertexType(getVDims(Par))

# function BubbleType(VDims::Tuple, type=Float64)
#     return BubbleType(
#         zeros(type, VDims),
#         zeros(type, VDims),
#         zeros(type, VDims),

#         zeros(type, VDims),
#         zeros(type, VDims),
#         zeros(type, VDims),
#         zeros(type, VDims)
#     )
# end
# BubbleType(Par) = BubbleType(getVDims(Par))

function SigmaType(NUnique::Int, N::Int, type=Float64)
    return SigmaType(
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
        zeros(type, NUnique, N)
    )
end
SigmaType(Par) = SigmaType(Par.System.NPairs, Par.NumericalParams.N)

function StateType(NUnique::Int, N::Int, VDims::Tuple, type=Float64)
    return StateType(
        zeros(type, NUnique),
        SigmaType(tpye, NUnique, N),
        zeros(type, VDims)
    )
end
StateType(Par) = StateType(Par.System.NUnique, Par.NumericalParams.N, getVDims(Par), _getFloatType(Par))
StateType(f_int, iSigma_x, iSigma_y, iSigma_z, Gamma) = StateType(f_int, SigmaType(iSigma_x, iSigma_y, iSigma_z), Gamma)
RecursiveArrayTools.ArrayPartition(x) = ArrayPartition(x.f_int, x.iSigma.x, x.iSigma.y, x.iSigma.z, x.Gamma)
StateType(Arr::ArrayPartition) = StateType(Arr.x...)

function NumericalParams(;
    T::Real = 0.5, # Temperature
    N::Integer = 24,

    accuracy = 1e-6,
    lambda_min = exp(-10.),
    lambda_max = exp(10.),

    lenIntw::Int = N,
    lenIntw_acc::Int = 2*maximum((N, lenIntw))
    )

    return NumericalParams(
        T,
        N,

        accuracy,
        lambda_min,
        lambda_max,

        lenIntw,
        lenIntw_acc
    )
end

OptionParams(;use_symmetry::Bool = true,MinimalOutput::Bool = false,kwargs...) = OptionParams(use_symmetry,MinimalOutput)
Params(System;kwargs...) = OneLoopParams(System,NumericalParams(;kwargs...),OptionParams(;kwargs...))

function OneLoopWorkspace(Deriv,State,X,Par)
    setZero!(Deriv)
    setZero!(X)
    return OneLoopWorkspace(
        StateType(State.x...),
        StateType(Deriv.x...),
        X,
        Par
    )
end

#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

### Propagators will depend on an additional flavor
### Instead of modifying the propagators, I will simply use them
### as V_, by doing iG_(iSigma.x, ...)

function get_w(nw, T)
    return pi * T * (2 * nw + 1)
end

function get_sign_iw(nw::Integer,N::Integer)
# s = sign(nw)
nw_bounds = min(nw, N - 1)  ### used to be min(abs(nw),...), but nw is set positive in gamma
return nw_bounds + 1        ### used to be s * ...
end

### Sigma inputted as State.iSigma, which is Array{T, 2}
function iSigma_(iSigma::AbstractArray, x::Integer, nw::Integer)
    N = size(iSigma, 2)
    s = 1
    if nw < 0
        nw = -nw - 1
        s = -1
    end
    iw = get_sign_iw(nw, N)
    return s * iSigma[x, iw]
end

function iG_(iSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw,T)
    return w / (w^2 + w * iSigma_(iSigma, x, nw) + Lam^2)
end

function iS_(iSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw, T)
    return -iG_(iSigma, x, Lam, nw, T)^2 * 2 * Lam / w
end

function iSKat_(iSigma::AbstractArray, DSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw, T)
    return -iG_(iSigma, x, Lam, nw, T)^2 * (2 * Lam / w + iSigma_(DSigma, x, nw))
end

####################################################
######### VERTICES ## VERTICES ## VERTICES #########
####################################################

### Symmetries: 
###     s <--> -s
###     t <--> -t, i <--> j
###     u <--> -u, i <--> j
function ConvertFreqArgs(ns, nt, nu, Nw)
    swapsites = nt * nu < 0
    ns, nt, nu = abs.((ns, nt, nu))

    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2) ### weird cutoff, idk why
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)

    return ns, nt, nu, swapsites
end

using LinearAlgebra

function V_(Vertex::AbstractArray, ns::Int, nt::Int, nu::Int, Rij::Integer, Rji::Integer, N::Integer)

    ### Compute the matrices, lets say M is found
    ### Then return M * Vertex[:, Rij, ns + 1, nt + 1, nu + 1]

    MS_ = [
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### xy1
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### xz1
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### yx1
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### yz1
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; ### zx1
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; ### zy1
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; ### xy2
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; ### xz2
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; ### yx2
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; ### yz2
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; ### zx2
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; ### zy2
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; ### xy3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ### xz3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; ### yx3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; ### yz3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; ### zx3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0  ### zy3
    ]

    MT_ = [
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### xy1
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; ### xz1
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### yx1
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; ### yz1
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### zx1
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### zy1
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; ### xy2
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; ### xz2
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; ### yx2
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; ### yz2
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; ### zx2
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; ### zy2
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; ### xy3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ### xz3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; ### yx3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; ### yz3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; ### zx3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0  ### zy3
    ]

    MU_ = [
        1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### xy1
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; ### xz1
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### yx1
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; ### yz1
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### zx1
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ### zy1
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; ### xy2
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; ### xz2
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; ### yx2
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; ### yz2
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; ### zx2
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; ### zy2
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; ### xy3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; ### xz3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; ### yx3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; ### yz3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; ### zx3
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1; ### zy3
    ]

    FlavorTransform = Matrix{Float64}(I, 21, 21)

    if ns < 0
        FlavorTransform = MS_ * FlavorTransform
    end
    if nt < 0
        FlavorTransform = MT_ * FlavorTransform
    end
    if nu < 0
        FlavorTransform = MU_ * FlavorTransform
    end

    ns, nt, nu, swapsites = ConvertFreqArgs(ns, nt, nu, N)
    Rij = ifelse(swapsites, Rji, Rij)
    return FlavorTransform * Vertex[:, Rij, ns+1, nt+1, nu+1]
end

function mixedFrequencies(ns,nt,nu,nwpr)
	nw1=Int((ns + nt + nu - 1) / 2)
    nw2=Int((ns - nt - nu - 1) / 2)
    nw3=Int((-ns + nt - nu - 1) / 2)
    nw4=Int((-ns - nt + nu - 1) / 2)

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

fd_ = Dict(
    "xx" => 1,
    "yy" => 2,
    "zz" => 3,
    "xy1" => 4,
    "xz1" => 5,
    "yx1" => 6,
    "yz1" => 7,
    "zx1" => 8,
    "zy1" => 9,
    "xy2" => 10,
    "xz2" => 11,
    "yx2" => 12,
    "yz2" => 13,
    "zx2" => 14,
    "zy2" => 15,
    "xy3" => 16,
    "xz3" => 17,
    "yx3" => 18,
    "yz3" => 19,
    "zx3" => 20,
    "zy3" => 21
)

### The X bubble containing the k sum
function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
	(; State, X, Par) = Workspace
	N = Par.NumericalParams.N
	(; Npairs, Nsum, siteSum, invpairs) = Par.System

    Vert(Rij, s, t, u) = V_(State.Gamma, s, t, u, Rij, invpairs[Rij], N)
	ns = is - 1 ### because julia indexes from 1
	nt = it - 1
	nu = iu - 1
	wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

	# get fields of siteSum struct as Matrices for better use of LoopVectorization
	S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m

	for Rij in 1:Npairs
		#loop over all left hand side inequivalent pairs Rij
        X_sum = zeros(42)
		for k_spl in 1:Nsum[Rij]
			#loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
			ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]
			Ptm = Props[xk,xk,:,:]*m ### Props now contains two flavor indices

            V12 = Vert(ki, ns, wpw1, -wpw2)
            V34 = Vert(kj, ns, -wmw3, -wmw4)

            X_sum[fd_["xx"]] += -V12[fd_["xx"]] * V34[fd_["xx"]] * Ptm[1, 1] - V12[fd_["xy1"]] * V34[fd_["yx1"]] * Ptm[2, 2] - V12[fd_["xz1"]] * V34[fd_["zx1"]] * Ptm[3, 3]
            X_sum[fd_["yy"]] += -V12[fd_["yy"]] * V34[fd_["yy"]] * Ptm[2, 2] - V12[fd_["yz1"]] * V34[fd_["zy1"]] * Ptm[3, 3] - V12[fd_["yx1"]] * V34[fd_["xy1"]] * Ptm[1, 1]
            X_sum[fd_["zz"]] += -V12[fd_["zz"]] * V34[fd_["zz"]] * Ptm[3, 3] - V12[fd_["zx1"]] * V34[fd_["xz1"]] * Ptm[1, 1] - V12[fd_["zy1"]] * V34[fd_["yz1"]] * Ptm[2, 2]

            ### Xab1 = -Vaa Vab1 - Vab1 Vbb - Vac1 Vcb1
            X_sum[fd_["xy1"]] += -V12[fd_["xx"]] * V34[fd_["xy1"]] * Ptm[1, 1] - V12[fd_["xy1"]] * V34[fd_["yy"]] * Ptm[2, 2] - V12[fd_["xz1"]] * V34[fd_["zy1"]] * Ptm[3, 3]
            X_sum[fd_["xz1"]] += -V12[fd_["xx"]] * V34[fd_["xz1"]] * Ptm[1, 1] - V12[fd_["xz1"]] * V34[fd_["zz"]] * Ptm[3, 3] - V12[fd_["xy1"]] * V34[fd_["yz1"]] * Ptm[2, 2]
            X_sum[fd_["yx1"]] += -V12[fd_["yy"]] * V34[fd_["yx1"]] * Ptm[2, 2] - V12[fd_["yx1"]] * V34[fd_["xx"]] * Ptm[1, 1] - V12[fd_["yz1"]] * V34[fd_["zx1"]] * Ptm[3, 3]
            X_sum[fd_["yz1"]] += -V12[fd_["yy"]] * V34[fd_["yz1"]] * Ptm[2, 2] - V12[fd_["yz1"]] * V34[fd_["zz"]] * Ptm[3, 3] - V12[fd_["yx1"]] * V34[fd_["xz1"]] * Ptm[1, 1]
            X_sum[fd_["zx1"]] += -V12[fd_["zz"]] * V34[fd_["zx1"]] * Ptm[3, 3] - V12[fd_["zx1"]] * V34[fd_["xx"]] * Ptm[1, 1] - V12[fd_["zy1"]] * V34[fd_["yx1"]] * Ptm[2, 2]
            X_sum[fd_["zy1"]] += -V12[fd_["zz"]] * V34[fd_["zy1"]] * Ptm[3, 3] - V12[fd_["zy1"]] * V34[fd_["yy"]] * Ptm[2, 2] - V12[fd_["zx1"]] * V34[fd_["xy1"]] * Ptm[1, 1]
            
            ### Xab2 = -Vab3 Vab2 - Vab2 Vba3
            X_sum[fd_["xy2"]] += -V12[fd_["xy3"]] * V34[fd_["xy2"]] * Ptm[2, 1] - V12[fd_["xy2"]] * V34[fd_["yx3"]] * Ptm[1, 2]
            X_sum[fd_["xz2"]] += -V12[fd_["xz3"]] * V34[fd_["xz2"]] * Ptm[3, 1] - V12[fd_["xz2"]] * V34[fd_["zx3"]] * Ptm[1, 3]
            X_sum[fd_["yx2"]] += -V12[fd_["yx3"]] * V34[fd_["yx2"]] * Ptm[1, 2] - V12[fd_["yx2"]] * V34[fd_["xy3"]] * Ptm[2, 1]
            X_sum[fd_["yz2"]] += -V12[fd_["yz3"]] * V34[fd_["yz2"]] * Ptm[3, 2] - V12[fd_["yz2"]] * V34[fd_["zy3"]] * Ptm[2, 3]
            X_sum[fd_["zx2"]] += -V12[fd_["zx3"]] * V34[fd_["zx2"]] * Ptm[1, 3] - V12[fd_["zx2"]] * V34[fd_["xz3"]] * Ptm[3, 1]
            X_sum[fd_["zy2"]] += -V12[fd_["zy3"]] * V34[fd_["zy2"]] * Ptm[2, 3] - V12[fd_["zy2"]] * V34[fd_["yz3"]] * Ptm[3, 2]

            ### Xab3 = -Vab3 Vab3 - Vab2 Vba2
            X_sum[fd_["xy3"]] += -V12[fd_["xy3"]] * V34[fd_["xy3"]] * Ptm[2, 1] - V12[fd_["xy2"]] * V34[fd_["yx2"]] * Ptm[1, 2]
            X_sum[fd_["xz3"]] += -V12[fd_["xz3"]] * V34[fd_["xz3"]] * Ptm[3, 1] - V12[fd_["xz2"]] * V34[fd_["zx2"]] * Ptm[1, 3]
            X_sum[fd_["yx3"]] += -V12[fd_["yx3"]] * V34[fd_["yx3"]] * Ptm[1, 2] - V12[fd_["yx2"]] * V34[fd_["xy2"]] * Ptm[2, 1]
            X_sum[fd_["yz3"]] += -V12[fd_["yz3"]] * V34[fd_["yz3"]] * Ptm[3, 2] - V12[fd_["yz2"]] * V34[fd_["zy2"]] * Ptm[2, 3]
            X_sum[fd_["zx3"]] += -V12[fd_["zx3"]] * V34[fd_["zx3"]] * Ptm[1, 3] - V12[fd_["zx2"]] * V34[fd_["xz2"]] * Ptm[3, 1]
            X_sum[fd_["zy3"]] += -V12[fd_["zy3"]] * V34[fd_["zy3"]] * Ptm[2, 3] - V12[fd_["zy2"]] * V34[fd_["yz2"]] * Ptm[3, 2]

			# X_sum .*= Ptm ### attention. Something could go wrong here
            # I KNEW IT!!!
		end

		X[:, Rij, is, it, iu] .+= X_sum ### same here
    end
    return
end

function addY!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
	(; State, X, Par) = Workspace
	N = Par.NumericalParams.N
	(; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    Vert(Rij, s, t, u) = V_(State.Gamma, s, t, u, Rij, invpairs[Rij], N)
	ns = is - 1 ### because julia indexes from 1
	nt = it - 1
	nu = iu - 1
	wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

	# Xtilde only defined for nonlocal pairs Rij != Rii
	for Rij in 1:Npairs
		Rij in OnsitePairs && continue
		# loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij) 
		(; xi, xj) = PairTypes[Rij]

        function P_(n::Int, m::Int)
            return Props[xi, xj, n, m]
        end

        function PT_(n::Int, m::Int)
            return Props[xj, xi, m, n]
        end

        V13 = Vert(Rij, -wmw1, nt, wmw3)
        V24 = Vert(Rij, wpw2, -nt, -wpw4) ### potential problems due to -nt?

        V31 = Vert(Rij, wmw3, nt, -wmw1)
        V42 = Vert(Rij, -wpw4, -nt, wpw2)

        X_sum = zeros(42, Npairs, N, N, N)

        ### Yaa = Vaa Vaa + Vab2 Vab2 + Vac2 Vac2 + (w -- -w + t) 

        X_sum[21 + fd_["xx"], Rij, is, it, iu] += ( ### add 21 because this is X_sum[21 + n] = Y[n]
            (V13[fd_["xx"]] * V24[fd_["xx"]] * P_(1, 1)
            + V13[fd_["xy2"]] * V24[fd_["xy2"]] * P_(2, 2) 
            + V13[fd_["xz2"]] * V24[fd_["xz2"]] * P_(3, 3)) ### Blind-copied props, potentially wrong
        
            + (V31[fd_["xx"]] * V42[fd_["xx"]] * PT_(1, 1)
            + V31[fd_["xy2"]] * V42[fd_["xy2"]] * PT_(2, 2) 
            + V31[fd_["xz2"]] * V42[fd_["xz2"]] * PT_(2, 2))
        )

        X_sum[21 + fd_["yy"], Rij, is, it, iu] += (
            (V13[fd_["yy"]] * V24[fd_["yy"]] * P_(2, 2)
            + V13[fd_["yx2"]] * V24[fd_["yx2"]] * P_(1, 1)  
            + V13[fd_["yz2"]] * V24[fd_["yz2"]] * P_(3, 3)) ### Blind-copied props, potentially wrong
        
            + (V31[fd_["yy"]] * V42[fd_["yy"]] * PT_(2, 2)
            + V31[fd_["yx2"]] * V42[fd_["yx2"]] * PT_(1, 1) 
            + V31[fd_["yz2"]] * V42[fd_["yz2"]] * PT_(3, 3))
        )

        X_sum[21 + fd_["zz"], Rij, is, it, iu] += (
            (V13[fd_["zz"]] * V24[fd_["zz"]] * P_(3, 3)
            + V13[fd_["zx2"]] * V24[fd_["zx2"]] * P_(1, 1) 
            + V13[fd_["zy2"]] * V24[fd_["zy2"]] * P_(2, 2)) ### Blind-copied props, potentially wrong
        
            + (V31[fd_["zz"]] * V42[fd_["zz"]] * PT_(3, 3)
            + V31[fd_["zx2"]] * V42[fd_["zx2"]] * PT_(1, 1)
            + V31[fd_["zy2"]] * V42[fd_["zy2"]] * PT_(2, 2))
        )

        ### Yab1 = Vab3 Vab3 + Vab1 Vab1 + (w -- -w + t)

        X_sum[21 + fd_["xy1"], Rij, is, it, iu] += (
            (V13[fd_["xy3"]] * V24[fd_["xy3"]] * P_(2, 1)
            + V13[fd_["xy1"]] * V24[fd_["xy1"]] * P_(1, 2))

            + (V31[fd_["xy3"]] * V42[fd_["xy3"]] * PT_(2, 1)
            + V31[fd_["xy1"]] * V42[fd_["xy1"]] * PT_(1, 2))
        )

        X_sum[21 + fd_["xz1"], Rij, is, it, iu] += (
            (V13[fd_["xz3"]] * V24[fd_["xz3"]] * P_(3, 1)
            + V13[fd_["xz1"]] * V24[fd_["xz1"]] * P_(1, 3))

            + (V31[fd_["xz3"]] * V42[fd_["xz3"]] * PT_(3, 1)
            + V31[fd_["xz1"]] * V42[fd_["xz1"]] * PT_(1, 3))
        )

        X_sum[21 + fd_["yx1"], Rij, is, it, iu] += (
            (V13[fd_["yx3"]] * V24[fd_["yx3"]] * P_(1, 2)
            + V13[fd_["yx1"]] * V24[fd_["yx1"]] * P_(2, 1))

            + (V31[fd_["yx3"]] * V42[fd_["yx3"]] * PT_(1, 2)
            + V31[fd_["yx1"]] * V42[fd_["yx1"]] * PT_(2, 1))
        )

        X_sum[21 + fd_["yz1"], Rij, is, it, iu] += (
            (V13[fd_["yz3"]] * V24[fd_["yz3"]] * P_(3, 2)
            + V13[fd_["yz1"]] * V24[fd_["yz1"]] * P_(2, 3))

            + (V31[fd_["yz3"]] * V42[fd_["yz3"]] * PT_(3, 2)
            + V31[fd_["yz1"]] * V42[fd_["yz1"]] * PT_(2, 3))
        )

        X_sum[21 + fd_["zx1"], Rij, is, it, iu] += (
            (V13[fd_["zx3"]] * V24[fd_["zx3"]] * P_(1, 3)
            + V13[fd_["zx1"]] * V24[fd_["zx1"]] * P_(3, 1))

            + (V31[fd_["zx3"]] * V42[fd_["zx3"]] * PT_(1, 3)
            + V31[fd_["zx1"]] * V42[fd_["zx1"]] * PT_(3, 1))
        )

        X_sum[21 + fd_["zy1"], Rij, is, it, iu] += (
            (V13[fd_["zy3"]] * V24[fd_["zy3"]] * P_(2, 3)
            + V13[fd_["zy1"]] * V24[fd_["zy1"]] * P_(3, 2))

            + (V31[fd_["zy3"]] * V42[fd_["zy3"]] * PT_(2, 3)
            + V31[fd_["zy1"]] * V42[fd_["zy1"]] * PT_(3, 2))
        )

        ### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

        X_sum[21 + fd_["xy2"], Rij, is, it, iu] += (
            (V13[fd_["xx"]] * V24[fd_["yx2"]] * P_(1, 1)
            + V13[fd_["xy2"]] * V24[fd_["yy"]] * P_(2, 2)
            + V13[fd_["xz2"]] * V24[fd_["yz2"]] + P_(3, 3))

            + (V31[fd_["xx"]] * V42[fd_["yx2"]] * PT_(1, 1)
            + V31[fd_["xy2"]] * V42[fd_["yy"]] * PT_(2, 2)
            + V31[fd_["xz2"]] * V42[fd_["yz2"]] * PT_(3, 3))
        )

        X_sum[21 + fd_["xz2"], Rij, is, it, iu] += (
            (V13[fd_["xx"]] * V24[fd_["zx2"]] * P_(1, 1)
            + V13[fd_["xz2"]] * V24[fd_["zz"]] * P_(3, 3)
            + V13[fd_["xy2"]] * V24[fd_["zy2"]] * P_(2, 2))

            + (V31[fd_["xx"]] * V42[fd_["zx2"]] * PT_(1, 1)
            + V31[fd_["xz2"]] * V42[fd_["zz"]] * PT_(3, 3)
            + V31[fd_["xy2"]] * V42[fd_["zy2"]] * PT_(2, 2))
        )

        X_sum[21 + fd_["yx2"], Rij, is, it, iu] += (
            (V13[fd_["yy"]] * V24[fd_["xy2"]] * P_(2, 2)
            + V13[fd_["yx2"]] * V24[fd_["xx"]] * P_(1, 1)
            + V13[fd_["yz2"]] * V24[fd_["xz2"]] * P_(3, 3))

            + (V31[fd_["yy"]] * V42[fd_["xy2"]] * PT_(2, 2)
            + V31[fd_["yx2"]] * V42[fd_["xx"]] * PT_(1, 1)
            + V31[fd_["yz2"]] * V42[fd_["xz2"]] * PT_(3, 3))
        )

        X_sum[21 + fd_["yz2"], Rij, is, it, iu] += (
            (V13[fd_["yy"]] * V24[fd_["zy2"]] * P_(2, 2)
            + V13[fd_["yz2"]] * V24[fd_["zz"]] * P_(3, 3)
            + V13[fd_["yx2"]] * V24[fd_["zx2"]] * P_(1, 1))

            + (V31[fd_["yy"]] * V42[fd_["zy2"]] * PT_(2, 2)
            + V31[fd_["yz2"]] * V42[fd_["zz"]] * PT_(3, 3)
            + V31[fd_["yx2"]] * V42[fd_["zx2"]] * PT_(1, 1))
        )

        X_sum[21 + fd_["zx2"], Rij, is, it, iu] += (
            (V13[fd_["zz"]] * V24[fd_["xz2"]] * P_(3, 3)
            + V13[fd_["zx2"]] * V24[fd_["xx"]] * P_(1, 1)
            + V13[fd_["zy2"]] * V24[fd_["xy2"]] * P_(2, 2))

            + (V31[fd_["zz"]] * V42[fd_["xz2"]] * PT_(3, 3)
            + V31[fd_["zx2"]] * V42[fd_["xx"]] * PT_(1, 1)
            + V31[fd_["zy2"]] * V42[fd_["xy2"]] * PT_(2, 2))
        )

        X_sum[21 + fd_["zy2"], Rij, is, it, iu] += (
            (V13[fd_["zz"]] * V24[fd_["yz2"]] * P_(3, 3)
            + V13[fd_["zy2"]] * V24[fd_["yy"]] * P_(2, 2)
            + V13[fd_["zx2"]] * V24[fd_["yx2"]] * P_(1, 1))

            + (V13[fd_["zz"]] * V24[fd_["yz2"]] * PT_(3, 3)
            + V13[fd_["zy2"]] * V24[fd_["yy"]] * PT_(2, 2)
            + V13[fd_["zx2"]] * V24[fd_["yx2"]] * PT_(1, 1))
        )

        ### Yab3 = Vab3 Vba1 + Vab1 Vba3 + (w -- -w + t)

        X_sum[21 + fd_["xy3"], Rij, is, it, iu] += (
            (V13[fd_["xy3"]] * V24[fd_["yx1"]] * P_(2, 1)
            + V13[fd_["xy1"]] * V24[fd_["yx3"]] * P_(1, 2))

            + (V31[fd_["xy3"]] * V42[fd_["yx1"]] * PT_(2, 1)
            + V31[fd_["xy1"]] * V42[fd_["yx3"]] * PT_(1, 2)) 
        )

        X_sum[21 + fd_["xz3"], Rij, is, it, iu] += (
            (V13[fd_["xz3"]] * V24[fd_["zx1"]] * P_(3, 1)
            + V13[fd_["xz1"]] * V24[fd_["zx3"]] * P_(1, 3)) 

            + (V31[fd_["xz3"]] * V42[fd_["zx1"]] * PT_(3, 1)
            + V31[fd_["xz1"]] * V42[fd_["zx3"]] * PT_(1, 3)) 
        )

        X_sum[21 + fd_["yx3"], Rij, is, it, iu] += (
            (V13[fd_["yx3"]] * V24[fd_["xy1"]] * P_(1, 2)
            + V13[fd_["yx1"]] * V24[fd_["xy3"]] * P_(2, 1)) 

            + (V31[fd_["yx3"]] * V42[fd_["xy1"]] * PT_(1, 2)
            + V31[fd_["yx1"]] * V42[fd_["xy3"]] * PT_(2, 1)) 
        )

        X_sum[21 + fd_["yz3"], Rij, is, it, iu] += (
            (V13[fd_["yz3"]] * V24[fd_["zy1"]] * P_(3, 2)
            + V13[fd_["yz1"]] * V24[fd_["zy3"]] * P_(2, 3)) 

            + (V31[fd_["yz3"]] * V42[fd_["zy1"]] * PT_(3, 2)
            + V31[fd_["yz1"]] * V42[fd_["zy3"]] * PT_(2, 3)) 
        )

        X_sum[21 + fd_["zx3"], Rij, is, it, iu] += (
            (V13[fd_["zx3"]] * V24[fd_["xz1"]] * P_(1, 3)
            + V13[fd_["zx1"]] * V24[fd_["xz3"]] * P_(3, 1)) 

            + (V31[fd_["zx3"]] * V42[fd_["xz1"]] * PT_(1, 3)
            + V31[fd_["zx1"]] * V42[fd_["xz3"]] * PT_(3, 1)) 
        )

        X_sum[21 + fd_["zy3"], Rij, is, it, iu] += (
            (V13[fd_["zy3"]] * V24[fd_["yz1"]] * P_(2, 3)
            + V13[fd_["zy1"]] * V24[fd_["yz3"]] * P_(3, 2)) 

            + (V31[fd_["zy3"]] * V42[fd_["yz1"]] * PT_(2, 3)
            + V31[fd_["zy1"]] * V42[fd_["yz3"]] * PT_(3, 2)) 
        )

        X .+= X_sum
    end
end

function getXBubble!(Workspace, Lam)
	Par = Workspace.Par
    (; T, N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System
	 
	iGx(x, nw) = iG_(Workspace.State.iSigma.x, x, Lam, nw, T)
    iGy(x, nw) = iG_(Workspace.State.iSigma.y, x, Lam, nw, T)
    iGz(x, nw) = iG_(Workspace.State.iSigma.z, x, Lam, nw, T)

    ### Does flavor type commute with differentiation in 1 / (G^-1 + iSigma) ? Probably
	iSKatx(x,nw) = iSKat_(Workspace.State.iSigma.x, Workspace.Deriv.iSigma.x, x, Lam, nw, T)
	iSKaty(x,nw) = iSKat_(Workspace.State.iSigma.y, Workspace.Deriv.iSigma.y, x, Lam, nw, T)
	iSKatz(x,nw) = iSKat_(Workspace.State.iSigma.z, Workspace.Deriv.iSigma.z, x, Lam, nw, T)

	function getKataninProp!(BubbleProp,nw1,nw2)
		for i in 1:Par.System.NUnique, j in 1:Par.System.NUnique
			BubbleProp[i, j, 1, 1] = iSKatx(i, nw1) * iGx(j, nw2) * T
			BubbleProp[i, j, 1, 2] = iSKatx(i, nw1) * iGy(j, nw2) * T
			BubbleProp[i, j, 1, 3] = iSKatx(i, nw1) * iGz(j, nw2) * T
			BubbleProp[i, j, 2, 1] = iSKaty(i, nw1) * iGx(j, nw2) * T
			BubbleProp[i, j, 2, 2] = iSKaty(i, nw1) * iGy(j, nw2) * T
			BubbleProp[i, j, 2, 3] = iSKaty(i, nw1) * iGz(j, nw2) * T
			BubbleProp[i, j, 3, 1] = iSKatz(i, nw1) * iGx(j, nw2) * T
			BubbleProp[i, j, 3, 2] = iSKatz(i, nw1) * iGy(j, nw2) * T
			BubbleProp[i, j, 3, 3] = iSKatz(i, nw1) * iGz(j, nw2) * T
		end

        return BubbleProp
		# return SMatrix{NUnique, NUnique, 3, 3}(BubbleProp)
        ### SMatrix can only create 2d array (according to ChatGPT). Use SArray instead
	end

	for is in 1:N, it in 1:N
        BubbleProp = zeros(NUnique, NUnique, 3, 3)
        ns = is - 1
        nt = it - 1
        for nw in -lenIntw:lenIntw-1 # Matsubara sum
            sprop = getKataninProp!(BubbleProp,nw,nw+ns)
            for iu in 1:N
                nu = iu - 1
                if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
                    continue
                end
                addY!(Workspace, is, it, iu, nw, sprop) # add to XTilde-type bubble functions

                ### If no u--t symmetry, then add all the bubbles
                ### If use u--t symmetry, then only add for nu smaller then nt (all other obtained by symmetry)
                # if(!Par.Options.use_symmetry || nu<=nt)
                addX!(Workspace, is, it, iu, nw, sprop)
                # end
            end
        end
	end
end

function symmetrizeBubble!(X::Array{T, 5}, Par) where {T}
    N = Par.NumericalParams.N
    (;Npairs,OnsitePairs) = Par.System
    use_symmetry = Par.Options.use_symmetry
    # use the u <--> t symmetry
    if(use_symmetry)
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
    for iu in 1:N
		for it in 1:N, is in 1:N, R in OnsitePairs
            X[21 + 1, R, is, it, iu] = X[1, R, it, is, iu]  ###
            X[21 + 2, R, is, it, iu] = X[2, R, it, is, iu]  ### Yaa = Xaa
            X[21 + 3, R, is, it, iu] = X[3, R, it, is, iu]  ###
            for n in 1:6
                X[21 + 3 + n, R, is, it, iu] = X[9 + n, R, it, is, iu]      ### Yab1 = Xab2
                X[21 + 9 + n, R, is, it, iu] = X[3 + n, R, it, is, iu]      ### Yab2 = Xab1
                X[21 + 15 + n, R, is, it, iu] = X[15 + n, R, it, is, iu]    ### Yab3 = Xab3
            end
			# X.Ta[R,is,it,iu] = X.a[R,is,it,iu]
			# X.Tb[R,is,it,iu] = X.b[R,is,it,iu]
			# X.Tc[R,is,it,iu] = X.c[R,is,it,iu]
			# X.Td[R,is,it,iu] = -X.c[R,is,iu,it]
		end
    end
    # X.Td .= X.Ta .- X.Tb .- X.Tc
end

### Gamma = Gamma[flavortype, Rij, is, it, iu]
function addToVertexFromBubble!(Gamma::Array{T, 5}, X::Array{T, 5}) where {T}
    for iu in axes(Gamma,5)
        for it in axes(Gamma,4), is in axes(Gamma,3), Rij in axes(Gamma,2)
            for n in 1:9 ### Zaa(s,t,u) = -Yaa(s,u,t) ; Zab1(s,t,u) = -Yab1(s,u,t)
                Gamma[n, Rij, is, it, iu] += (
                    X[n, Rij, is, it, iu] 
                    + X[21 + n, Rij, is, it, iu] 
                    - X[21 + n, Rij, is, iu, it]
                )
            end
            for n in 1:6 ### Zab2(s,t,u) = -Yab3(s,u,t) ; Zab3(s,t,u) = -Yab2(s,u,t)
                Gamma[9 + n, Rij, is, it, iu] += (
                    X[9 + n, Rij, is, it, iu] 
                    + X[21 + 9 + n, Rij, is, it, iu] 
                    - X[21 + 15 + n, Rij, is, iu, it]
                )
                Gamma[15 + n, Rij, is, it, iu] += (
                    X[15 + n, Rij, is, it, iu] 
                    + X[21 + 15 + n, Rij, is, it, iu] 
                    - X[21 + 9 + n, Rij, is, iu, it]
                )
            end
            # Gamma.a[Rij,is,it,iu] += X.a[Rij,is,it,iu] - X.Ta[Rij,it,is,iu] + X.Ta[Rij,iu,is,it]
            # Gamma.b[Rij,is,it,iu] += X.b[Rij,is,it,iu] - X.Tc[Rij,it,is,iu] + X.Tc[Rij,iu,is,it]
            # Gamma.c[Rij,is,it,iu] += X.c[Rij,is,it,iu] - X.Tb[Rij,it,is,iu] + X.Td[Rij,iu,is,it]
        end
    end 
    return Gamma
end

function symmetrizeVertex!(Gamma::Array{T, 5},Par) where {T}
	N = Par.NumericalParams.N
	for iu in 1:N
		for it in 1:N, is in 1:N, R in Par.System.OnsitePairs
			Gamma[:, R, is, it, iu] .= -Gamma[:, R, it, is, iu]
		end
	end
end

######################################################################
######### FLOW EQUATIONS ## FLOW EQUATIONS ## FLOW EQUATIONS #########
######################################################################

function getDFint!(Workspace, Lam::Real)
    (; State, Deriv, Par) = Workspace
    (; T, lenIntw_acc) = Par.NumericalParams
    NUnique = Par.System.NUnique
	
	iSigmax(x, nw) = iSigma_(State.iSigma.x, x, nw)
	iSigmay(x, nw) = iSigma_(State.iSigma.y, x, nw)
	iSigmaz(x, nw) = iSigma_(State.iSigma.z, x, nw)

	iGx(x, nw) = iG_(State.iSigma.x, x, Lam, nw, T)
	iGy(x, nw) = iG_(State.iSigma.y, x, Lam, nw, T)
	iGz(x, nw) = iG_(State.iSigma.z, x, Lam, nw, T)

	iSx(x, nw) = iS_(State.iSigma.x, x, Lam, nw, T)
	iSy(x, nw) = iS_(State.iSigma.y, x, Lam, nw, T)
	iSz(x, nw) = iS_(State.iSigma.z, x, Lam, nw, T)

	Theta(Lam, w) = w^2 / (w^2 + Lam^2)
	
	for x in 1:NUnique
		sumres = 0.
		for nw in -lenIntw_acc:lenIntw_acc-1
			w = get_w(nw,T) ### is computed in iS, iG and iSigma too.
			sumres += iSx(x, nw) / iGy(x, nw) * Theta(Lam, w) * iSigmax(x, nw) / w # 1 / w is simply G_0 without Lambda (right?)
            sumres += iSy(x, nw) / iGy(x, nw) * Theta(Lam, w) * iSigmay(x, nw) / w
            sumres += iSz(x, nw) / iGz(x, nw) * Theta(Lam, w) * iSigmaz(x, nw) / w
        end
		Deriv.f_int[x] = -0.5 * T * sumres ### I divided the original by 3, because 3 flavors
	end
end

function get_Self_Energy!(Workspace, Lam)
	Par = Workspace.Par
	@inline iSx(x, nw) = iS_(Workspace.State.iSigma.x, x, Lam, nw, Par.NumericalParams.T) / 2
	@inline iSy(x, nw) = iS_(Workspace.State.iSigma.y, x, Lam, nw, Par.NumericalParams.T) / 2
	@inline iSz(x, nw) = iS_(Workspace.State.iSigma.z, x, Lam, nw, Par.NumericalParams.T) / 2
	compute1PartBubble!(Workspace.Deriv.iSigma, Workspace.State.Gamma, [iSx, iSy, iSz], Par)
end

function compute1PartBubble!(Dgamma::SigmaType, Gamma::Array{T, 5}, Props, Par) where {T}
    invpairs = Par.System.invpairs

	setZero!(Dgamma)

    #### WARUM HIER T, U, S ???????
	# @inline Gamma_a(Rij,s,t,u) = V_(Gamma.a, t, u, s, Rij, invpairs[Rij], Par.NumericalParams.N) # Tilde-type can be obtained by permutation of vertices
	# @inline Gamma_b(Rij,s,t,u) = V_(Gamma.b, t, u, s, Rij, invpairs[Rij], Par.NumericalParams.N) # cTilde corresponds to b type vertex!
    @inline Gamma_(Rij, s, t, u) = V_(Gamma, s, t, u, Rij, invpairs[Rij], Par.NumericalParams.N)

    addTo1PartBubble!(Dgamma, Gamma_, Props, Par)
end

function addTo1PartBubble!(Dgamma::SigmaType, Gamma_::Function, Props, Par)

    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; siteSum, Nsum, OnsitePairs) = Par.System

	Threads.@threads for iw1 in 1:N
		nw1 = iw1 - 1
    	for (x, Rx) in enumerate(OnsitePairs)
			for nw in -lenIntw_acc:lenIntw_acc-1
				jsum = zeros(3)
				wpw1 = nw1 + nw + 1 #w + w1: Adding two fermionic Matsubara frequencies gives a +1 for the bosonic index
				wmw1 = nw - nw1
				for k_spl in 1:Nsum[Rx]
					(; m, ki, xk) = siteSum[k_spl, Rx]
                    gam = Gamma_(ki, 0, -wmw1, -wpw1)
                    jsum[fd_["xx"]] += (
                        gam[fd_["xx"]] * Props[1](xk, nw)
                        + gam[fd_["yx1"]] * Props[2](xk, nw)
                        + gam[fd_["zx1"]] * Props[3](xk, nw)
                    ) * m ### different m for all flavors? Probably not
                    jsum[fd_["yy"]] += (
                        gam[fd_["xy1"]] * Props[1](xk, nw)
                        + gam[fd_["yy"]] * Props[2](xk, nw)
                        + gam[fd_["zy1"]] * Props[3](xk, nw)
                    ) * m
                    jsum[fd_["zz"]] += (
                        gam[fd_["xz1"]] * Props[1](xk, nw)
                        + gam[fd_["yz1"]] * Props[2](xk, nw)
                        + gam[fd_["zz"]] * Props[3](xk, nw)
                    ) * m
				end
				Dgamma.x[x, iw1] += -T * jsum[1] #For the self-energy derivative, the factor of 1/2 must be in the propagator
                Dgamma.y[x, iw1] += -T * jsum[2]
                Dgamma.z[x, iw1] += -T * jsum[3]
            end
		end
	end
    # return Dgamma
end

function getDeriv!(Deriv, State, setup, Lam)
    (X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(State, Deriv, X, Par)

    println("\n============ getDFint ============")
    getDFint!(Workspace, Lam)

    println("======== get_Self_Energy =========")
    get_Self_Energy!(Workspace, Lam)

    println("=========== getXBubble ===========")
    getXBubble!(Workspace, Lam)

    println("======== symmetrizeBubble ========")
    symmetrizeBubble!(Workspace.X, Par)

    println("===== addToVertexFromBubble ======")
    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.X)

    println("======== symmetrizeVertex ========")
    symmetrizeVertex!(Workspace.Deriv.Gamma, Par)
    return
end

####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################

t_to_Lam(t) = exp(t)
Lam_to_t(t) = log(t)

function AllocateSetup(Par::OneLoopParams)
    println("One Loop: T= ",Par.NumericalParams.T)
    ## Allocate Memory:
    floattype = _getFloatType(Par)
    X = zeros(floattype, getBubbleVDims(Par))
    return (X,Par)
end

function InitializeState(Par, multiCouplings)

    N = Par.NumericalParams.N;
    ( ; couplings, NUnique) = Par.System;

    VDims = getVDims(Par);
    floattype = _getFloatType(Par)

    State = ArrayPartition(
        zeros(floattype, NUnique),          ### f_int
        zeros(floattype, NUnique, N),       ### iSigma_x
        zeros(floattype, NUnique, N),       ### iSigma_y
        zeros(floattype, NUnique, N),       ### iSigma_z
        zeros(floattype, VDims),            ### Gamma
    );

    Gamma = State.x[5]
    setToBareVertex!(Gamma, multiCouplings)      ### sets initial conditions to couplings
    return State

end

function launchPMFRG!(State, setup, Deriv!::Function;
    method = DP5()
    )

    Par = setup[end]
    (; lambda_max, lambda_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(lambda_max)
    tend = get_t_min(lambda_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)
    
    println(typeof(Deriv_subst!))
    println(typeof(State))
    println(typeof(setup))

    problem = ODEProblem(Deriv_subst!, State, (t0, tend), setup) # function, initial state, timespan, ??
    sol = solve(problem, method, reltol = accuracy, abstol = accuracy, save_everystep = true, dt=Lam_to_t(0.2 * lambda_max))
    return sol
end

SolveFRG(Par, couplings; kwargs...) = launchPMFRG!(InitializeState(Par, couplings),AllocateSetup(Par),getDeriv!; kwargs...)

function get_t_min(Lam)
    Lam < exp(-30) && @warn "lambda_min too small! Set to exp(-30) instead."
    max(Lam_to_t(Lam),-30.)
end

function generateSubstituteDeriv(getDeriv!::Function)
    
    function DerivSubs!(Deriv,State,par,t)
        Lam = t_to_Lam(t)
        a = getDeriv!(Deriv,State,par,Lam)
        Deriv .*= Lam
        a
    end
    
end

function setToBareVertex!(Gamma::AbstractArray{T,5}, couplings::AbstractVector) where T
    for Rj in axes(Gamma,2)
        Gamma[fd_["yz2"], Rj, :, :, :] .= -couplings[Rj][1]
        Gamma[fd_["zy2"], Rj, :, :, :] .= -couplings[Rj][1]
        Gamma[fd_["zx2"], Rj, :, :, :] .= -couplings[Rj][2]
        Gamma[fd_["xz2"], Rj, :, :, :] .= -couplings[Rj][2]
        Gamma[fd_["xy2"], Rj, :, :, :] .= -couplings[Rj][3]
        Gamma[fd_["yx2"], Rj, :, :, :] .= -couplings[Rj][3]

        Gamma[fd_["yz3"], Rj, :, :, :] .= couplings[Rj][1]
        Gamma[fd_["zy3"], Rj, :, :, :] .= couplings[Rj][1]
        Gamma[fd_["zx3"], Rj, :, :, :] .= couplings[Rj][2]
        Gamma[fd_["xz3"], Rj, :, :, :] .= couplings[Rj][2]
        Gamma[fd_["xy3"], Rj, :, :, :] .= couplings[Rj][3]
        Gamma[fd_["yx3"], Rj, :, :, :] .= couplings[Rj][3]
    end
    return Gamma
end

#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################

struct Observables
    Chi
    gamma
end

getChi(State::ArrayPartition, Lam::Real, Par) = getChi(State.x[2], State.x[3], State.x[5], Lam, Par)

### highly unlikely to be correct. Need to check again
function getChi(iSigmaX::AbstractArray, iSigmaY::AbstractArray, Gamma::AbstractArray, Lam::Real, Par)
	(;T,N,lenIntw_acc) = Par.NumericalParams
	(;Npairs,invpairs,PairTypes,OnsitePairs) = Par.System

	iGx(x,w) = iG_(iSigmaX, x, Lam, w, T)
    iGy(x,w) = iG_(iSigmaY, x, Lam, w, T)
	Vxy2(Rij,s,t,u) = V_(Gamma,s,t,u,Rij,invpairs[Rij],N)[fd_["xy2"]]

	Chi = zeros(_getFloatType(Par),Npairs)

	for Rij in 1:Npairs
		(;xi,xj) = PairTypes[Rij]
		for nK in -lenIntw_acc:lenIntw_acc-1
			if Rij in OnsitePairs
				Chi[Rij,1] += T * iGx(xi,nK) * iGy(xi, nK)
			end
			for nK2 in -lenIntw_acc:lenIntw_acc-1
				npwpw2 = nK+nK2+1
				w2mw = nK2-nK
				#use that Vc_0 is calculated from Vb
				GGGG = iGx(xi,nK)^2 * iGy(xj,nK2)^2
				Chi[Rij] += T^2 * GGGG * Vxy2(Rij,0,npwpw2,w2mw)
			end
        end
    end
	return(Chi)
end


Pkg.add("StructArrays")
Pkg.add("JLD2")

using JLD2

##########################################################
######### DIMER SUSC ## DIMER SUSC ## DIMER SUSC #########
##########################################################

using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using CairoMakie

System = getPolymer(2) # create a structure that contains all information about the geometry of the problem.
couplings = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    T = 0.5, # Temperature for the simulation.
    N = 8, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-5,
    lambda_max = 100.,
    lambda_min = .01
)

@time sol = SolveFRG(Par, couplings, method = DP5());

## Evaluation Dimer Susc

tr = LinRange(3,-2,20)
save_object("dimer_flow_noah.jld2", [(sol(t), exp(t), Par) for t in tr])

chiR = [getChi(sol(t),exp(t),Par) for t in tr] # getting susceptibility
fig = Figure()
ax = Axis(fig[1,1], ylabel = L"χ",xlabel = L"Λ")

scatterlines!(ax,exp.(tr),getindex.(chiR,1))
# scatterlines!(ax,exp.(tr),getindex.(chiR,2))
display("image/png", fig)

######################################################################
######### SQUARE LATTICE ## SQUARE LATTICE ## SQUARE LATTICE #########
######################################################################

using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using SpinFRGLattices.SquareLattice

NLen = 5 # Number of nearest neighbor bonds up to which correlations are treated in the lattice. For NLen = 5, all correlations C_{ij} are zero if sites i and j are separated by more than 5 nearest neighbor bonds.
J1 = 1
J2 = 0.1
couplings = [J1,J2] # Construct a vector of couplings: nearest neighbor coupling is J1 (J2) and further couplings to zero. For finite further couplings simply provide a longer array, i.e [J1,J2,J3,...]

System = getSquareLattice(NLen,couplings) # create a structure that contains all information about the geometry of the problem. 
println(System)

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    T=0.5, # Temperature for the simulation.
    lambda_max = exp(10.),
    lambda_min = exp(-10.),
    N = 8, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-3,
)

@time sol = SolveFRG(Par,method = DP5());

# ## Evaluation Square lattice
# @time begin
    
#     using PMFRGEvaluation
#     using CairoMakie #for plotting. You can use whatever plotting package you like of course

#     System = SquareLattice.getSquareLattice(NLen)
#     Lattice = LatticeInfo(System,SquareLattice)
#     let 
#         chi_R = getChi(sol[end],Par.NumericalParams.lambda_min,Par)
        
#         chi = getFourier(chi_R,Lattice)
        
#         k = LinRange(-2pi,2pi,300)
        
#         chik = [chi(x,y) for x in k, y in k]
        
#         fig, ax, hm = heatmap(k,k,chik,axis = (;aspect = 1))
#         Colorbar(fig[1,2],hm)
#         fig
#     end

# end