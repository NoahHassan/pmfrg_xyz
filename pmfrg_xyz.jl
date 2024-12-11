### pinned SpinFRGLattices at v0.5.2
import Pkg;
Pkg.add("SpinFRGLattices")

Pkg.add("RecursiveArrayTools")


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

struct VertexType{T}
    a::Array{T, 4}
    b::Array{T, 4}
    c::Array{T, 4}
end

struct BubbleType{T}
    a::Array{T, 4}
    b::Array{T, 4}
    c::Array{T, 4}

    Ta::Array{T, 4}
    Tb::Array{T, 4}
    Tc::Array{T, 4}
    Td::Array{T, 4}
end

struct StateType{T}
    f_int::Array{T}         ### additional index in f, Sigma and Gamma for inequivalent sites (x in geometry package)
    iSigma::Array{T, 2}
    Gamma::VertexType{T}
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
    Params::NumericalParams
    Options::OptionParams
end

struct OneLoopWorkspace
    State::StateType    ### Stores the current state
    Deriv::StateType    ### Stores the derivative
    X::BubbleType       ### Stores the bubble function X and XTilde
    Par                 ### Params
end

getVDims(Par) = (Par.System.Npairs, Par.Params.N, Par.Params.N, Par.Params.N)
_getFloatType(Par) = typeof(Par.NumericalParams.T)

function VertexType(VDims::Tuple)
    return VertexType(
        zeros(VDims),
        zeros(VDims),
        zeros(VDims)
    )
end
VertexType(Par) = VertexType(getVDims(Par))

function BubbleType(VDims::Tuple, type=Float64)
    return BubbleType(
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),

        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims),
        zeros(type, VDims)
    )
end
BubbleType(Par) = BubbleType(getVDims(Par))

function StateType(NUnique::Int, N::Int, VDims::Tuple, type=Float64)
    return StateType(
        zeros(type, NUnique),
        zeros(type, NUnique, N),
        VertexType(VDims)
    )
end
StateType(Par) = StateType(Par.System.NUnique, Par.Params.N, getVDims(Par), _getFloatType(Par))
StateType(f_int, iSigma, Gamma_a, Gamma_b, Gamma_c) = StateType(f_int, iSigma, VertexType(Gamma_a, Gamma_b, Gamma_c))
RecursiveArrayTools.ArrayPartition(x) = ArrayPartition(x.f_int, x.iSigma, x.Gamma.a, x.Gamma.b, x.Gamma.c)
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

OptionParams(;usesymmetry::Bool = true,MinimalOutput::Bool = false,kwargs...) = OptionParams(usesymmetry,MinimalOutput)
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
    return -iG_(iSigma, x, Lam, nw, T)^2 * (2 * Lam / w + iSigma(DSigma, x, nw))
end

####################################################
######### VERTICES ## VERTICES ## VERTICES #########
####################################################