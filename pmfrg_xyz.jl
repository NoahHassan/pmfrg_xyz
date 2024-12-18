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
    NumericalParams::NumericalParams
    Options::OptionParams
end

struct OneLoopWorkspace
    State::StateType    ### Stores the current state
    Deriv::StateType    ### Stores the derivative
    X::BubbleType       ### Stores the bubble function X and XTilde
    Par                 ### Params
end

getVDims(Par) = (Par.System.Npairs, Par.NumericalParams.N, Par.NumericalParams.N, Par.NumericalParams.N)
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
StateType(Par) = StateType(Par.System.NUnique, Par.NumericalParams.N, getVDims(Par), _getFloatType(Par))
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

function V_(Vertex::AbstractArray, ns::Int, nt::Int, nu::Int, Rij::Integer, Rji::Integer, N::Integer)
    ns, nt, nu, swapsites = ConvertFreqArgs(ns, nt, nu, N)
    Rij = ifelse(swapsites, Rji, Rij)
    return Vertex[Rij, ns+1, nt+1, nu+1]
end

function mixedFrequencies(ns,nt,nu,nwpr)
	nw1=Int((ns + nt + nu - 1) / 2)
    nw2=Int((ns - nt - nu - 1) / 2)
    nw3=Int((-ns + nt - nu - 1) / 2)
    nw4=Int((-ns - nt + nu - 1) / 2)

	wpw1 = nwpr + nw1 + 1
    wpw2 = nwpr + nw2 + 1
    wmw3 = nwpr - nw3
    wmw4 = nwpr - nw4

	return wpw1, wpw2, wmw3, wmw4
end

function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
	(; State, X, Par) = Workspace
	N = Par.NumericalParams.N
	(; Npairs, Nsum, siteSum, invpairs) = Par.System

    Va_(Rij,s,t,u) = V_(State.Gamma.a, s, t, u, Rij, invpairs[Rij], N)
	Vb_(Rij,s,t,u) = V_(State.Gamma.b, s, t, u, Rij, invpairs[Rij], N)
	Vc_(Rij,s,t,u) = V_(State.Gamma.c, s, t, u, Rij, invpairs[Rij], N)
	ns = is - 1 ### because julia indexes from 1
	nt = it - 1
	nu = iu - 1
	wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

	# get fields of siteSum struct as Matrices for better use of LoopVectorization
	S_ki = siteSum.ki
	S_kj = siteSum.kj
	S_xk = siteSum.xk
	S_m = siteSum.m

	for Rij in 1:Npairs
		#loop over all left hand side inequivalent pairs Rij
		Xa_sum = 0. #Perform summation on this temp variable before writing to State array as Base.setindex! proved to be a bottleneck!
		Xb_sum = 0.
		Xc_sum = 0.
		for k_spl in 1:Nsum[Rij]
			#loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
			ki,kj,m,xk = S_ki[k_spl,Rij],S_kj[k_spl,Rij],S_m[k_spl,Rij],S_xk[k_spl,Rij]
			Ptm = Props[xk,xk]*m

            Va12 = Va_(ki, ns, wpw1, wpw2)
            Vb12 = Vb_(ki, ns, wpw1, wpw2)
            Vc12 = Vc_(ki, ns, wpw1, wpw2)

            Va34 = Va_(kj, ns, wmw3, wmw4)
            Vb34 = Vb_(kj, ns, wmw3, wmw4)
            Vc34 = Vc_(kj, ns, wmw3, wmw4)

            Vc21 = Vc_(ki, ns, wpw2, wpw1)
            Vc43 = Vc_(kj, ns, wmw4, wmw3)

			Xa_sum += (
				+Va12 * Va34 
				+Vb12 * Vb34 * 2
			)* Ptm

			Xb_sum += (
				+Va12 * Vb34
				+Vb12 * Va34
				+Vb12 * Vb34
			)* Ptm
			
			Xc_sum += (
				+Vc12 * Vc34
				+Vc21 * Vc43
			)* Ptm
		end
		X.a[Rij,is,it,iu] += Xa_sum
		X.b[Rij,is,it,iu] += Xb_sum
		X.c[Rij,is,it,iu] += Xc_sum
    end
    return
end

function addXTilde!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
	(; State, X, Par) = Workspace
	N = Par.NumericalParams.N
	(; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    Va_(Rij,s,t,u) = V_(State.Gamma.a, s, t, u, Rij, invpairs[Rij], N)
	Vb_(Rij,s,t,u) = V_(State.Gamma.b, s, t, u, Rij, invpairs[Rij], N)
	Vc_(Rij,s,t,u) = V_(State.Gamma.c, s, t, u, Rij, invpairs[Rij], N)
	ns = is - 1 ### because julia indexes from 1
	nt = it - 1
	nu = iu - 1
	wpw1, wpw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

	# Xtilde only defined for nonlocal pairs Rij != Rii
	for Rij in 1:Npairs
		Rij in OnsitePairs && continue
		# loop over all left hand side inequivalent pairs Rij
		Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
		(; xi, xj) = PairTypes[Rij]

		# These values are used several times so they are saved locally
		Va12 = Va_(Rji, wpw1, ns, wpw2)
		Va21 = Va_(Rij, wpw2, ns, wpw1)
		Va34 = Va_(Rji, wmw3, ns, wmw4)
		Va43 = Va_(Rij, wmw4, ns, wmw3)

		Vb12 = Vb_(Rji, wpw1, ns, wpw2)
		Vb21 = Vb_(Rij, wpw2, ns, wpw1)
		Vb34 = Vb_(Rji, wmw3, ns, wmw4)
		Vb43 = Vb_(Rij, wmw4, ns, wmw3)

		Vc12 = Vc_(Rji, wpw1, ns, wpw2)
		Vc21 = Vc_(Rij, wpw2, ns, wpw1)
		Vc34 = Vc_(Rji, wmw3, ns, wmw4)
		Vc43 = Vc_(Rij, wmw4, ns, wmw3)

	    X.Ta[Rij,is,it,iu] += (
			(+Va21 * Va43
			+2*Vc21 * Vc43) * Props[xi,xj]
			+(Va12 * Va34
			+2*Vc12 * Vc34)* Props[xj,xi]
		)
		
	    X.Tb[Rij,is,it,iu] += (
			(+Va21 * Vc43
			+Vc21 * Vc43
			+Vc21 * Va43) * Props[xi,xj]

			+(Va12 * Vc34
			+Vc12 * Vc34
			+Vc12 * Va34)* Props[xj,xi]
		)
		Vb12 = Vb_(Rji, wpw1, wpw2, ns)
		Vb21 = Vb_(Rij, wpw2, wpw1, ns)
		Vb34 = Vb_(Rji, wmw3, wmw4, ns)
		Vb43 = Vb_(Rij, wmw4, wmw3, ns)

		Vc12 = Vc_(Rji, wpw1, wpw2, ns)
		Vc21 = Vc_(Rij, wpw2, wpw1, ns)
		Vc34 = Vc_(Rji, wmw3, wmw4, ns)
		Vc43 = Vc_(Rij, wmw4, wmw3, ns)


	    X.Tc[Rij,is,it,iu] += (
			(+Vb21 * Vb43
			+Vc21 * Vc43
			) * Props[xi,xj]
			+(Vb12 * Vb34
			+Vc12 * Vc34
	    	)* Props[xj,xi]
		)
    end
end

function getXBubble!(Workspace, Lam)
	Par = Workspace.Par
    (; T, N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System
	 
	iG(x, nw) = iG_(Workspace.State.iSigma, x, Lam, nw, T)
	iSKat(x,nw) = iSKat_(Workspace.State.iSigma, Workspace.Deriv.iSigma, x, Lam, nw, T)

	function getKataninProp!(BubbleProp,nw1,nw2)
		for i in 1:Par.System.NUnique, j in 1:Par.System.NUnique
			BubbleProp[i,j] = iSKat(i,nw1) * iG(j,nw2)* T
		end
		return SMatrix{NUnique,NUnique}(BubbleProp)
	end
	for is in 1:N, it in 1:N
        BubbleProp = zeros(NUnique,NUnique)
        ns = is - 1
        nt = it - 1
        for nw in -lenIntw:lenIntw-1 # Matsubara sum
            sprop = getKataninProp!(BubbleProp,nw,nw+ns)
            for iu in 1:N
                nu = iu - 1
                if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
                    continue
                end
                addXTilde!(Workspace,is,it,iu,nw,sprop) # add to XTilde-type bubble functions
                if(!Par.Options.use_symmetry || nu<=nt)
                    addX!(Workspace,is,it,iu,nw,sprop)# add to X-type bubble functions
                end
            end
        end
	end
end

function symmetrizeBubble!(X::BubbleType,Par)
    N = Par.NumericalParams.N
    (;Npairs,OnsitePairs) = Par.System
    use_symmetry = Par.Options.use_symmetry
    # use the u <--> t symmetry
    if(use_symmetry)
        for it in 1:N
            for iu in it+1:N, is in 1:N, Rij in 1:Npairs
                X.a[Rij,is,it,iu] = -X.a[Rij,is,iu,it]
                X.b[Rij,is,it,iu] = -X.b[Rij,is,iu,it]
                X.c[Rij,is,it,iu] = (
                + X.a[Rij,is,it,iu]+
                - X.b[Rij,is,it,iu]+
                + X.c[Rij,is,iu,it])
            end
        end
    end
    #local definitions of X.Tilde vertices
    for iu in 1:N
		for it in 1:N, is in 1:N, R in OnsitePairs
			X.Ta[R,is,it,iu] = X.a[R,is,it,iu]
			X.Tb[R,is,it,iu] = X.b[R,is,it,iu]
			X.Tc[R,is,it,iu] = X.c[R,is,it,iu]
			X.Td[R,is,it,iu] = -X.c[R,is,iu,it]
		end
    end
    X.Td .= X.Ta .- X.Tb .- X.Tc
end

function addToVertexFromBubble!(Gamma::VertexType,X::BubbleType)
    for iu in axes(Gamma.a,4)
        for it in axes(Gamma.a,3), is in axes(Gamma.a,2), Rij in axes(Gamma.a,1)
            Gamma.a[Rij,is,it,iu] += X.a[Rij,is,it,iu] - X.Ta[Rij,it,is,iu] + X.Ta[Rij,iu,is,it]
            Gamma.b[Rij,is,it,iu] += X.b[Rij,is,it,iu] - X.Tc[Rij,it,is,iu] + X.Tc[Rij,iu,is,it]
            Gamma.c[Rij,is,it,iu] += X.c[Rij,is,it,iu] - X.Tb[Rij,it,is,iu] + X.Td[Rij,iu,is,it]
        end
    end 
    return Gamma
end

function symmetrizeVertex!(Gamma::VertexType,Par)
	N = Par.NumericalParams.N
	for iu in 1:N
		for it in 1:N, is in 1:N, R in Par.System.OnsitePairs
			Gamma.c[R,is,it,iu] = -Gamma.b[R,it,is,iu]
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
	
	iSigma(x, nw) = iSigma_(State.iSigma, x, nw)
	iG(x, nw) = iG_(State.iSigma, x, Lam, nw, T)
	iS(x, nw) = iS_(State.iSigma, x, Lam, nw, T)

	Theta(Lam, w) = w^2 / (w^2 + Lam^2)
	
	for x in 1:NUnique
		sumres = 0.
		for nw in -lenIntw_acc:lenIntw_acc-1
			w = get_w(nw,T) ### is computed in iS, iG and iSigma too.
			sumres += iS(x, nw) / iG(x, nw) * Theta(Lam, w) * iSigma(x, nw) / w # 1 / w is simply G_0 without Lambda (right?)
		end
		Deriv.f_int[x] = -1.5 * T * sumres
	end
end

function get_Self_Energy!(Workspace, Lam)
	Par = Workspace.Par
	@inline iS(x, nw) = iS_(Workspace.State.iSigma, x, Lam, nw, Par.NumericalParams.T) / 2
	compute1PartBubble!(Workspace.Deriv.iSigma, Workspace.State.Gamma, iS, Par)
end

function compute1PartBubble!(Dgamma::AbstractArray, Gamma::VertexType, Prop, Par)
    invpairs = Par.System.invpairs

	setZero!(Dgamma)

    #### WARUM HIER T, U, S ???????
	@inline Gamma_a(Rij,s,t,u) = V_(Gamma.a, t, u, s, Rij, invpairs[Rij], Par.NumericalParams.N) # Tilde-type can be obtained by permutation of vertices
	@inline Gamma_b(Rij,s,t,u) = V_(Gamma.b, t, u, s, Rij, invpairs[Rij], Par.NumericalParams.N) # cTilde corresponds to b type vertex!

    addTo1PartBubble!(Dgamma, Gamma_a, Gamma_b, Prop, Par)
end

function addTo1PartBubble!(Dgamma::AbstractArray, Gamma_a::Function, Gamma_b::Function, Prop, Par)

    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; siteSum, Nsum, OnsitePairs) = Par.System

	Threads.@threads for iw1 in 1:N
		nw1 = iw1 - 1
    	for (x, Rx) in enumerate(OnsitePairs)
			for nw in -lenIntw_acc:lenIntw_acc-1
				jsum = 0.
				wpw1 = nw1 + nw + 1 #w + w1: Adding two fermionic Matsubara frequencies gives a +1 for the bosonic index
				wmw1 = nw - nw1
				for k_spl in 1:Nsum[Rx]
					(; m, ki, xk) = siteSum[k_spl, Rx]
					jsum += (Gamma_a(ki, wpw1, 0, wmw1) + 2 * Gamma_b(ki, wpw1, 0, wmw1)) * Prop(xk, nw) * m
				end
				Dgamma[x,iw1] += -T * jsum #For the self-energy derivative, the factor of 1/2 must be in the propagator
			end
		end
	end
    # return Dgamma
end

function getDeriv!(Deriv, State, setup, Lam)
    (X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(Deriv, State, X, Par)

    getDFint!(Workspace, Lam)

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

function AllocateSetup(Par::OneLoopParams)
    println("One Loop: T= ",Par.NumericalParams.T)
    ## Allocate Memory:
    X = BubbleType(Par)
    return (X,Par)
end

function InitializeState(Par)

    N = Par.NumericalParams.N;
    ( ; couplings, NUnique) = Par.System;

    VDims = getVDims(Par);
    floattype = _getFloatType(Par)

    State = ArrayPartition(
        zeros(floattype, NUnique),          ### f_int
        zeros(floattype, NUnique, N),       ### iSigma (self energy)
        zeros(floattype, VDims),            ### Va                  This has dimensions (NPairs, n_matsubara, ...),
        zeros(floattype, VDims),            ### Vb                  because V = V(s,t,u) (hence 3 matsubara indices)
        zeros(floattype, VDims)             ### Vc
    );

    Gamma_c = State.x[5]
    setToBareVertex!(Gamma_c, couplings)         ### sets initial conditions to -couplings (see eq. 57)
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

SolveFRG(Par;kwargs...) = launchPMFRG!(InitializeState(Par),AllocateSetup(Par),getDeriv!; kwargs...)

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

function setToBareVertex!(Gamma_c::AbstractArray{T,4}, couplings::AbstractVector) where T
    for Rj in axes(Gamma_c,1)
        Gamma_c[Rj,:,:,:] .= -couplings[Rj]
    end
    return Gamma_c
end

#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################

struct Observables
    Chi
    gamma
end

getChi(State::ArrayPartition, Lam::Real, Par) = getChi(State.x[2], State.x[5], Lam, Par)

function getChi(iSigma::AbstractArray, Gamma_c::AbstractArray, Lam::Real, Par)
	(;T,N,lenIntw_acc) = Par.NumericalParams
	(;Npairs,invpairs,PairTypes,OnsitePairs) = Par.System

	iG(x,w) = iG_(iSigma, x, Lam, w, T)
	Vc_(Rij,s,t,u) = V_(Gamma_c,s,t,u,Rij,invpairs[Rij],N)

	Chi = zeros(_getFloatType(Par),Npairs)

	for Rij in 1:Npairs
		(;xi,xj) = PairTypes[Rij]
		for nK in -lenIntw_acc:lenIntw_acc-1
			if Rij in OnsitePairs
				Chi[Rij,1] += T * iG(xi,nK) ^2
			end
			for nK2 in -lenIntw_acc:lenIntw_acc-1
				npwpw2 = nK+nK2+1
				w2mw = nK2-nK
				#use that Vc_0 is calculated from Vb
				GGGG = iG(xi,nK)^2 * iG(xj,nK2)^2
				Chi[Rij] += T^2 * GGGG *Vc_(Rij,0,npwpw2,w2mw)
			end
        end
    end
	return(Chi)
end


Pkg.add("StructArrays")

##########################################################
######### DIMER SUSC ## DIMER SUSC ## DIMER SUSC #########
##########################################################

using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using CairoMakie

System = getPolymer(2) # create a structure that contains all information about the geometry of the problem. 

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    T = 0.5, # Temperature for the simulation.
    N = 8, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-5,
    lambda_max = 100.,
    lambda_min = .01
)

@time sol = SolveFRG(Par,method = DP5());


## Evaluation Dimer Susc

tr = LinRange(3,-2,20)
chiR = [getChi(sol(t),exp(t),Par) for t in tr] # getting susceptibility
fig = Figure()
ax = Axis(fig[1,1], ylabel = L"χ",xlabel = L"Λ")

scatterlines!(ax,exp.(tr),getindex.(chiR,1))
scatterlines!(ax,exp.(tr),getindex.(chiR,2))
fig

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

## Evaluation Square lattice
@time begin
    
    using PMFRGEvaluation
    using CairoMakie #for plotting. You can use whatever plotting package you like of course

    System = SquareLattice.getSquareLattice(NLen)
    Lattice = LatticeInfo(System,SquareLattice)
    let 
        chi_R = getChi(sol[end],Par.NumericalParams.lambda_min,Par)
        
        chi = getFourier(chi_R,Lattice)
        
        k = LinRange(-2pi,2pi,300)
        
        chik = [chi(x,y) for x in k, y in k]
        
        fig, ax, hm = heatmap(k,k,chik,axis = (;aspect = 1))
        Colorbar(fig[1,2],hm)
        fig
    end

end