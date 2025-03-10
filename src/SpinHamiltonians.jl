struct BkqParam
    Bkq::Dict{Int64, Dict{Int64, ComplexF64}}
end

Base.getindex(bkqparam::BkqParam, k) = bkqparam.Bkq[k]  # enable direct access to the dictionary
Base.keys(bkqparam::BkqParam) = keys(bkqparam.Bkq)  # enable direct access to the dictionary

"""
mult:    Spin multiplicity (=2S+1, where S is the total spin quantum number)
gtensor: Zeeman tensor parametrizing the electronic magnetic moment
Dtensor: Zero-field splitting tensor
Atensors: Hyperfine coupling tensors parametrizing the magnetic field created by the electrons
gammas: Gyromagnetic ratios
"""
struct SHParam
    mult::Int64
    gtensor::Matrix{Float64}
    Bkq::BkqParam
    Atensors::Vector{Matrix{Float64}}
    gammas::Vector{Float64}
end

"""
Constructors with default values for Atensors and gammas: empty lists.
"""
SHParam(mult::Integer, gtensor::Matrix, Dtensor::Matrix) = SHParam(mult, gtensor, Dtensor, Vector{Matrix{Float64}}(), Vector{Float64}())
SHParam(mult::Integer, gtensor::Matrix, Bkq::BkqParam) = SHParam(mult, gtensor, Bkq, Vector{Matrix{Float64}}(), Vector{Float64}())

"""
Constructor taking D-tensor instead of Bkq as argument
"""
function SHParam(mult::Integer, gtensor::Matrix, Dtensor::Matrix, Atensors::Vector{Matrix{Float64}}, gammas::Vector{Float64})
    Bkq = trafo_Dtensor_WybourneBkq(Dtensor)
    return SHParam(mult, gtensor, Bkq, Atensors, gammas)
end

function SHParam_lanthanoid(filename::String, Ln::String, format="Wyb_real")
    # XXXLucasXXX: we should have a single function that can be used
    Bkq = read_Bkq(filename, Ln, format)
    J = ground_J[Ln]
    mult = Int64(2J+1)
    gtensor = Lande_gJ[Ln]*Matrix(1.0I, 3, 3)
    return SHParam(mult, gtensor, Bkq)
end

"""
gtensor: Spin Hamiltonian parameter describing the electronic magnetic moment operator
gammas: Gyromagnetic ratios of the nuclei of interest
R: Nuclear positions relative to the paramagnetic center (metal)

gammas and R must have the same length (number of nuclei for which we calculate the Atensors).
"""
function calc_Atensors_PDA(gtensor::Matrix{Float64}, gammas::Vector{Float64}, R::Vector{Vector{Float64}})
    Nnuc = length(gammas)
    @assert Nnuc == length(R)
    return [alpha^2*gammas[A]/2*calc_dipole_matrix(R[A])*gtensor for A in 1:Nnuc]
end

"""
S::(Pseudo)spin quantum number
H_fieldfree: The field-free Hamiltonian
Mel_trafo:: Matrix applied to vector of base operators to get magnetic moment operators (here: using g-tensor)
BHF_trafo:: Matrix applied to vector of base operators to get hyperfine field operators (here: using HFC A-tensors)
base_op:: vector of base operators (here: spin operators)
"""
struct SpinHamiltonian <: CompModel
    H_fieldfree::HermMat
    Mel_trafo::Matrix{Float64}
    BHF_trafo::Vector{Matrix{Float64}}
    base_op::Vector{Matrix{ComplexF64}}
end

function get_S(sh::SpinHamiltonian)
    mult = size(sh.H_fieldfree)[1]
    S = (mult-1)/2
    return S
end

function calc_Sop(S)
    S_singlemat = calc_soperators(S)   # this is a single array with three dimensions
    return [S_singlemat[:, :, 1], S_singlemat[:, :, 2], S_singlemat[:, :, 3]]
end

"""
Calculate H_fieldfree using Wybourne parametrization.
"""
function calc_H_fieldfree(shparam::SHParam)
    Bkq = shparam.Bkq
    J = (shparam.mult-1)/2
    Tkq = calc_STOs_WE(J)
    H_fieldfree = sum([sum([Bkq[k][q]*Tkq[(k,q)] for q in -k:k]) for k in keys(Bkq)])
    return Hermitian(H_fieldfree)
end

"""
Fallback constructor.
"""
function SpinHamiltonian(args...)
    shparam = SHParam(args...)
    return SpinHamiltonian(shparam)
end

function SpinHamiltonian(shparam::SHParam)
    H_fieldfree = calc_H_fieldfree(shparam)
    S = 0.5*(shparam.mult - 1)    # Spin quantum number
    Sop = calc_Sop(S)
    Mel_trafo = -0.5*shparam.gtensor
    Nnuc = length(shparam.Atensors)
    BHF_trafo = [-(1/shparam.gammas[i])*shparam.Atensors[i] for i in 1:Nnuc]
    return SpinHamiltonian(H_fieldfree, Mel_trafo, BHF_trafo, Sop)
end

function calc_magneticmoment_operator(shparam::SHParam)
    S = 0.5*(shparam.mult - 1)    # Spin quantum number
    Sop = calc_soperators(S)
    Mel = [0.5*sum(shparam.gtensor[i,j] * Sop[:, :, j] for j in 1:3) for i in 1:3]
    return Mel
end

function calc_operators(shparam::SHParam)
    H_fieldfree = calc_H_fieldfree(shparam)
    Mel = calc_magneticmoment_operator(shparam)
    return H_fieldfree, Mel
end

"""
D-tensor has to be provided in atomic units! (not the more common cm-1)
"""
function calc_dyadic(sh::SpinHamiltonian, T::Real)
    solution = eigen(sh.H_fieldfree)
    energies = solution.values
    states = solution.vectors

    SS = calc_F_deriv2(energies, states, sh.base_op, T)
    return SS
end

"""
D-tensor has to be provided in atomic units! (not the more common cm-1)
"""
function calc_tetradic(sh::SpinHamiltonian, T::Real)
    solution = eigen(sh.H_fieldfree)
    energies = solution.values
    states = solution.vectors

    SSSS = calc_F_deriv4(energies, states, sh.base_op, T)
    return SSSS
end

function calc_contactshift_fielddep_Br(s::Float64, Aiso::Matrix{Float64}, g::Matrix{Float64}, D::Matrix{Float64}, T::Real, B0::Float64, gfactor::Float64, direct::Bool=false, selforient::Bool=false)

    gammaI = 2.6752e8*1e-6 
    gammaI *= 2.35051756758e5

    beta = 1/(kB*T)

    SS = calc_dyadics(s, D, T, false)

    #Br = Brillouin(s, T, B0)
    Br = Brillouin_truncated(s, T, B0, gfactor)
    sigma = zeros(length(Aiso), 3, 3)

    chi = pi*alpha^2 * g * SS * g'

    shiftcon = Float64[]

    for (i, Aiso_val) in enumerate(Aiso)

        for l in 1:3, k in 1:3, o in 1:3, p in 1:3
            if k == p
                sigma[i, l, k] += -(1/2) * g[l, o] * (Aiso_val*2pi) *(1/gammaI) * SS[o, p]
            end
        end

        con = 0

        if selforient
            con = ((1/45 * beta/mu0 *tr(sigma[i,:,:])*tr(chi) - 1/15 * beta/mu0 * tr(sigma[i,:,:]*chi))*B0^2) * 1e6
        end

        if direct
            sigma[i,:,:] *= Br
        end

        con += -1/3 * tr(sigma[i, :, :]) * 1e6
        
        push!(shiftcon, con)
    end

    return shiftcon
end

"""
Returns l0, l+1, l-1 (irreducible components of) angular momentum operators with quantum number l
"""
function calc_lops_irreducible(l)
    l0 = calc_sz(l)
    lp1 = -(1/sqrt(2))*calc_splusminus(l,+1)
    lm1 = (1/sqrt(2))*calc_splusminus(l,-1)
    return l0, lp1, lm1
end

"""
Calculates spherical tensor operators J^(k)_q recursively.
"""
function calc_STOs_recursive(l)
    dim = Int64(2l+1)
    T00 = Matrix(1.0I, dim, dim)
    T10, T11, T1m1 = calc_lops_irreducible(l)
    # initialize elements with k=0 and k=1
    T = Dict((0,0) => T00, (1,1) => T11, (1,0) => T10, (1,-1) => T1m1)
    for k in 2:(2l+1), m in -k:k   # k=2l+1 should already be exactly zero; just doing this for testing XXXLucasXXX
        Tkm = zeros(dim, dim)
        for q in -(k-1):(k-1), q_prime in -1:1
            Tkm += clebschgordan(k-1, q, 1, q_prime, k, m)*T[(k-1, q)]*T[(1, q_prime)]
        end
        T[(k, m)] = Tkm
    end
    return T
end

function RME_Lucas(l, k)
    return sqrt((factorial(big(k)))^2 * factorial(big(Int64(2l+k+1))) / Int64(2l+1) / 2^k / factorial(big(2k)) / factorial(big(Int64(2l-k))))
end

function RME_Wybourne(l, k)
    return sqrt(factorial(big(Int64(2l+k+1))) / Int64(2l+1) / factorial(big(Int64(2l-k)))) / 2^k
end

# use Wigner-Eckart theorem to calculate the STOs
# default normalization (as given by the function that calculates the RME)
# is Wybourne (same as Buckmaster and Smith/Thornley)
function calc_STOs_WE(l, RME_func=RME_Wybourne)
    dim = Int64(2l+1)
    offset = l+1
    T = Dict()
    for k in 0:Int64(2l), m in -k:k
        Tkm = zeros(dim, dim)
        RME = RME_func(l, k)
        for M in -l:l, M_prime in -l:l
                Tkm[Int64(offset-M), Int64(offset-M_prime)] = clebschgordan(l, M_prime, k, m, l, M) * RME
        end
        T[(k, m)] = Tkm
    end
    return T
end

"""
Transformation of spherical components (l=0 and l=2) to the Cartesian form
of a symmetric 3x3 tensor.
"""
function symtensor_trafo_sph_Cart(tensor_0_0, tensor_2)
    tensor = Matrix{ComplexF64}(undef, 3, 3)
    tensor[1,1] = 0.5*(tensor_2[2] + tensor_2[-2]) - tensor_2[0]/sqrt(6) - tensor_0_0/sqrt(3)
    tensor[2,2] = -0.5*(tensor_2[2] + tensor_2[-2]) - tensor_2[0]/sqrt(6) - tensor_0_0/sqrt(3)
    tensor[3,3] = 2*tensor_2[0]/sqrt(6) - tensor_0_0/sqrt(3)
    tensor[1,2] = tensor[2,1] = (-im/2)*(tensor_2[2] - tensor_2[-2])
    tensor[1,3] = tensor[3,1] = -0.5*(tensor_2[1] - tensor_2[-1])
    tensor[2,3] = tensor[3,2] = (im/2)*(tensor_2[1] + tensor_2[-1])
    @assert norm(imag(tensor))/norm(real(tensor)) <1e-10
    return real(tensor)
end

"""
Transformation of Cartesian form
of a symmetric 3x3 tensor to spherical components (l=0 and l=2).
"""
function symtensor_trafo_Cart_sph(tensor::Matrix{T}) where T <: Real
    @assert norm(tensor-tensor') < 1e-10
    tensor_0_0 = (-1/sqrt(3))*(tensor[1,1]+tensor[2,2]+tensor[3,3])
    tensor_2 = Dict{Int64,ComplexF64}()
    tensor_2[2] = 0.5*(tensor[1,1]-tensor[2,2] +2*im*tensor[1,2])
    tensor_2[1] = -tensor[1,3]-im*tensor[2,3]
    tensor_2[0] = (2*tensor[3,3]-tensor[1,1]-tensor[2,2])/sqrt(6)
    tensor_2[-1] = tensor[1,3]-im*tensor[2,3]
    tensor_2[-2] = 0.5*(tensor[1,1]-tensor[2,2]-2*im*tensor[1,2])
    return tensor_0_0, tensor_2
end

"""
First derivative of the spin dyadic with respect to beta.
"""
function JJbeta(shparam::SHParam)
    J = (shparam.mult-1)/2
    JJderiv_0_0 = J*(J+1)/sqrt(3)
    JJderiv_2 = Dict(q => 0.0 for q in -2:2)
    return symtensor_trafo_sph_Cart(JJderiv_0_0, JJderiv_2)
end

"""
Second derivative of the spin dyadic with respect to beta.
"""
function JJbeta2(shparam::SHParam)
    J = (shparam.mult-1)/2
    Bkq = shparam.Bkq
    JJderiv_0_0 = 0.0
    # The Wybourne parameters are not proper spherical tensors.
    # Therefore, we have to take the complex conjugate.
    JJderiv_2 = Dict(q => J*(J+1)*(2J+3)*(2J-1)/5/sqrt(6) * conj(Bkq[2][q]) for q in -2:2)
    return symtensor_trafo_sph_Cart(JJderiv_0_0, JJderiv_2)
end

"""
Couple two sets of Wybourne parameters to a new spherical tensor of order K.
"""
function Bk_otimes_Bk(Bkq::BkqParam, k, ktilde, K)
    result = Dict(Q => 0.0*im for Q in -K:K)
    for q in -k:k
        for qtilde in -ktilde:ktilde
            for Q in -K:K
                # The Wybourne ligand field parameters are not proper spherical tensors.
                # Therefore, we have to take the complex conjugate.
                result[Q] += clebschgordan(k, q, ktilde, qtilde, K, Q)*conj(Bkq[k][q])*conj(Bkq[ktilde][qtilde])
            end
        end
    end
    return result
end

"""
Third derivative of the spin dyadic with respect to beta.
"""
function JJbeta3(shparam::SHParam)
    J = (shparam.mult-1)/2
    Bkq = shparam.Bkq
    JJderiv_0_0 = 0.0
    for k in keys(Bkq)
        JJderiv_0_0 += (-1/2/sqrt(3))*k*(k+1)/sqrt(2k+1) * RME_Wybourne(J, k)^2 * Bk_otimes_Bk(Bkq, k, k, 0)[0]
    end
    JJderiv_2 = Dict(Q => 0.0*im for Q in -2:2)
    # first: k=ktilde contribution to anisotropic part of dyadic
    for k in keys(Bkq)
        for Q in -2:2
            JJderiv_2[Q] += sqrt(1/30)*sqrt(k*(k+1)/(2k-1)/(2k+1)/(2k+3))*(6J*(J+1)-2.5*k*(k+1)+3)*RME_Wybourne(J, k)^2 * Bk_otimes_Bk(Bkq, k, k, 2)[Q]
        end
    end
    # second: |k-ktilde|=2 contribution to anisotropic part of dyadic
    for k in keys(Bkq)
        if k==2
            continue
        end
        for Q in -2:2
            JJderiv_2[Q] += -6/sqrt(5) * sqrt(k*(k-1)/(2k+1)/(2k-1)/(2k-3))*RME_Wybourne(J, k)^2 * Bk_otimes_Bk(Bkq, k, k-2, 2)[Q]
        end
    end
    return symtensor_trafo_sph_Cart(JJderiv_0_0, JJderiv_2)
end

function trafo_Dtensor_WybourneBkq(Dtensor::Matrix{T}) where T<:Real
    @assert norm(Dtensor-Dtensor') < 1e-10
    trace = tr(Dtensor)
    if trace>norm(Dtensor)*1e-10
        @warn "Non-traceless D-tensor provided. Only its traceless part will be preserved."
    end
    tensor_0_0, tensor_2 = symtensor_trafo_Cart_sph(Dtensor)   # tensor_0_0 is discarded in the following
    R2 = sqrt(3/2)   # ratio of Wybourne and Koster/Statz normalization for k=2
    B2q = Dict(q => conj(tensor_2[q])/R2 for q in -2:2)
    return BkqParam(Dict(2 => B2q))
end

"""
calc_dyadic_orderx: Calculate angular momentum dyadic truncated at a certain order in beta (1,2, or 3).
"""
function calc_dyadic_order1(shparam::SHParam, T::Real)
    JJderiv1 = JJbeta(shparam)
    beta = 1/kB/T
    return JJderiv1*beta
end

function calc_dyadic_order2(shparam::SHParam, T::Real)
    JJderiv2 = JJbeta2(shparam)
    beta = 1/kB/T
    return calc_dyadic_order1(shparam, T) + JJderiv2*beta^2/2
end

function calc_dyadic_order3(shparam::SHParam, T::Real)
    JJderiv3 = JJbeta3(shparam)
    beta = 1/kB/T
    return calc_dyadic_order2(shparam, T) + JJderiv3*beta^3/6
end

function calc_susceptibility_fromdyadic(dyadic::Matrix{Float64}, gtensor::Matrix{Float64})
    return -pi*alpha^2 * gtensor * dyadic * gtensor'
end