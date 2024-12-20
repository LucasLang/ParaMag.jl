"""
mult:    Spin multiplicity (=2S+1, where S is the total spin quantum number)
gtensor: Zeeman tensor parametrizing the electronic magnetic moment
Dtensor: Zero-field splitting tensor
Atensors: Hyperfine coupling tensors parametrizing the magnetic field created by the electrons
"""
struct SHParam
    mult::Int64
    gtensor::Matrix{Float64}
    Dtensor::Matrix{Float64}
    Atensors::Vector{Matrix{Float64}}
    gammas::Vector{Float64}
end

"""
Constructor with default value for Atensors: empty list.
"""
SHParam(mult::Integer, gtensor::Matrix, Dtensor::Matrix, Atensors=Vector{Matrix{Float64}}(), gammas=Vector{Float64}()) = SHParam(mult, gtensor, Dtensor, Atensors, gammas)

"""
H_fieldfree: The field-free Hamiltonian
Mel_trafo:: Matrix applied to base operators to get magnetic moment operators (here: using g-tensor)
BHF_trafo:: Matrix applied to base operators to get hyperfine field operators (here: using HFC A-tensors)
base_op:: base operators (here: spin operators)
"""
struct SpinHamiltonian <: CompModel
    H_fieldfree::HermMat
    Mel_trafo::Matrix{Float64}
    BHF_trafo::Vector{Matrix{Float64}}
    base_op::Vector{Matrix{ComplexF64}}
end

function calc_Sop(S)
    S_singlemat = calc_soperators(S)   # this is a single array with three dimensions
    return [S_singlemat[:, :, 1], S_singlemat[:, :, 2], S_singlemat[:, :, 3]]
end

"""
Implements S^T D S
"""
function calc_H_fieldfree(shparam::SHParam)
    S = 0.5*(shparam.mult - 1)    # Spin quantum number
    Sop = calc_soperators(S)
    StDS = sum(shparam.Dtensor[i, j] * Sop[:, :, i] * Sop[:, :, j] for i in 1:3, j in 1:3)
    return Hermitian(StDS)
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
function calc_dyadics(s::Float64, D::Matrix{Float64}, T::Real, quadruple::Bool=false)
    S = calc_soperators(s)

    Hderiv = [S[:,:,1], S[:,:,2], S[:,:,3]]

    StDS = sum(D[i, j] * S[:, :, i] * S[:, :, j] for i in 1:3, j in 1:3)

    solution = eigen(StDS)
    energies = solution.values
    states = solution.vectors

    SS = calc_F_deriv2(energies, states, Hderiv, T)

    if quadruple

        SSSS = calc_F_deriv4(energies, states, Hderiv, T)

        return SS, SSSS
    else
        return SS
    end

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

function calc_contactshift_fielddep(s::Float64, Aiso::Matrix{Float64}, g::Matrix{Float64}, D::Matrix{Float64}, T::Real, B0::Float64, direct::Bool=false, selforient::Bool=false)

    gammaI = 2.6752e8*1e-6 
    gammaI *= 2.35051756758e5

    beta = 1/(kB*T)

    SS, SSSS = calc_dyadics(s, D, T, true)
    
    chi = pi*alpha^2 * g * SS * g'

    sigma1 = zeros(length(Aiso), 3, 3)
    sigma3 = zeros(length(Aiso), 3, 3, 3, 3)

    shiftcon = Float64[]

    for (i, Aiso_val) in enumerate(Aiso)

        for l in 1:3, m in 1:3, n in 1:3, k in 1:3, o in 1:3, p in 1:3, q in 1:3, r in 1:3
            if k == p && m == 1 && n == 1 && q == 1 && r == 1
                sigma1[i, l, k] += -(1/2) * g[l, o] * (Aiso_val*2pi) * (1/gammaI) * SS[o, p]
            end

            if k == r
                sigma3[i, l, m, n, k] += -(1/8) * g[l, o] * g[m, p] * g[n, q] * (Aiso_val*2pi) *(1/gammaI) * SSSS[o, p, q, r]
            end
        end

        con = -1/3 * tr(sigma1[i, :, :]) * 1e6

        if direct
            con += - 1/30 * trace_ord2(sigma3[i, :, :, :, :]) * B0^2 * 1e6
        end

        if selforient
            con += ((1/45 * beta/mu0 *tr(sigma1[i,:,:])*tr(chi) - 1/15 * beta/mu0 * tr(sigma1[i,:,:]*chi))*B0^2) * 1e6
        end

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
Calculate H_fieldfree using Wybourne parametrization.
Later: remove Wyb and use dispatch based on type.
Bkq: dictionary with the ligand field parameters
J: Total angular momentum (or spin, for TM complexes) quantum number
"""
function calc_H_fieldfree_Wyb(Bkq::Dict{Tuple{Int, Int}, ComplexF64}, J)
    Tkq = calc_STOs_WE(J)
    H_fieldfree = sum([Bkq[key]*Tkq[key] for key in keys(Bkq)])
    return Hermitian(H_fieldfree)
end

function calc_H_fieldfree_Wyb(filename::String, Ln::String)
    J = ground_J[Ln]
    Bkq = read_Bkq(filename, Ln)
    return calc_H_fieldfree_Wyb(Bkq, J)
end

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

function calc_dyadics_Wyb(J::Real, Bkq::Dict{Tuple{Int, Int}, ComplexF64}, T::Real)
    Jop = calc_soperators(J)
    Hderiv = [Jop[:,:,1], Jop[:,:,2], Jop[:,:,3]]
    H_fieldfree = calc_H_fieldfree_Wyb(Bkq, J)
    solution = eigen(H_fieldfree)
    energies = solution.values
    states = solution.vectors
    return calc_F_deriv2(energies, states, Hderiv, T)
end

"""
First derivative of the spin dyadic with respect to beta.
"""
function JJbeta(J)
    JJderiv_0_0 = J*(J+1)/sqrt(3)
    JJderiv_2 = Dict(q => 0.0 for q in -2:2)
    return symtensor_trafo_sph_Cart(JJderiv_0_0, JJderiv_2)
end

function JJbeta2(Bkq::Dict{Tuple{Int, Int}, ComplexF64}, J)
    JJderiv_0_0 = 0.0
    # The Wybourne ligand field parameters are not proper spherical tensors.
    # Therefore, we have to take the complex conjugate.
    JJderiv_2 = Dict(q => J*(J+1)*(2J+3)*(2J-1)/5/sqrt(6) * conj(Bkq[(2,q)]) for q in -2:2)
    return symtensor_trafo_sph_Cart(JJderiv_0_0, JJderiv_2)
end