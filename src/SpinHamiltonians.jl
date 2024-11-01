"""
gtensor: Zeeman tensor parametrizing the electronic magnetic moment
Atensor: Hyperfine coupling tensor parametrizing the magnetic field created by the electrons
Dtensor: Zero-field splitting tensor
mult:    Spin multiplicity (=2S+1, where S is the total spin quantum number)
"""
struct SpinHamiltonian <: CompModel
    gtensor::Matrix{Float64}
    Atensor::Matrix{Float64}
    Dtensor::Matrix{Float64}
    mult::Int64
end

"""
Implements S^T D S
"""
function calc_H_fieldfree(sh::SpinHamiltonian)
    S = 0.5*(sh.mult - 1)    # Spin quantum number
    Sop = calc_soperators(S)
    StDS = sum(sh.Dtensor[i, j] * Sop[:, :, i] * Sop[:, :, j] for i in 1:3, j in 1:3)
    return StDS
end

function calc_magneticmoment_operator(sh::SpinHamiltonian)
    S = 0.5*(sh.mult - 1)    # Spin quantum number
    Sop = calc_soperators(S)
    Mel = [0.5*sum(sh.gtensor[i,j] * Sop[:, :, j] for j in 1:3) for i in 1:3]
    return Mel
end

function calc_operators(sh::SpinHamiltonian)
    H_fieldfree = calc_H_fieldfree(sh)
    Mel = calc_magneticmoment_operator(sh)
    return H_fieldfree, Mel
end

function calc_dyadics(s::Float64, D::Matrix{Float64}, T::Real, quadruple::Bool)
    S = calc_soperators(s)

    Hderiv = [S[:,:,1], S[:,:,2], S[:,:,3]]

    StDS = sum(D[i, j] * S[:, :, i] * S[:, :, j] for i in 1:3, j in 1:3)

    solution = eigen(StDS)
    energies = solution.values
    states = solution.vectors

    SS = -calc_F_deriv2(energies, states, Hderiv, T)

    if quadruple

        SSSS = -calc_F_deriv4(energies, states, Hderiv, T)

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

# use Wigner-Eckart theorem to calculate the STOs
function calc_STOs_WE(l)
    dim = Int64(2l+1)
    offset = l+1
    T = Dict()
    for k in 0:Int64(2l), m in -k:k
        Tkm = zeros(dim, dim)
        RME = sqrt((factorial(k))^2 * factorial(Int64(2l+k+1)) / Int64(2l+1) / 2^k / factorial(2k) / factorial(Int64(2l-k)))
        for M in -l:l, M_prime in -l:l
                Tkm[Int64(offset-M), Int64(offset-M_prime)] = clebschgordan(l, M_prime, k, m, l, M) * RME
        end
        T[(k, m)] = Tkm
    end
    return T
end
