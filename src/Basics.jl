"""
Calculate lz operator in basis of complex atomic orbitals.
"""
calc_lz(l::Int) = diagm(l:-1:-l)

function calc_lplusminus(l::Int, sign::Int)
    @assert Int64(abs(sign)) == 1       # sign may only be +1 or -1
    dim = 2l+1
    mvalues = l:-1:-l
    op = zeros(dim,dim)
    for i_prime in 1:dim
        m_prime = mvalues[i_prime]
        for i in 1:dim
            m = mvalues[i]
            if m_prime == (m + sign)
                op[i_prime, i] = sqrt(l*(l+1) - m*(m+sign))
            end
        end
    end
    return op
end

"""
Calculate sz operator in basis of complex atomic orbitals.
"""
calc_sz(l::Float64) = diagm(l:-1:-l)

function calc_splusminus(l::Float64, sign::Int)
    @assert Int64(abs(sign)) == 1       # sign may only be +1 or -1
    dim = Int(2l+1)
    mvalues = l:-1:-l
    op = zeros(dim,dim)
    for i_prime in 1:dim
        m_prime = mvalues[i_prime]
        for i in 1:dim
            m = mvalues[i]
            if m_prime == (m + sign)
                op[i_prime, i] = sqrt(l*(l+1) - m*(m+sign))
            end
        end
    end
    return op
end

"""
Returns lz, lplus, lminus in basis of complex atomic orbitals of angular momentum l.
"""
function calc_lops_complex(l::Int)
    return calc_lz(l), calc_lplusminus(l,+1), calc_lplusminus(l,-1)
end

function Tesla2au(B::Real)
    return B/2.35051756758e5
end

"""
nu: Resonance frequency in MHz
"""
function MHz2Tesla(nu::Real)
    gamma = 42.577478518     # gyromagnetic ratio of the proton in MHz/T
    return nu/gamma
end

"""
nu: Resonance frequency in MHz
"""
function MHz2au(nu::Real)
    return Tesla2au(MHz2Tesla(nu))
end
