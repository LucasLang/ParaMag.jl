# The functions in this file should not be implemented for specific computational models
# (like LFT or Spin Hamiltonians), but for the CompModel abstract type.
# All functions that are specific to a certain concrete computational model, should go to
# the respective file.

function xyz2spher(x::Real, y::Real, z::Real)
    theta = acos(z)
    phi = atan(y,x)   # 2-argument atan (also known as atan2)
    return theta, phi
end

function setup_Lebedev_grids()
    pkgpath = dirname(pathof(MagFieldLFT))
    gridsizes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]
    lebedev_grids = Vector{Vector{Tuple{Float64, Float64, Float64}}}(undef, 0)
    for N in gridsizes
        xyzgrid = readdlm("$pkgpath/grids/grid_$N")
        grid = Vector{Tuple{Float64, Float64, Float64}}(undef, 0)
        for i in 1:(size(xyzgrid)[1])
            theta, phi = xyz2spher(xyzgrid[i,1], xyzgrid[i,2], xyzgrid[i,3])
            weight = 4pi * xyzgrid[i,4]
            push!(grid, (theta, phi, weight))
        end
        push!(lebedev_grids, grid)
    end
    return lebedev_grids
end

lebedev_grids = setup_Lebedev_grids()

const HermMat = Hermitian{T, Matrix{T}} where T <: Number  # Complex Hermitian (or real symmetric) matrix

"""
Note: this function assumes that the field dependent Hamiltonian is of the form
H(B) = H0 - M*B,
where M is the electronic magnetic moment operator.
This means that Hamiltonians that e.g. include some kind of diamagnetic operator
(quadratic in the external field) are not covered by this function!
"""
function calc_H_magfield(H_fieldfree::HermMat, Mel::Vector{Matrix{ComplexF64}}, B::Vector{Float64})
    H_magfield = H_fieldfree - (Mel[1]*B[1] + Mel[2]*B[2] + Mel[3]*B[3])
    return Hermitian(H_magfield)    # this is important for the eigensolver to yield orthogonal eigenvectors
end

function fieldfree_GS_energy(H_fieldfree::HermMat)
    energies = eigvals(H_fieldfree)
    @assert isa(energies, Vector{Float64})  # energies need to be real
    return energies[1]
end

function calc_free_energy(H_fieldfree::HermMat, Mel::Vector{Matrix{ComplexF64}}, B0_mol::Vector{Float64}, T::Real)
    energies, states = calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    kB = 3.166811563e-6    # Boltzmann constant in Eh/K
    beta = 1/(kB*T)
    energies_exp = exp.(-beta*energies)
    Z = sum(energies_exp)   # canonical partition function
    return -log(Z)/beta
end

"""
H_fieldfree: Hamiltonian in the absence of a magnetic field
H: Hamiltonian whose eigenstates and eigenenergies we want to calculate (relative to fieldfree ground state energy E0)
"""
function calc_solutions(H_fieldfree::HermMat, H::HermMat)
    E0 = fieldfree_GS_energy(H_fieldfree)
    solution = eigen(H)
    @assert isa(solution.values, Vector{Float64})  # energies need to be real
    energies = solution.values .- E0         # All energies relative to fieldfree GS energy
    states = solution.vectors
    return energies, states
end

calc_solutions(H_fieldfree::HermMat) = calc_solutions(H_fieldfree, H_fieldfree)

function calc_solutions_magfield(H_fieldfree::HermMat, Mel::Vector{Matrix{ComplexF64}}, B0_mol::Vector{Float64})
    H_magfield = calc_H_magfield(H_fieldfree, Mel, B0_mol)
    return calc_solutions(H_fieldfree, H_magfield)
end

function calc_Zel(energies::Vector{Float64}, T::Real)
    beta = 1/(kB*T)
    energies_exp = exp.(-beta*energies)
    Z = sum(energies_exp)   # canonical partition function
end

function calc_average_magneticmoment(energies::Vector{Float64}, states::Matrix{ComplexF64}, Mel::Vector{Matrix{ComplexF64}}, T::Real)
    beta = 1/(kB*T)
    energies_exp = exp.(-beta*energies)
    Zel = sum(energies_exp)   # canonical partition function
    Mel_eigenbasis = [states'*Melcomp*states for Melcomp in Mel]
    Mel_avg = [sum(energies_exp .* diag(Mel_eigenbasis_comp))/Zel for Mel_eigenbasis_comp in Mel_eigenbasis]
    @assert norm(imag(Mel_avg)) < 1e-12    # energies need to be real
    return real(Mel_avg)
end

function calc_F_deriv1(energies::Vector{Float64}, states::Matrix{ComplexF64}, Hderiv::Vector{Matrix{ComplexF64}}, T::Real)
    beta = 1/(kB*T)
    energies_exp = exp.(-beta*energies)
    Zel = sum(energies_exp)   # canonical partition function
    Hderiv_eigenbasis = [states'*Hderiv_comp*states for Hderiv_comp in Hderiv]
    Fderiv1 = [sum(energies_exp .* diag(Hderiv_eigenbasis_comp))/Zel for Hderiv_eigenbasis_comp in Hderiv_eigenbasis]
    @assert norm(imag(Fderiv1)) < 1e-12    # need to be real
    return real(Fderiv1)
end

"""
N_polar: Number of grid points for polar angle
N_azimuthal: Number of grid points for azimuthal angle
"""
function spherical_product_grid(N_polar::Integer, N_azimuthal::Integer)
    dtheta = pi/N_polar
    dphi = 2pi/N_azimuthal
    grid = Vector{Tuple{Float64, Float64, Float64}}(undef, 0)
    for theta in (dtheta/2):dtheta:(pi-dtheta/2)
        for phi in (dphi/2):dphi:(2pi-dphi/2)
            weight = sin(theta)*dtheta*dphi
            push!(grid, (theta, phi, weight))
        end
    end
    return grid
end

"""
The function f is assumed to return a list of values (vector-valued).
f is a function of points on the unit sphere, parametrized via a polar angle theta and an azimuthal angle phi.
"""
function integrate_spherical(f::Function, grid::Vector{Tuple{Float64, Float64, Float64}})
    dim = length(f(0,0))    # number of components in the output
    integrals = zeros(dim)
    for (theta, phi, weight) in grid
        integrals += weight*f(theta, phi)
    end
    return integrals
end

function calc_dipole_matrix(R::Vector{T}) where T<:Real
    R_length = norm(R)
    idmat = Matrix(1.0I, 3, 3)
    return ((3*R*R')/(R_length^2) - idmat)/(R_length^3)
end

"""
Calculate magnetic field created by a magnetic dipole moment (everything in atomic units).
"""
function dipole_field(M::Vector{T1}, R::Vector{T2}) where {T1<:Real, T2<:Real}
    return alpha^2 * (calc_dipole_matrix(R)*M)
end

"""
B0: field strength (scalar quantity).
The B0 vector in the laboratory frame is (0,0,B0)
"""
function calc_B0_mol(B0::Real, chi::Real, theta::Real)
    Bx = -sin(theta)*cos(chi)*B0
    By = sin(theta)*sin(chi)*B0
    Bz = cos(theta)*B0
    return [Bx, By, Bz]
end

function calc_Bind_avg_lab_z(chi::Real, theta::Real, Bind_avg_mol::Vector{Float64})
    return -sin(theta)*cos(chi)*Bind_avg_mol[1] +
            sin(theta)*sin(chi)*Bind_avg_mol[2] +
            cos(theta)*         Bind_avg_mol[3]
end

"""
theta: Euler angle: Angle between z axis (molecular frame) and Z axis (lab frame)
chi: Euler angle describing rotations of the molecule around its molecular frame z axis
H_fieldfree: Hamiltonian in the absence of a magnetic field (Slater determinant basis)
L: Orbital angular momentum operators (Slater determinant basis)
S: Total spin operators (Slater determinant basis)
Mel: Electronic magnetic dipole moment operators (Slater determinant basis)
R: Vectors from the points at which we want to know the induced field (typically nuclear positions) to the paramagnetic center (atomic units = Bohr)
B0: Magnitude of the external magnetic field (atomic units)
T: Temperature (Kelvin)
"""
function calc_integrands(theta::Real, chi::Real, H_fieldfree::HermMat, Mel::Vector{Matrix{ComplexF64}}, R::Vector{Vector{Float64}}, B0::Real, T::Real)
    B0_mol = calc_B0_mol(B0, chi, theta)
    energies, states = calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Zel = calc_Zel(energies, T)
    Mel_avg_mol = calc_average_magneticmoment(energies, states, Mel, T)
    Bind_avg_mol = [dipole_field(Mel_avg_mol, R_i) for R_i in R]
    Bind_avg_lab_z = [calc_Bind_avg_lab_z(chi, theta, B_i) for B_i in Bind_avg_mol]
    return [Zel*Bind_avg_lab_z; Zel]
end

# XXXLucasXXX: Generalize this function such that it only depends on magnetic moment operators, not on L, S, or ligand field parameters
#              Goal: Also use this together with spin Hamiltonian
"""
R: Vectors from the points at which we want to know the induced field (typically nuclear positions) to the paramagnetic center (atomic units = Bohr)
B0: Magnitude of the external magnetic field (atomic units)
T: Temperature (Kelvin)
"""
function calc_Bind(param::LFTParam, R::Vector{Vector{Float64}}, B0::Real, T::Real, grid::Vector{Tuple{Float64, Float64, Float64}})
    H_fieldfree, Mel = calc_operators_SDbasis(param)
    integrands(theta, chi) = calc_integrands(theta, chi, H_fieldfree, Mel, R, B0, T)
    integrals = integrate_spherical(integrands, grid)
    numerators = integrals[1:(end-1)]
    D = integrals[end]   # denominator (normalization)
    return [N/D for N in numerators]
end

function estimate_shifts_finitefield(param::LFTParam, R::Vector{Vector{Float64}}, B0::Real, T::Real, grid::Vector{Tuple{Float64, Float64, Float64}})
    Bind_values = calc_Bind(param, R, B0, T, grid)
    return (Bind_values / B0) * 1e6    # convert to ppm
end

struct DegenerateSet
    E::Float64              # energy of the states belonging to the set
    states::Vector{Int64}   # indices of the states belonging to the set
end

function determine_degenerate_sets(energies::Vector{Float64})
    degenerate_sets = Vector{DegenerateSet}(undef, 0)
    current_energy = energies[1]
    current_states = [1]
    for i in 2:length(energies)
        if abs(energies[i]-current_energy) < 1e-10                  # threshold is arbitrarily chosen. Need to think more about this
            push!(current_states, i)
        else
            push!(degenerate_sets, DegenerateSet(current_energy, current_states))
            current_energy = energies[i]
            current_states = [i]
        end
    end
    push!(degenerate_sets, DegenerateSet(current_energy, current_states))
    return degenerate_sets
end

"""
Determine all indices that are distinct from the given M_indices.
"""
function calc_N_indices(number_indices::Real, M_indices::Vector{Int64})
    return [i for i in 1:number_indices if !(i in M_indices)]
end

function calc_Hderiv_part(Hderiv_eigenbasis::Vector{Matrix{ComplexF64}}, I_indices::Vector{Int64}, J_indices::Vector{Int64})
    return [Hderivcomp[I_indices, J_indices] for Hderivcomp in Hderiv_eigenbasis]
end

function calc_Hderiv_eigenbasis_denom(Hderiv_eigenbasis::Vector{Matrix{ComplexF64}}, denoms::Vector{Float64})
    Hderiv_eigenbasis_denom = deepcopy(Hderiv_eigenbasis)
    for n in 1:length(denoms)
        for i in 1:3
            Hderiv_eigenbasis_denom[i][:, n] /= denoms[n]
        end
    end
    return Hderiv_eigenbasis_denom
end

function calc_F_deriv2(energies::Vector{Float64}, states::Matrix{ComplexF64}, Hderiv::Vector{Matrix{ComplexF64}}, T::Real)
    beta = 1/(kB*T)
    Zel = calc_Zel(energies, T)
    Hderiv_eigenbasis = [states'*Hderivcomp*states for Hderivcomp in Hderiv]
    degenerate_sets = determine_degenerate_sets(energies)
    Fderiv2 = im*zeros(3, 3)
    for M in degenerate_sets
        M_indices = M.states
        N_indices = calc_N_indices(length(energies), M_indices)
        Boltzmann_factor = exp(-beta*M.E)
        Hderiv_MM = calc_Hderiv_part(Hderiv_eigenbasis, M_indices, M_indices)
        Hderiv_MN = calc_Hderiv_part(Hderiv_eigenbasis, M_indices, N_indices)
        Hderiv_eigenbasis_denom = calc_Hderiv_eigenbasis_denom(Hderiv_eigenbasis, energies .- M.E)
        Hderiv_MN_denom = calc_Hderiv_part(Hderiv_eigenbasis_denom, M_indices, N_indices)
        for k in 1:3
            for l in 1:3
                part1 = 0.5*beta*tr(Hderiv_MM[k] * Hderiv_MM[l]')
                part2 = tr(Hderiv_MN[k] * Hderiv_MN_denom[l]')
                Fderiv2[k, l] += Boltzmann_factor * (part1 + part2)
            end
        end
    end
    Fderiv2 *= -1/Zel
    Fderiv1 = calc_F_deriv1(energies, states, Hderiv, T)
    Fderiv2 += 0.5*beta* Fderiv1*Fderiv1'
    Fderiv2 += transpose(Fderiv2)  # symmetrization
    @assert norm(imag(Fderiv2)) < 1e-5
    return real(Fderiv2)
end

function calc_F_deriv3(energies::Vector{Float64}, states::Matrix{ComplexF64}, Hderiv::Vector{Matrix{ComplexF64}}, T::Real)
    beta = 1/(kB*T)
    Zel = calc_Zel(energies, T)
    Hderiv_eigenbasis = [states'*Hderivcomp*states for Hderivcomp in Hderiv]
    degenerate_sets = determine_degenerate_sets(energies)
    Fderiv3 = im*zeros(3, 3, 3)
    for M in degenerate_sets
        M_indices = M.states
        N_indices = [i for i in 1:length(energies) if !(i in M_indices)]
        Boltzmann_factor = exp(-beta*M.E)
        Hderiv_MM = calc_Hderiv_part(Hderiv_eigenbasis, M_indices, M_indices)
        Hderiv_MN = calc_Hderiv_part(Hderiv_eigenbasis, M_indices, N_indices)
        Hderiv_NN = calc_Hderiv_part(Hderiv_eigenbasis, N_indices, N_indices)
        Hderiv_eigenbasis_denom = calc_Hderiv_eigenbasis_denom(Hderiv_eigenbasis, energies .- M.E)
        Hderiv_MN_denom = calc_Hderiv_part(Hderiv_eigenbasis_denom, M_indices, N_indices)
        for l in 1:3
            for k in 1:3
                for j in 1:3
                    part1 = (1/6)*beta^2*tr(Hderiv_MM[l] * Hderiv_MM[k] * Hderiv_MM[j])
                    part2 = beta*tr(Hderiv_MN_denom[l] * Hderiv_MN[k]' * Hderiv_MM[j])
                    part3 = - tr(Hderiv_MN_denom[l] * Hderiv_MN_denom[k]' * Hderiv_MM[j])
                    part4 = tr(Hderiv_MN_denom[l] * Hderiv_NN[k] * Hderiv_MN_denom[j]')
                    Fderiv3[l,k,j] += Boltzmann_factor * (part1 + part2 + part3 + part4)
                end
            end
        end
    end
    Fderiv3 *= 1/Zel
    Fderiv1 = calc_F_deriv1(energies, states, Hderiv, T)
    Fderiv2 = calc_F_deriv2(energies, states, Hderiv, T)
    for l in 1:3
        for k in 1:3
            for j in 1:3
                Fderiv3[l,k,j] += -(1/6)*beta^2*Fderiv1[l]*Fderiv1[k]*Fderiv1[j]
                Fderiv3[l,k,j] += 0.5*beta*Fderiv1[l]*Fderiv2[k,j]
            end
        end
    end
    # symmetrization:
    nindices = 3
    Fderiv3_symmetrized = im*zeros(3,3,3)
    for k in 1:factorial(nindices)   # loop over all permutations of three indices
        Fderiv3_symmetrized += permutedims(Fderiv3, Permutation(nindices,k))
    end
    @assert norm(imag(Fderiv3_symmetrized)) < 1e-4
    return real(Fderiv3_symmetrized)
end

function calc_F_deriv4(energies::Vector{Float64}, states::Matrix{ComplexF64}, Hderiv::Vector{Matrix{ComplexF64}}, T::Real)
    beta = 1/(kB*T)
    Zel = calc_Zel(energies, T)
    Hderiv_eigenbasis = [states'*Hderivcomp*states for Hderivcomp in Hderiv]
    degenerate_sets = determine_degenerate_sets(energies)
    Fderiv4 = im*zeros(3, 3, 3, 3)
    for M in degenerate_sets
        M_indices = M.states
        N_indices = [i for i in 1:length(energies) if !(i in M_indices)]
        Boltzmann_factor = exp(-beta*M.E)
        Hderiv_MM = calc_Hderiv_part(Hderiv_eigenbasis, M_indices, M_indices)
        Hderiv_MN = calc_Hderiv_part(Hderiv_eigenbasis, M_indices, N_indices)
        Hderiv_NN = calc_Hderiv_part(Hderiv_eigenbasis, N_indices, N_indices)
        Hderiv_eigenbasis_denom = calc_Hderiv_eigenbasis_denom(Hderiv_eigenbasis, energies .- M.E)
        Hderiv_eigenbasis_denom2 = calc_Hderiv_eigenbasis_denom(Hderiv_eigenbasis, (energies .- M.E) .* (energies .- M.E))
        Hderiv_MN_denom = calc_Hderiv_part(Hderiv_eigenbasis_denom, M_indices, N_indices)
        Hderiv_MN_denom2 = calc_Hderiv_part(Hderiv_eigenbasis_denom2, M_indices, N_indices)
        Hderiv_NN_denom = calc_Hderiv_part(Hderiv_eigenbasis_denom, N_indices, N_indices)
        for l in 1:3, k in 1:3, j in 1:3, i in 1:3
            part1 = (1/24)*beta^3*tr(Hderiv_MM[l] * Hderiv_MM[k] * Hderiv_MM[j] * Hderiv_MM[i])
            part2 = (1/2)*beta^2*tr(Hderiv_MN_denom[l] * Hderiv_MN[k]' * Hderiv_MM[j] * Hderiv_MM[i])
            part3 = -beta*tr(Hderiv_MN_denom2[l] * Hderiv_MN[k]' * Hderiv_MM[j] * Hderiv_MM[i])
            part4 = tr(Hderiv_MN_denom2[l] * Hderiv_MN_denom[k]' * Hderiv_MM[j] * Hderiv_MM[i])
            part5 = 0.5*beta*tr(Hderiv_MN_denom[l] * Hderiv_MN[k]' * Hderiv_MN_denom[j] * Hderiv_MN[i]')
            part6 = beta*tr(Hderiv_MN_denom[l] * Hderiv_NN[k] * Hderiv_MN_denom[j]' * Hderiv_MM[i])
            part7 = -tr(Hderiv_MN_denom2[l] * Hderiv_MN[k]' * Hderiv_MN_denom[j] * Hderiv_MN[i]')
            part8 = -tr(Hderiv_MN_denom[l] * Hderiv_NN[k] * Hderiv_MN_denom2[j]' * Hderiv_MM[i])
            part9 = part8'
            part10 = tr(Hderiv_MN_denom[l] * Hderiv_NN_denom[k] * Hderiv_NN_denom[j] * Hderiv_MN[i]')
            Fderiv4[l,k,j,i] += Boltzmann_factor * (part1 + part2 + part3 + part4 + part5 +
                                                    part6 + part7 + part8 + part9 + part10)
        end
    end
    Fderiv4 *= -1/Zel
    Fderiv1 = calc_F_deriv1(energies, states, Hderiv, T)
    Fderiv2 = calc_F_deriv2(energies, states, Hderiv, T)
    Fderiv3 = calc_F_deriv3(energies, states, Hderiv, T)
    for l in 1:3, k in 1:3, j in 1:3, i in 1:3
        Fderiv4[l,k,j,i] += (1/24)*beta^3*Fderiv1[l]*Fderiv1[k]*Fderiv1[j]*Fderiv1[i] -
                        (1/4)*beta^2*Fderiv1[l]*Fderiv1[k]*Fderiv2[j,i] +
                        (1/6)*beta*Fderiv1[l]*Fderiv3[k,j,i] +
                        (1/8)*beta*Fderiv2[l,k]*Fderiv2[j,i]
    end
    # symmetrization:
    nindices = 4
    Fderiv4_symmetrized = im*zeros(3,3,3,3)
    for k in 1:factorial(nindices)   # loop over all permutations of three indices
        Fderiv4_symmetrized += permutedims(Fderiv4, Permutation(nindices,k))
    end
    @assert norm(imag(Fderiv4_symmetrized))/norm(real(Fderiv4_symmetrized)) < 1e-8    #1e-10
    return real(Fderiv4_symmetrized)
end


# XXXLucasXXX: The following two functions should be modified such that they can be used with any model Hamiltonian
#              (not just full multiplet LFT)
function calc_susceptibility_vanVleck(param::LFTParam, T::Real)
    B = [0.0, 0.0, 0.0]
    Fderiv2 = calc_F_deriv2(param, T, B)
    return -4pi*alpha^2 * Fderiv2
end

"""
Returns the chemical shifts in ppm calculated according to the Kurland-McGarvey equation (point-dipole approximation)
R: Vectors from the points at which we want to know the induced field (typically nuclear positions) to the paramagnetic center (atomic units = Bohr)
T: Temperature (Kelvin)
"""
function calc_shifts_KurlandMcGarvey(param::LFTParam, R::Vector{Vector{Float64}}, T::Real)
    chi = calc_susceptibility_vanVleck(param, T)
    shifts = Vector{Float64}(undef, 0)
    for Ri in R
        D = calc_dipole_matrix(Ri)
        sigma = -(1/(4pi)) * chi * D
        shift = -(1/3)*tr(sigma)
        push!(shifts, shift)
    end
    shifts *= 1e6    # convert to ppm
    return shifts
end

function Brillouin(S::Float64, T::Float64, B0::Float64)

    muB = 0.5
    ge = 2.0

    k_B = 3/2*kB*T /(ge*muB*S*(S+1)*B0)
    a = muB*ge*B0/(2*kB*T)
    if B0 != 0
        Br = k_B*((2*S+1)/tanh((2*S+1)*a) - 1/tanh(a))
    else
        Br = 1.0
    end

    return Br
end

function Brillouin_truncated(S::Float64, T::Float64, B0::Float64, gfactor::Float64=2.0)
    Br = 1 + B0^2*0.5^2*gfactor^2/(240*(kB*T)^2*S*(S+1)) * (1-(2*S+1)^4)
    return Br
end

function orientation_tensor(B0::Float64, T::Float64, chi::Array{Float64, 2})
    #field-induced self-orientation tensor 
    #(eq 112 from Parigi, G. et al, Progress in Nuclear Magnetic Resonance Spectroscopy 114–115 (2019) 211–236) 

    w,v = eigen(chi)
    a = (B0^2)/(5*mu0*kB*T)
    chiiso = (1/3)*tr(Diagonal(w))
    Pw = zeros(3,3)
    for i in 1:3
    	Pw[i,i]= (1 + a*(w[i]-chiiso))
    end
    P = v * Pw * v'
    return P
end

function trace_ord2(tensor::Array{Float64, 4})
    # computes the order 2 trace of a supersymmetric fourth order tensor

    trace = @tensor begin    # XXXLucasLangXXX: trace should be deleted either in this or the next line
        trace = tensor[i, i, j, j]
    end

    return trace
end

function product_ord3(tensor::Array{Float64, 4}, Dip::Array{Float64, 2})
    #dot product between two fourth order tensors

    sigma = zeros(Float64, 3, 3, 3, 3)

    @tensor begin
        sigma[l, m, n, k] := tensor[l, m, n, q] * Dip[q, k]
    end

    return sigma
end

function calc_shifts_KurlandMcGarvey_ord4(chi::Array{Float64, 2}, chi3::Array{Float64, 4}, R::Vector{Vector{Float64}}, T::Real, B0::Real, direct::Bool=false, indirect::Bool=false)
    #pcs calculation with Kurland-McGarvey equation 
    #the saturation effect is accounted for with fourth order tensor (determined via analytical equation)

    beta = 1/(kB*T)
  
    shifts = Vector{Float64}(undef, 0)
    for Ri in R
        Dip = calc_dipole_matrix(Ri)
        sigma = -(1/(4pi)) * chi * Dip
        shift = -(1/3)*tr(sigma)

        if direct
            tau = -(1/(4pi)) * (1/6) * product_ord3(chi3,Dip)
            shift += - (1/5)*trace_ord2(tau)*B0^2
        end

        if indirect
            shift += (1/45 * beta/mu0 *tr(sigma)*tr(chi) - 1/15 * beta/mu0 * tr(sigma*chi))*B0^2
        end

        push!(shifts, shift)
    end
    shifts *= 1e6    # convert to ppm
    return shifts
end

function calc_shifts_KurlandMcGarvey_Br(chi::Array{Float64, 2}, R::Vector{Vector{Float64}}, T::Real, B0::Float64, S::Float64, gfactor::Float64, direct::Bool=false, indirect::Bool=false)
    #pcs calculation with Kurland-McGarvey equation
    #the saturation effect is accounted for with Brillouin equation

    beta = 1/(kB*T)

    #Br = Brillouin(S, T, B0)
    Br = Brillouin_truncated(S, T, B0, gfactor)

    shifts = Vector{Float64}(undef, 0)
    for Ri in R
        Dip = calc_dipole_matrix(Ri)
        sigma = -(1/(4pi)) * chi * Dip
        shift = -(1/3)*tr(sigma)

        if direct
            shift = -(1/3)*tr(sigma .* Br)
        end

        if indirect
            shift += (1/45 * beta/mu0 *tr(sigma)*tr(chi) - 1/15 * beta/mu0 * tr(sigma*chi))*B0^2 
        end

        push!(shifts, shift)
    end
    shifts *= 1e6    # convert to ppm
    return shifts
end
