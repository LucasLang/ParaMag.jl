using Test, LinearAlgebra

using MagFieldLFT

function test_createSDs()
    ref = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]]
    return MagFieldLFT.create_SDs(3,2) == ref
end

function test_createSDs2()
    nel = 1
    norb = 2
    ref = [[1], [2], [3], [4]]
    return MagFieldLFT.create_SDs(nel, norb) == ref
end

function test_U_complex2real()
    ref = [1/sqrt(2) 0 0 0 -im/sqrt(2);
           0 -1/sqrt(2) 0 im/sqrt(2) 0;
           0 0 1 0 0;
           0 1/sqrt(2) 0 im/sqrt(2) 0;
           1/sqrt(2) 0 0 0 im/sqrt(2)]
    return MagFieldLFT.U_complex2real(2) ≈ ref
end

function test_calc_lops()
    l=2
    ref_lz = [2 0 0 0 0;
              0 1 0 0 0;
              0 0 0 0 0;
              0 0 0 -1 0;
              0 0 0 0 -2]
    ref_lplus = [0 2 0 0 0;
                 0 0 sqrt(6) 0 0;
                 0 0 0 sqrt(6) 0;
                 0 0 0 0 2;
                 0 0 0 0 0]
    ref_lminus = [0 0 0 0 0;
                  2 0 0 0 0;
                  0 sqrt(6) 0 0 0;
                  0 0 sqrt(6) 0 0;
                  0 0 0 2 0]
    lz, lplus, lminus = MagFieldLFT.calc_lops_complex(l)
    return lz ≈ ref_lz && lplus ≈ ref_lplus && lminus ≈ ref_lminus
end

function test_calcERIs_complex()
    F = Dict(0=>3.2, 2=>(1/7/7)*1.72, 4=>(1/21/21)*2.20)
    l = 2
    ERIs = MagFieldLFT.calcERIs_complex(l, F)
    return ERIs[2,1,4,5] ≈ (-2.12/9) && ERIs[3,2,5,5] ≈ 0.0 && ERIs[2,2,5,5] ≈ (195.92/63)
end

function test_calcERIs_real()
    F = Dict(0=>3.2, 2=>(1/7/7)*1.72, 4=>(1/21/21)*2.20)
    l = 2
    ERIs = MagFieldLFT.calcERIs_real(l, F)
    return ERIs[2,2,5,5] ≈ (195.92/63) && ERIs[1,4,2,5] ≈ (-0.64/21)
end

function test_ERIs_symmetries()
    F = Dict(0=>3.2, 2=>1.72, 4=>2.20)
    l = 2
    ERIs = MagFieldLFT.calcERIs_real(l, F)
    @assert abs(ERIs[1,2,4,5]) > 1e-4
    t1 = ERIs[1,2,4,5] ≈ ERIs[2,1,4,5]
    t2 = ERIs[1,2,4,5] ≈ ERIs[1,2,5,4]
    t3 = ERIs[1,2,4,5] ≈ ERIs[2,1,5,4]
    t4 = ERIs[1,2,4,5] ≈ ERIs[4,5,1,2]
    t5 = ERIs[1,2,4,5] ≈ ERIs[5,4,1,2]
    t6 = ERIs[1,2,4,5] ≈ ERIs[4,5,2,1]
    t7 = ERIs[1,2,4,5] ≈ ERIs[5,4,2,1]
    return t1 && t2 && t3 && t4 && t5 && t6 && t7
end

function test_spinorb2orbindex()
    return MagFieldLFT.spinorb2orbindex(2) == (1, 'β') && MagFieldLFT.spinorb2orbindex(5) == (3, 'α')
end

function test_occ_list()
    SD = [1,3,4,7]
    occ_alpha_list = MagFieldLFT.occ_list(SD, 'α')
    occ_beta_list = MagFieldLFT.occ_list(SD, 'β')
    return occ_alpha_list == [(1, 1),(2, 2),(4, 4)] && occ_beta_list == [(3,2)]
end

function test_orb2spinorb_and_back()
    test1 = MagFieldLFT.spinorb2orbindex(MagFieldLFT.orb2spinorbindex(3,'α')) == (3,'α')
    test2 = MagFieldLFT.orb2spinorbindex(MagFieldLFT.spinorb2orbindex(17)...) == 17
    return test1 && test2
end

function test_unocc_list()
    SD = [1,2,4,5,6,11]
    norb = 7
    unocc_alpha_list = MagFieldLFT.unocc_list(SD, norb, 'α')
    unocc_beta_list = MagFieldLFT.unocc_list(SD, norb, 'β')
    return unocc_alpha_list == [2,4,5,7] && unocc_beta_list == [4,5,6,7]
end

function test_SD2index()
    N = 3
    norb = 7
    M = 2norb
    SDs = MagFieldLFT.create_SDs(N,norb)
    test1 = MagFieldLFT.SD2index(SDs[51], M) == 51
    test2 = MagFieldLFT.SD2index(SDs[32], M) == 32
    return test1 && test2
end

function test_Z_summand()
    M = 14
    N = 3
    i = 2
    P = 3
    return MagFieldLFT.Z_summand(i,P,N,M) == 11
end

function test_calc_exc_equal()
    SD = [1,2]
    norb = 3
    sigma_p = 'α'
    sigma_q = 'α'
    exc = MagFieldLFT.calc_exc_occ2unocc(SD, norb, sigma_p, sigma_q)
    return exc == [(6,2,1,-1), (8,3,1,-1)]
end

function test_calc_exc_minus()
    SD = [3,4]
    norb = 3
    sigma_p = 'β'
    sigma_q = 'α'
    exc = MagFieldLFT.calc_exc_occ2unocc(SD, norb, sigma_p, sigma_q)
    return exc == [(7,1,2,1), (14,3,2,-1)]
end

function test_calc_exc_occ2self()
    SD = [1,2,4]
    norb = 5
    exc_alpha_self = MagFieldLFT.calc_exc_occ2self(SD, norb, 'α')
    exc_beta_self = MagFieldLFT.calc_exc_occ2self(SD, norb,'β')
    test_alpha = exc_alpha_self == [(2,1,1,1)]
    test_beta = exc_beta_self == [(2,1,1,1), (2,2,2,1)]
    return test_alpha && test_beta
end

function test_calc_exclists()
    l = 1
    N = 2
    exc = MagFieldLFT.calc_exclists(l,N)
    test1 = exc.alpha[1] == [(1,1,1,1), (6,2,1,-1), (8,3,1,-1)]
    test2 = exc.minus[10] == [(7,1,2,1), (14,3,2,-1)]
    return test1 && test2
end

"""
For 2 electrons in a shell of p orbitals, there are the following terms with their energies:
E(3P) = F0 - 5*F2
E(1D) = F0 + F2
E(1S) = F0 + 10*F2
(see Griffith (the theory of transition-metal ions) Chapter 4.5)
"""
function test_calc_H_nonrel1()
    l=1
    N=2
    exc = MagFieldLFT.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    F = Dict(0 => 1.0, 2 => 2.0)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = MagFieldLFT.calc_H_nonrel(param, exc)
    E = eigvals(H)
    return E[1]≈-9 && E[9]≈-9 && E[10]≈3 && E[14]≈3 && E[15]≈21
end

"""
For 2 electrons in a shell of d orbitals, there are the following terms with their energies:
E(3F) = A - 8*B
E(1D) = A - 3*B + 2*C
E(3P) = A + 7*B
E(1G) = A + 4*B + 2*C
E(1S) = A + 14*B + 7*C
(see Griffith (the theory of transition-metal ions) Table 4.6)
"""
function test_calc_H_nonrel2()
    l=2
    N=2
    exc = MagFieldLFT.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = MagFieldLFT.Racah2F(A,B,C)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = MagFieldLFT.calc_H_nonrel(param, exc)
    E = eigvals(H)
    return E[1]≈-15 && E[22]≈3 && E[27]≈15 && E[36]≈17 && E[45]≈57
end

function test_calc_H_fieldfree()
    l=2
    N=2
    exc = MagFieldLFT.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = MagFieldLFT.Racah2F(A,B,C)
    zeta = 0.5
    param = LFTParam(N, norb, l, hLFT, F, zeta)
    H = MagFieldLFT.calc_H_fieldfree(param, exc)
    E = eigvals(H)
    return E[1]≈E[5] && E[6]≈E[12] && E[13]≈E[21] && !(E[1]≈E[6]) && !(E[6]≈E[13])
end

"""
This test checks for the total orbital angular momentum expectation value L(L+1)
for the different terms.
"""
function test_calc_L()
    l=2
    N=2
    exc = MagFieldLFT.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = MagFieldLFT.Racah2F(A,B,C)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = MagFieldLFT.calc_H_nonrel(param, exc)
    C = eigvecs(H)
    Lx, Ly, Lz = MagFieldLFT.calc_L(l, exc)
    L2 = Lx*Lx + Ly*Ly + Lz*Lz
    L2val = diag(C'*L2*C)
    return L2val[1]≈12 && L2val[22]≈6 && L2val[27]≈2 && L2val[36]≈20 && (L2val[45]+1)≈1
end

"""
This test checks for the spin orbital angular momentum expectation value S(S+1)
for the different terms.
"""
function test_calc_S()
    l=2
    N=2
    exc = MagFieldLFT.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = MagFieldLFT.Racah2F(A,B,C)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = MagFieldLFT.calc_H_nonrel(param, exc)
    C = eigvecs(H)
    Sx, Sy, Sz = MagFieldLFT.calc_S(l, exc)
    S2 = Sx*Sx + Sy*Sy + Sz*Sz
    S2val = diag(C'*S2*C)
    return S2val[1]≈2 && (S2val[22]+1)≈1 && S2val[27]≈2 && (S2val[36]+1)≈1 && (S2val[45]+1)≈1
end

"""
This test checks for the total angular momentum expectation value J(J+1)
for the lowest three spin-orbit-coupled terms (originating from 3F term).
"""
function test_total_J()
    l=2
    N=2
    exc = MagFieldLFT.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = MagFieldLFT.Racah2F(A,B,C)
    zeta = 0.5
    param = LFTParam(N, norb, l, hLFT, F, zeta)
    H = MagFieldLFT.calc_H_fieldfree(param, exc)
    C = eigvecs(H)
    Lx, Ly, Lz = MagFieldLFT.calc_L(l, exc)
    Sx, Sy, Sz = MagFieldLFT.calc_S(l, exc)
    Jx = Lx+Sx
    Jy = Ly+Sy
    Jz = Lz+Sz
    J2 = Jx*Jx + Jy*Jy + Jz*Jz
    J2val = diag(C'*J2*C)
    return J2val[1]≈6 && J2val[6]≈12 && J2val[13]≈20
end

function test_read_AILFT_params_ORCA()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H = MagFieldLFT.calc_H_nonrel(param, exc)
    E = eigvals(H)
    E = (E .- E[1])*27.211    # take ground state energy as reference and convert to eV
    # println(round(E[5], digits=3))
    # println(round(E[13], digits=3))
    # println(round(E[17], digits=3))
    # println(round(E[end], digits=3))
    test1 = round(E[5], digits=3) == 1.638    # first excited quartet state
    test2 = round(E[13], digits=3) == 1.645   # third excited quartet state
    test3 = round(E[17], digits=3) == 2.398   # lowest doublet state
    test4 = round(E[end], digits=3) == 9.930  # highest doublet state
    return test1 && test2 && test3 && test4
end

# because of precision issues in the printed LFT parameters, the energies do not coincide exactly
# with what is printed in the ORCA output file!
# TO DO: change later.
function test_Ercomplex()
    param = read_AILFT_params_ORCA("ErCl63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H = MagFieldLFT.calc_H_nonrel(param, exc)
    E = eigvals(H)
    E = (E .- E[1])*27.211  # take ground state energy as reference and convert to eV
    test1 = round(real(E[29]), digits=3) == 0.020
    test2 = round(real(E[53]), digits=3) == 2.194
    test3 = round(real(E[75]), digits=3) == 2.238
    return test1 && test2 && test3
end

function test_Ercomplex_SOC()
    param = read_AILFT_params_ORCA("ErCl63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H = MagFieldLFT.calc_H_fieldfree(param, exc)
    E = eigvals(H)
    E = (E .- E[1])*219474.63  # take ground state energy as reference and convert to cm-1
    test1 = round(real(E[17]), digits=3) ==   5862.094
    test2 = round(real(E[43]), digits=3) ==   12235.462
    test3 = round(real(E[364]), digits=3) ==  118616.544
    return test1 && test2 && test3
end

function test_calc_free_energy()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H_fieldfree = MagFieldLFT.calc_H_fieldfree(param, exc)
    S = MagFieldLFT.calc_S(param.l, exc)
    L = MagFieldLFT.calc_L(param.l, exc)
    Mel = MagFieldLFT.calc_magneticmoment_operator(L, S)
    B = [0.0, 0.0, 1.0e-5]
    T = 298.0
    F1 = MagFieldLFT.calc_free_energy(H_fieldfree, Mel, B, T)
    return F1 ≈ -0.0013082816934478216
end

function test_average_magnetic_moment()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H_fieldfree = MagFieldLFT.calc_H_fieldfree(param, exc)
    S = MagFieldLFT.calc_S(param.l, exc)
    L = MagFieldLFT.calc_L(param.l, exc)
    Mel = MagFieldLFT.calc_magneticmoment_operator(L, S)
    B0_mol = [0.0, 0.0, 0.0]
    T = 298.0
    energies, states = MagFieldLFT.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Hderiv = -Mel
    Fderiv1 = MagFieldLFT.calc_F_deriv1(energies, states, Hderiv, T)
    Mel_avg = -Fderiv1
    return Mel_avg + [1.0, 1.0, 1.0] ≈ [1.0, 1.0, 1.0]    # magnetization is zero in absence of field
end

# At low temperature, magnetization should be that of the ground state (approximately MS=-3/2)
function test_average_magnetic_moment2()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H_fieldfree = MagFieldLFT.calc_H_fieldfree(param, exc)
    S = MagFieldLFT.calc_S(param.l, exc)
    L = MagFieldLFT.calc_L(param.l, exc)
    Mel = MagFieldLFT.calc_magneticmoment_operator(L, S)
    B0_mol = [0.0, 0.0, 1.0e-4]
    T = 1.0
    energies, states = MagFieldLFT.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Mel = MagFieldLFT.calc_magneticmoment_operator(L,S)
    Mel_avg = MagFieldLFT.calc_average_magneticmoment(energies, states, Mel, T)
    return Mel_avg ≈ [-0.0003645214756898332, -1.2322563787476262e-13, 1.4631881898029349]
end

# At low temperature, magnetization should be that of the ground state (approximately MS=-1/2,
# having magnetization of <-1/2| Mz | -1/2> = - <-1/2|Sz|-1/2> = +1/2)
function test_average_magnetic_moment3()
    nel = 9
    norb = 5
    l = 2
    hLFT = diagm([0.3, 0.05, 0.0, 0.05, 0.1])   # energies of x2-y2, xz, z2, yz, xy
    F = Dict(0 => 0.0, 2 => 0.0, 4 => 0.0)    # does not matter for d9 system
    zeta = 0.0
    param = MagFieldLFT.LFTParam(nel, norb, l, hLFT, F, zeta)

    T = 0.0001
    B0_mol = [0, 0, 1.0e-7]

    H_fieldfree, Mel = MagFieldLFT.calc_operators(param)
    energies, states = MagFieldLFT.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Mel_avg_finitefield = MagFieldLFT.calc_average_magneticmoment(energies, states, Mel, T)

    return norm(Mel_avg_finitefield - [0.0, 0.0, 0.5]) < 1.0e-4
end

function test_integrate_spherical()
    Y_1_m1(theta,phi) = 0.5*sqrt(3/(2pi)) * exp(-im*phi)*sin(theta)
    Y_1_0(theta,phi) = 0.5*sqrt(3/pi) * cos(theta)
    grid = MagFieldLFT.spherical_product_grid(50,50)
    # Results of integration: inner products [<Y_1_m1|Y_1_m1>, <Y_1_m1|Y_1_0>, <Y_1_0|Y_1_0>]
    f(x,y) = [Y_1_m1(x,y)'*Y_1_m1(x,y), Y_1_m1(x,y)'*Y_1_0(x,y), Y_1_0(x,y)'*Y_1_0(x,y)]
    integrals = MagFieldLFT.integrate_spherical(f, grid)
    integrals = [round(x, digits=3) for x in integrals]
    return integrals ≈ [1.000, 0.000, 1.000]
end

function test_integrate_spherical_Lebedev()
    Y_1_m1(theta,phi) = 0.5*sqrt(3/(2pi)) * exp(-im*phi)*sin(theta)
    Y_1_0(theta,phi) = 0.5*sqrt(3/pi) * cos(theta)
    grid = lebedev_grids[25]
    # Results of integration: inner products [<Y_1_m1|Y_1_m1>, <Y_1_m1|Y_1_0>, <Y_1_0|Y_1_0>]
    f(x,y) = [Y_1_m1(x,y)'*Y_1_m1(x,y), Y_1_m1(x,y)'*Y_1_0(x,y), Y_1_0(x,y)'*Y_1_0(x,y)]
    integrals = MagFieldLFT.integrate_spherical(f, grid)
    integrals = [round(x, digits=10) for x in integrals]
    return integrals ≈ [1.000, 0.000, 1.000]
end

function test_dipole_matrix()
    R = [1,5,-2.0]
    dipmat = MagFieldLFT.calc_dipole_matrix(R)
    ref = [-0.005477225575051661 0.003042903097250923 -0.0012171612389003691;
    0.003042903097250923 0.009128709291752768 -0.006085806194501846;
    -0.0012171612389003691 -0.006085806194501846 -0.0036514837167011074]
    return dipmat ≈ ref
end

function test_dipole_field()
    R = [5,7,9]
    M = [-4.0, 3, 7.5]
    B_dip = MagFieldLFT.dipole_field(M, R)
    ref = [2.9330997775211083e-7, 1.7331548609510163e-7, 1.2230892547235214e-7]
    return B_dip ≈ ref
end

function test_determine_degenerate_sets()
    energies = [0.0, 1e-12, 1e-11, 2.2, 2.2+1e-11, 2.2+1e-9, 7, 9, 9, 9]
    degenerate_sets = MagFieldLFT.determine_degenerate_sets(energies)
    D = MagFieldLFT.DegenerateSet
    ref = [D(0.0, [1,2,3]), D(2.2, [4,5]), D(2.2+1e-9, [6]), D(7.0, [7]), D(9.0, [8,9,10])]
    passed = true
    for i in 1:length(ref)
        passed = passed && (degenerate_sets[i].E == ref[i].E)
        passed = passed && (degenerate_sets[i].states == ref[i].states)
    end
    return passed
end

function test_determine_degenerate_sets2()
    energies = zeros(5)
    degenerate_sets = MagFieldLFT.determine_degenerate_sets(energies)
    D = MagFieldLFT.DegenerateSet
    ref = [D(0.0, [1,2,3,4,5])]
    passed = true
    for i in 1:length(ref)
        passed = passed && (degenerate_sets[i].E == ref[i].E)
        passed = passed && (degenerate_sets[i].states == ref[i].states)
    end
    return passed
end

# in weak-field / high-temperature limit, finite-field magnetization should be linear in external field
function test_calc_susceptibility_vanVleck()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 298
    B0_mol = [0, 0, 1.0e-7]

    # version 1
    chi = MagFieldLFT.calc_susceptibility_vanVleck(lft, T)
    Mel_avg_linear = (1/(4pi*MagFieldLFT.alpha^2))*chi*B0_mol

    # version 2
    H_fieldfree, Mel = MagFieldLFT.calc_operators(param)
    energies, states = MagFieldLFT.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Mel_avg_finitefield = MagFieldLFT.calc_average_magneticmoment(energies, states, Mel, T)

    return norm(Mel_avg_finitefield - Mel_avg_linear) < 1.0e-10
end

bohrinangstrom = 0.529177210903
# atom counting starting from 1 (total number of atoms is 49, NH proton is the last one)
r_Ni = [0.000,   0.000,   0.000]                      # atom 33
r_NH = [0.511,  -2.518,  -0.002]                      # atom 49
r_CH1 = [1.053,   1.540,   3.541]                     # atom 23
r_CH2 = [-0.961,  -1.048,  -3.741]                    # atom 32
r_alpha1_alpha2prime_1 = [-1.500,  -3.452,   1.130]   # atom 44
r_alpha1_alpha2prime_2 = [0.430,  -3.104,   2.402]    # atom 45
R_NH                   = r_Ni - r_NH
R_CH1                  = r_Ni - r_CH1
R_CH2                  = r_Ni - r_CH2
R_alpha1_alpha2prime_1 = r_Ni - r_alpha1_alpha2prime_1
R_alpha1_alpha2prime_2 = r_Ni - r_alpha1_alpha2prime_2
R_selected_NiSAL = [R_NH, R_CH1, R_CH2, R_alpha1_alpha2prime_1, R_alpha1_alpha2prime_2] / bohrinangstrom

function test_KurlandMcGarvey()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 298   # I actually did not find in the paper at which temperature they recorded it!?
    shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey(lft, R_selected_NiSAL, T)
    ref = [-96.1951957001265, 30.53047905687625, 30.175351679457314, -22.804411834183288, -19.50459482031643]
    # calculated PCS at CASSCF/NEVPT2/QDPT level according to SI of paper:
    # [-61.5, 21.7, 21.3, -15.5, -13.6]
    # My shifts are larger in magnitude by around 50%, but the relative size and sign is correct
    # This could well be an artifact of only using CASSCF for determining the AILFT parameters (and/or of the LFT approximation)
    return shifts ≈ ref
end

function test_KurlandMcGarvey_vs_finitefield()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 298   # I actually did not find in the paper at which temperature they recorded it!?
    KMcG_shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey(lft, R_selected_NiSAL, T)
    grid = MagFieldLFT.spherical_product_grid(100,100)
    B0 = 1.0e-6
    finitefield_shifts = MagFieldLFT.estimate_shifts_finitefield(lft, R_selected_NiSAL, B0, T, grid)
    return norm(KMcG_shifts - finitefield_shifts) < 0.1
end

function test_KurlandMcGarvey_vs_finitefield_Lebedev()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 298   # I actually did not find in the paper at which temperature they recorded it!?
    KMcG_shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey(lft, R_selected_NiSAL, T)
    grid = lebedev_grids[20]
    B0 = 1.0e-7
    finitefield_shifts = MagFieldLFT.estimate_shifts_finitefield(lft, R_selected_NiSAL, B0, T, grid)
    return norm(KMcG_shifts - finitefield_shifts) < 1.0e-6
end

function test_Fderiv2_numeric_vs_analytic()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    B0_mol = [0.0, 0.0, 1.0e-4]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv1_0 = MagFieldLFT.calc_F_deriv1(lft, T, B0_mol)
    Fderiv1_x = MagFieldLFT.calc_F_deriv1(lft, T, B0_mol+[h, 0, 0])
    Fderiv1_y = MagFieldLFT.calc_F_deriv1(lft, T, B0_mol+[0, h, 0])
    Fderiv1_z = MagFieldLFT.calc_F_deriv1(lft, T, B0_mol+[0, 0, h])
    Fderiv2 = MagFieldLFT.calc_F_deriv2(lft, T, B0_mol)
    Fderiv2_numeric = zeros(3,3)
    Fderiv2_numeric[1, :] = (1/h)*(Fderiv1_x - Fderiv1_0)
    Fderiv2_numeric[2, :] = (1/h)*(Fderiv1_y - Fderiv1_0)
    Fderiv2_numeric[3, :] = (1/h)*(Fderiv1_z - Fderiv1_0)
    return norm(Fderiv2-Fderiv2_numeric) < 0.1
end

function test_Fderiv3_numeric_vs_analytic()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    B0_mol = [0.0, 0.0, 1.0e-4]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv2_0 = MagFieldLFT.calc_F_deriv2(lft, T, B0_mol)
    Fderiv2_x = MagFieldLFT.calc_F_deriv2(lft, T, B0_mol+[h, 0, 0])
    Fderiv2_y = MagFieldLFT.calc_F_deriv2(lft, T, B0_mol+[0, h, 0])
    Fderiv2_z = MagFieldLFT.calc_F_deriv2(lft, T, B0_mol+[0, 0, h])
    Fderiv3 = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol)
    Fderiv3_numeric = zeros(3,3,3)
    Fderiv3_numeric[1, :, :] = (1/h)*(Fderiv2_x - Fderiv2_0)
    Fderiv3_numeric[2, :, :] = (1/h)*(Fderiv2_y - Fderiv2_0)
    Fderiv3_numeric[3, :, :] = (1/h)*(Fderiv2_z - Fderiv2_0)
    # note: elements of the tensor have magnitudes that are all larger than 1e5
    return norm(Fderiv3-Fderiv3_numeric) < 8000
end

function test_Fderiv4_numeric_vs_analytic()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    B0_mol = [0.0, 0.0, 1.0e-4]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv3_0 = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol)
    Fderiv3_x = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol+[h, 0, 0])
    Fderiv3_y = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol+[0, h, 0])
    Fderiv3_z = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol+[0, 0, h])
    Fderiv4 = MagFieldLFT.calc_F_deriv4(lft, T, B0_mol)
    Fderiv4_numeric = zeros(3,3,3,3)
    Fderiv4_numeric[1, :, :, :] = (1/h)*(Fderiv3_x - Fderiv3_0)
    Fderiv4_numeric[2, :, :, :] = (1/h)*(Fderiv3_y - Fderiv3_0)
    Fderiv4_numeric[3, :, :, :] = (1/h)*(Fderiv3_z - Fderiv3_0)
    return rel_diff_norm(Fderiv4, Fderiv4_numeric) < 1e-4
end

function test_Fderiv4_numeric_vs_analytic_zerofield()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    B0_mol = [0.0, 0.0, 0.0]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv3_0 = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol)
    Fderiv3_x = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol+[h, 0, 0])
    Fderiv3_y = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol+[0, h, 0])
    Fderiv3_z = MagFieldLFT.calc_F_deriv3(lft, T, B0_mol+[0, 0, h])
    Fderiv4 = MagFieldLFT.calc_F_deriv4(lft, T, B0_mol)
    Fderiv4_numeric = zeros(3,3,3,3)
    Fderiv4_numeric[1, :, :, :] = (1/h)*(Fderiv3_x - Fderiv3_0)
    Fderiv4_numeric[2, :, :, :] = (1/h)*(Fderiv3_y - Fderiv3_0)
    Fderiv4_numeric[3, :, :, :] = (1/h)*(Fderiv3_z - Fderiv3_0)
    return rel_diff_norm(Fderiv4, Fderiv4_numeric) < 1e-4
end

function rel_diff_norm(value, ref)
    return norm(value-ref)/norm(ref)
end

function test_print_composition()
    C = [1,2,3]
    C = C/norm(C) # normalize
    U = [1 0 0; 0 1 0; 0 0 1]
    labels = ["ex", "ey", "ez"]
    thresh = 0.98
    buf = IOBuffer()
    MagFieldLFT.print_composition(C, U, labels, thresh, buf)
    printed_string = String(take!(buf))
    ref = """
     64.29%  ez
     28.57%  ey
      7.14%  ex
    """
    return printed_string == ref
end

function test_group_eigenvalues()
    values = [1,1,2,5,7,7,7,10]
    unique_values, indices = MagFieldLFT.group_eigenvalues(values)
    ref_values = [1,2,5,7,10]
    ref_indices = [[1,2], [3], [4], [5,6,7], [8]]
    return (unique_values == ref_values) && (indices == ref_indices)
end

function test_print_composition2()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H_fieldfree = MagFieldLFT.calc_H_fieldfree(param, exc)
    energies_rel, states_rel = MagFieldLFT.calc_solutions(H_fieldfree)
    C_rel_ground = states_rel[:,1]
    C_rel_first = states_rel[:,2]
    C_list, labels_list = MagFieldLFT.adapt_basis_Hnonrel_S2_Sz(param)
    thresh = 0.99
    buf = IOBuffer()
    MagFieldLFT.print_composition(C_rel_ground, C_list, labels_list, thresh, buf)
    printed_string_ground = String(take!(buf))
    MagFieldLFT.print_composition(C_rel_first, C_list, labels_list, thresh, buf)
    printed_string_first = String(take!(buf))

    ref_ground = """
     45.52%  E = -30.339968, S =  1.0, M_S = -1.0
     45.52%  E = -30.339968, S =  1.0, M_S =  1.0
      5.65%  E = -30.339968, S =  1.0, M_S =  0.0
      0.87%  E = -30.323346, S =  1.0, M_S =  0.0
      0.84%  E = -30.323346, S =  1.0, M_S = -1.0
      0.84%  E = -30.323346, S =  1.0, M_S =  1.0
    """
    ref_first = """
     39.81%  E = -30.339968, S =  1.0, M_S =  0.0
     28.52%  E = -30.339968, S =  1.0, M_S =  1.0
     28.52%  E = -30.339968, S =  1.0, M_S = -1.0
      1.26%  E = -30.323346, S =  1.0, M_S =  1.0
      1.26%  E = -30.323346, S =  1.0, M_S = -1.0
    """
    return (ref_ground == printed_string_ground) && (ref_first == printed_string_first)
end

function test_print_composition_ErIII()
    param = read_AILFT_params_ORCA("ErCl63-.out", "CASSCF")
    exc = MagFieldLFT.calc_exclists(param)
    H_fieldfree = MagFieldLFT.calc_H_fieldfree(param, exc)
    energies_rel, states_rel = MagFieldLFT.calc_solutions(H_fieldfree)
    C_rel_ground = states_rel[:,1]
    C_list, labels_list = MagFieldLFT.adapt_basis_L2_S2_J2_Jz(param)
    C_list_alt, labels_list_alt = MagFieldLFT.adapt_basis_L2_S2_J2_Jz(param, "term_symbols")
    thresh = 0.99
    buf = IOBuffer()
    MagFieldLFT.print_composition(C_rel_ground, C_list, labels_list, thresh, buf)
    printed_string_ground = String(take!(buf))
    MagFieldLFT.print_composition(C_rel_ground, C_list_alt, labels_list_alt, thresh, buf)
    printed_string_ground_alt = String(take!(buf))

    ref_ground = """
     35.74%  (L, S, J, M_J) = ( 6.0, 1.5, 7.5, 2.5)
     32.82%  (L, S, J, M_J) = ( 6.0, 1.5, 7.5, 6.5)
     19.61%  (L, S, J, M_J) = ( 6.0, 1.5, 7.5,-1.5)
      9.77%  (L, S, J, M_J) = ( 6.0, 1.5, 7.5,-5.5)
      0.69%  (L, S, J, M_J) = ( 7.0, 0.5, 7.5, 2.5)
      0.63%  (L, S, J, M_J) = ( 7.0, 0.5, 7.5, 6.5)
    """
    ref_ground_alt = """
     35.74%  Term: 4I15/2, M_J = 5/2
     32.82%  Term: 4I15/2, M_J = 13/2
     19.61%  Term: 4I15/2, M_J = -3/2
      9.77%  Term: 4I15/2, M_J = -11/2
      0.69%  Term: 2K15/2, M_J = 5/2
      0.63%  Term: 2K15/2, M_J = 13/2
    """
    return (ref_ground == printed_string_ground) && (ref_ground_alt == printed_string_ground_alt)
end

#at very low magnetic field the field dependent results and the zero field ones should be the same 
function test_KurlandMcGarvey_ord4_Br_field()

    bohrinangstrom = 0.529177210903
    # atom counting starting from 1 (total number of atoms is 49, NH proton is the last one)
    r_Ni = [0.000,   0.000,   0.000]                      # atom 33
    r_NH = [0.511,  -2.518,  -0.002]                      # atom 49  NH
    r_CH1 = [1.053,   1.540,   3.541]                     # atom 23  CH'
    r_CH2 = [-0.961,  -1.048,  -3.741]                    # atom 32  CH
    r_alpha1_alpha2prime_1 = [-1.500,  -3.452,   1.130]   # atom 44  alpha1
    r_alpha1_alpha2prime_2 = [0.430,  -3.104,   2.402]    # atom 45  alpha2'
    R_NH                   = r_Ni - r_NH
    R_CH1                  = r_Ni - r_CH1
    R_CH2                  = r_Ni - r_CH2
    R_alpha1_alpha2prime_1 = r_Ni - r_alpha1_alpha2prime_1
    R_alpha1_alpha2prime_2 = r_Ni - r_alpha1_alpha2prime_2
    R = [R_NH, R_CH1, R_CH2, R_alpha1_alpha2prime_1, R_alpha1_alpha2prime_2] / bohrinangstrom

    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 298. 
    B0_MHz = 10. 
    B0 = B0_MHz/42.577478518/2.35051756758e5
    S = 1.0

    shift_ord4 = MagFieldLFT.calc_shifts_KurlandMcGarvey_ord4(lft, R, T, B0, true, true)
    KMcG_shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey(lft, R, T)
    Br_shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey_Br(lft, R, T, B0, S, 2.0, true, true)

    return norm(KMcG_shifts - Br_shifts) < 1.0e-4  && norm(KMcG_shifts - shift_ord4) < 1.0e-4
end

#at very high temperature the difference in field dependent effects simulation should decrease
#due to the population of all the levels of the ground multiplet  
function test_KurlandMcGarvey_ord4_Br_temperature()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 400. 
    B0_MHz = 400. 
    B0 = B0_MHz/42.577478518/2.35051756758e5
    S = 1.0

    shift_ord4 = MagFieldLFT.calc_shifts_KurlandMcGarvey_ord4(lft, R_selected_NiSAL, T, B0, true, true)
    Br_shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey_Br(lft, R_selected_NiSAL, T, B0, S, 2.0, true, true)

    return norm((shift_ord4 - Br_shifts) .* 100 ./ shift_ord4) < 5.0e-2
end


function test_calc_contactshift()

    from_au = 27.2113834*8065.54477
    #from casscf-nevpt2 effective hamiltonian
    D = [3.021254 -5.482999 -16.364810; -5.482999 21.938653 -7.107412; -16.364810 -7.107412 -0.931482]
    D /= from_au  #from cm-1 in au

    #from casscf-nevpt2 effective hamiltonian
    g = [2.320356 0.036880 0.109631; 0.036810 2.180382 0.053365; 0.109499 0.053077 2.347002]
    #from DFT calc
    Aiso = [-0.2398	0.3677	-0.0953	0.1169	5.6311	0.8206	2.5466	0.9669	2.2237	-0.2058	0.3801	-0.0323	0.0943	5.7383	-0.1162	-0.1048	1.3015	4.0604	3.9516	0.929	-0.1677	-0.2015	-3.5469]
    T = 298.0

    B0_MHz = 400. 
    B0 = B0_MHz/42.577478518/2.35051756758e5

    shiftcon = MagFieldLFT.calc_contactshift_fielddep(1.0, Aiso, g, D, T, 0.0, false, false)
    
    indices = [1, 23, 5, 14]
    calc_shift = [shiftcon[i] for i in indices]

    #values from high-temperature limit equation (from paper Inorg. Chem. 2021, 60, 3, 2068–2075):
    # check = [-19.4159467515336, -265.947444074357, 455.934686207511, 464.614375497604]
    # the difference is due to D contribution 
    check = [-19.32946495357911, -285.90358316868117, 453.9038786492883, 462.54490718566717]

    return norm(check-calc_shift) < 1e-6

end

function test_KurlandMcGarvey_vs_finitefield_Lebedev_ord4()

    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = MagFieldLFT.LFT(param)
    T = 298.
    bohrinangstrom = 0.529177210903

    coordinate_matrice = []
    nomi_atomi = []
    coordinate_H = []

    open("nisalhdpt.xyz", "r") do file
        for ln in eachline(file)
            clmn = split(ln)
            push!(nomi_atomi, clmn[1])
            push!(coordinate_matrice, parse.(Float64, clmn[2:end]))
            if clmn[1]=="H"
                push!(coordinate_H, parse.(Float64, clmn[2:end]))
            end
        end
    end

    R = convert(Vector{Vector{Float64}}, coordinate_H) ./ bohrinangstrom

    shift_0 = MagFieldLFT.calc_shifts_KurlandMcGarvey_ord4(lft, R, T, 0.0, false, false)
    grid = lebedev_grids[20]
    B0_single = [10.]
    diff_list_ord4 = []
    diff_list_finitefield = []
    for B0_MHz in B0_single 
        B0 = B0_MHz/42.577478518/2.35051756758e5
        finitefield_shifts = MagFieldLFT.estimate_shifts_finitefield(lft, R, B0, T, grid)
        shift_ord4 = MagFieldLFT.calc_shifts_KurlandMcGarvey_ord4(lft, R, T, B0, true, true)
        push!(diff_list_ord4, shift_0 .- shift_ord4)
        push!(diff_list_finitefield, shift_0 .- finitefield_shifts)
    end

    return norm(diff_list_ord4[1]-diff_list_finitefield[1])<1e-7

end

function test_cubicresponse_spin()
    T = 298.0
    beta = 1/(MagFieldLFT.kB*T)
    S = 2
    Sp = MagFieldLFT.calc_lplusminus(S, +1)
    Sm = MagFieldLFT.calc_lplusminus(S, -1)
    Sz = MagFieldLFT.calc_lz(S)

    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5im * (Sp-Sm)
    Svec = [Sx, Sy, Sz]

    dim = 2S+1
    energies = zeros(dim)
    states = im*Matrix(1.0I, dim, dim)
    SiSjSkSl = MagFieldLFT.calc_F_deriv4(energies, states, Svec, 298.0)
    ref = zeros(3,3,3,3)
    d = Matrix(1.0I, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        ref[i,j,k,l] = (beta^3/45)*S*(S+1)*(2S^2+2S+1)*(d[i,j]*d[k,l] + d[i,k]*d[j,l] + d[i,l]*d[j,k])
    end
    return norm(SiSjSkSl - ref)/norm(ref) < 1.0e-10
end

function test_STOs()
    TWE = MagFieldLFT.calc_STOs_WE(3.0, MagFieldLFT.RME_Lucas)
    Trec = MagFieldLFT.calc_STOs_recursive(3.0)
    return norm(Trec[(4,-3)]-TWE[(4,-3)]) < 1e-10 && norm(Trec[(5,2)]-TWE[(5,2)]) < 1e-10
end

function test_PCS_PDA_finitefield_SH()
    mult = 3   # NiSAL has a triplet ground state

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*MagFieldLFT.cmm1_Hartree   # directly convert from cm-1 to Hartree

    gtensor = [2.1384111    0.0084976    0.0250646;
    0.0074791    2.0934328    0.0112682;
    0.0228213    0.0119502    2.1324169]

    sh = MagFieldLFT.SpinHamiltonian(mult, gtensor, Dtensor)

    T = 298.0
    grid = lebedev_grids[20]
    B0_MHz = 10.0
    B0 = B0_MHz/42.577478518/2.35051756758e5
    finitefield_shifts = MagFieldLFT.estimate_shifts_finitefield(sh, R_selected_NiSAL, B0, T, grid)
    KMcG_shifts = MagFieldLFT.calc_shifts_KurlandMcGarvey(sh, R_selected_NiSAL, T)
    return norm(finitefield_shifts-KMcG_shifts) < 1e-5
end

# Test against analytical expressions from Smith and Thornley (1966)
function test_Wybourne()
    J = 5
    Jz, Jp, Jm = MagFieldLFT.calc_lops_complex(J)
    T_Wyb = MagFieldLFT.calc_STOs_WE(J, MagFieldLFT.RME_Wybourne)
    ref_3_0 = 0.5*(5*Jz^3 - (3J*(J+1)-1)*Jz)
    ref_3_3 = -sqrt(5/16)*Jp^3
    return norm(T_Wyb[(3,0)]-ref_3_0)<1e-5 && norm(T_Wyb[(3,3)]-ref_3_3)<1e-5
end

function test_H_fieldfree_Wyb()
    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_Tb", "Tb", "Wyb_complex")
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    values = eigvals(sh.H_fieldfree)
    sort!(values)
    values = values .- values[1]
    ref = MagFieldLFT.cmm1_Hartree*[0.0, 0.18236297983901295, 72.82195951712006, 74.04535451377458,
    172.9690387781115, 191.64087572361137, 239.46631221328218, 245.49361237777234,
    263.3118254286973, 330.5619011015835, 339.36452977217255, 416.54655566755685, 418.04462983585125]
    return norm(values-ref)<1e-10
end

function test_JJbeta()
    J = 2.5
    Dtensor = zeros(3,3)
    gmatrix = 2*Matrix(1.0I, 3, 3)
    shparam = MagFieldLFT.SHParam(Int64(2J+1), gmatrix, Dtensor)
    JJ_deriv = MagFieldLFT.JJbeta(shparam)
    ref = -(J*(J+1)/3)*Matrix(1.0I, 3, 3)
    return norm(JJ_deriv-ref) <1e-10
end

# tested correctness of JJbeta2 by comparing with exact dyadic minus JJbeta at very high T
function test_JJbeta2()
    T = 298.0       # temperature in Kelvin
    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_Tb_real", "Tb")
    JJderiv2 = MagFieldLFT.JJbeta2(shparam)
    ref = [-0.008626405358262486 0.00031406199837191793 0.000584764690998459;
    0.00031406199837191793 -0.0037511527284775164 7.13592610231119e-5;
    0.000584764690998459 7.13592610231119e-5 0.012377558086740003]
    return norm(JJderiv2-ref) < 1e-10
end

# tested correctness of calc_dyadics_Wyb by comparing with sum of JJbeta and JJbeta2 terms (analytical)
# at very high T
function test_calc_dyadics_Wyb()
    T = 298.0       # temperature in Kelvin
    Ln = "Tb"
    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    exact_dyadic = MagFieldLFT.calc_dyadic(sh, T)
    ref = [-19983.04038430544 368.8536285344071 -131.44423778299375;
    368.8536285344071 -15564.84416103268 -273.97794073360933;
    -131.44423778299375 -273.97794073360933 -8543.79795610935]
    return norm(exact_dyadic - ref) < 1e-7
end

# we choose an unphysically high temperature here in order to
# suppress higher powers in beta in the exact dyadic
function test_JJbeta3()
    T = 100000      # temperature in Kelvin
    Ln = "Tb"
    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    dyadic_order2 = MagFieldLFT.calc_dyadic_order2(shparam, T)
    dyadic_order3 = MagFieldLFT.calc_dyadic_order3(shparam, T)
    exact_dyadic = MagFieldLFT.calc_dyadic(sh, T)
    residual = exact_dyadic - dyadic_order2
    beta3term = dyadic_order3 - dyadic_order2
    return norm(residual - beta3term)/norm(residual) <1e-2
end

function test_Bkq_real()
    Ln = "Tb"
    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    values = eigvals(sh.H_fieldfree)
    sort!(values)
    values = values .- values[1]
    ref = MagFieldLFT.cmm1_Hartree*[0.0, 0.18236297983901295, 72.82195951712006, 74.04535451377458,
    172.9690387781115, 191.64087572361137, 239.46631221328218, 245.49361237777234,
    263.3118254286973, 330.5619011015835, 339.36452977217255, 416.54655566755685, 418.04462983585125]
    return norm(values-ref)<1e-10
end

function test_symtensor_trafo()
    tensor = rand(3,3)
    tensor = 0.5*(tensor+tensor')   # symmetrize the random 3x3 tensor
    tensor_0_0, tensor_2 = MagFieldLFT.symtensor_trafo_Cart_sph(tensor)
    tensor_backtransformed = MagFieldLFT.symtensor_trafo_sph_Cart(tensor_0_0, tensor_2)
    return norm(tensor-tensor_backtransformed) < 1e-10
end

function test_dyadics()
    mult = 3   # NiSAL has a triplet ground state
    S =  (mult-1)/2

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*MagFieldLFT.cmm1_Hartree   # directly convert from cm-1 to Hartree

    sh = MagFieldLFT.SpinHamiltonian(mult, 2*Matrix(1.0I, 3, 3), Dtensor)

    T = 298.0

    dyadic = MagFieldLFT.calc_dyadic(sh, T)
    ref = [-713.6680177845467 -9.446126297670798 -28.19413400569129;
    -9.446126297670798 -681.0938139514958 -12.227914005431503;
    -28.19413400569129 -12.227914005431503 -720.4591260503914]
    return norm(dyadic-ref) < 1e-10
end

function test_SSbeta()
    mult = 3   # NiSAL has a triplet ground state
    S =  (mult-1)/2

    Dtensor = zeros(3,3)
    gmatrix = 2*Matrix(1.0I, 3, 3)
    shparam = MagFieldLFT.SHParam(mult, gmatrix, Dtensor)
    SS_beta = MagFieldLFT.JJbeta(shparam)
    SS_beta_ref = -S*(S+1)/3 *Matrix(1.0I, 3,3)
    return norm(SS_beta-SS_beta_ref) < 1e-15
end

function test_SSbeta2()
    mult = 3   # NiSAL has a triplet ground state
    S =  (mult-1)/2

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*MagFieldLFT.cmm1_Hartree   # directly convert from cm-1 to Hartree
    Dtensor = Dtensor - (tr(Dtensor)/3)*Matrix(1.0I, 3, 3)  # equation is only valid if D-tensor is traceless

    shparam = MagFieldLFT.SHParam(mult, 2*Matrix(1.0I, 3, 3), Dtensor)
    SS_beta2 = MagFieldLFT.JJbeta2(shparam)
    SS_beta2_ref = S*(S+1)*(2S+3)*(2S-1)/15 * Dtensor
    return norm(SS_beta2-SS_beta2_ref) < 1e-15
end

function test_SSbeta3()
    mult = 5   # For this test, use a fake multiplicity: anisotropic contribution is exactly zero for S=1
    S =  (mult-1)/2

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*MagFieldLFT.cmm1_Hartree   # directly convert from cm-1 to Hartree
    Dtensor = Dtensor - (tr(Dtensor)/3)*Matrix(1.0I, 3, 3)  # equation is only valid if D-tensor is traceless
    D2 = Dtensor * Dtensor
    trD2 = tr(D2)
    D2aniso = D2 - (trD2/3)*Matrix(1.0I, 3, 3)
    shparam = MagFieldLFT.SHParam(mult, 2*Matrix(1.0I, 3, 3), Dtensor)
    SS_beta3 = MagFieldLFT.JJbeta3(shparam)
    SS_beta3_ref =  -S*(S+1)*(2S+3)*(2S-1)*((S+2)*(S-1)/17.5*D2aniso - trD2/30*Matrix(1.0I, 3, 3))
    return norm(SS_beta3-SS_beta3_ref) < 1e-15
end

"""
For this test, I also compared with the CASSCF numbers in Table 1 of
Suturina et al., Angew. Chem. Int. Ed. 2017, 56, 12215-12218.
I got agreement for all six complexes in the table, but only test one of them here.
"""
function test_susceptibility_Ln()
    Ln = "Tb"

    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    T = 298.0
    susc = MagFieldLFT.calc_susceptibility_vanVleck(sh, T)
    susc = susc-tr(susc)/3*Matrix(1.0I, 3, 3)   # traceless part
    vals = eigvals(susc)

    au2angstrom3 = (MagFieldLFT.a0/MagFieldLFT.angstrom)^3
    vals *= au2angstrom3    # convert from atomic units to angstrom^3

    order = sortperm(vals, by=abs)   # convention: |chi_x| <= |chi_y| <= |chi_z|
    chi_x = vals[order[1]]
    chi_y = vals[order[2]]
    chi_z = vals[order[3]]

    chi_ax = chi_z - 0.5*(chi_x+chi_y)
    chi_rh = 0.5*(chi_x - chi_y)
    rhombicity = chi_rh/chi_ax

    ref_chi_ax = -0.5158916695847895
    ref_rhombicity = 0.24158864394051485

    return abs(chi_ax-ref_chi_ax)<1e-10 && abs(rhombicity-ref_rhombicity)<1e-10
end

function test_susceptibility_fromdyadic()
    Ln = "Tb"

    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    T = 298.0
    susc = MagFieldLFT.calc_susceptibility_vanVleck(sh, T)
    dyadic = MagFieldLFT.calc_dyadic(sh, T)
    susc_fromdyadic = MagFieldLFT.calc_susceptibility_fromdyadic(dyadic, shparam.gtensor)
    return norm(susc-susc_fromdyadic) < 1e-10
end

# The results differ by 3/2 for axial part and sqrt(3/2) for rhombic part (no idea why).
# When I insert parameters from Table S6 of Suturina et al. -> can reproduce experimental chi_ax, chi_rh
# When I insert exact dyadic instead of second-order one -> can reproduce CASSCF chi_ax, chi_rh
function test_Bleaney()
    T = 298.0

    # B20 and B22 constants of Tb complex from SI Table S8 of Suturina et al. (2017).
    B20 = -437*MagFieldLFT.cmm1_Hartree
    B22 = -142*MagFieldLFT.cmm1_Hartree
    #B20 = -314*MagFieldLFT.cmm1_Hartree   # Test: reproduce experimental susc
    #B22 = -250*MagFieldLFT.cmm1_Hartree   # Test: reproduce experimental susc
    mu0 = 4pi*MagFieldLFT.alpha^2
    muB = 0.5
    CJ = -157.5   # from Table 1 of Bleaney (1972)

    # Well-known equations for Bleaney's theory:
    chi_ax_Bleaney = -mu0*muB^2*CJ*B20/10/(MagFieldLFT.kB*T)^2
    chi_rh_Bleaney = -mu0*muB^2*CJ*B22/30/(MagFieldLFT.kB*T)^2

    # Susceptibility from dyadic approximated to second order in beta:
    Ln = "Tb"
    shparam = MagFieldLFT.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = MagFieldLFT.SpinHamiltonian(shparam)
    dyadic_order2 = MagFieldLFT.calc_dyadic_order2(shparam, T)
    #dyadic_order2 = MagFieldLFT.calc_dyadic(sh, T)   # For test: exact dyadic!
    susc_fromdyadic = MagFieldLFT.calc_susceptibility_fromdyadic(dyadic_order2, shparam.gtensor)
    susc_fromdyadic_aniso = susc_fromdyadic - tr(susc_fromdyadic)/3*Matrix(1.0I, 3, 3)
    vals = eigvals(susc_fromdyadic_aniso)
    order = sortperm(vals, by=abs)   # convention: |chi_x| <= |chi_y| <= |chi_z|
    chi_x = vals[order[1]]
    chi_y = vals[order[2]]
    chi_z = vals[order[3]]

    chi_ax = chi_z - 0.5*(chi_x+chi_y)
    chi_rh = 0.5*(chi_x - chi_y)

    au2angstrom3 = (MagFieldLFT.a0/MagFieldLFT.angstrom)^3
    #println(chi_ax_Bleaney*au2angstrom3)
    #println(chi_ax*au2angstrom3)
    #println(chi_rh_Bleaney/chi_ax_Bleaney)
    #println(chi_rh/chi_ax)

    return false
end

@testset "MagFieldLFT.jl" begin
    @test test_createSDs()
    @test test_createSDs2()
    @test test_U_complex2real()
    @test test_calc_lops()
    @test test_calcERIs_complex()
    @test test_calcERIs_real()
    @test test_spinorb2orbindex()
    @test test_occ_list()
    @test test_orb2spinorb_and_back()
    @test test_unocc_list()
    @test test_SD2index()
    @test test_Z_summand()
    @test test_ERIs_symmetries()
    @test test_calc_exc_equal()
    @test test_calc_exc_minus()
    @test test_calc_exc_occ2self()
    @test test_calc_exclists()
    @test test_calc_H_nonrel1()
    @test test_calc_H_nonrel2()
    @test test_calc_H_fieldfree()
    @test test_calc_L()
    @test test_calc_S()
    @test test_total_J()
    @test_broken test_read_AILFT_params_ORCA()
    @test test_Ercomplex_SOC()
    @test test_Ercomplex()
    @test test_calc_free_energy()
    @test test_average_magnetic_moment()
    @test test_average_magnetic_moment2()
    @test test_average_magnetic_moment3()
    @test test_integrate_spherical()
    @test test_integrate_spherical_Lebedev()
    @test test_dipole_matrix()
    @test test_dipole_field()
    @test test_determine_degenerate_sets()
    @test test_determine_degenerate_sets2()
    @test test_calc_susceptibility_vanVleck()
    @test test_KurlandMcGarvey()
    @test test_KurlandMcGarvey_vs_finitefield()
    @test test_KurlandMcGarvey_vs_finitefield_Lebedev()
    @test test_Fderiv2_numeric_vs_analytic()
    @test test_Fderiv3_numeric_vs_analytic()
    @test test_Fderiv4_numeric_vs_analytic()
    @test test_Fderiv4_numeric_vs_analytic_zerofield()
    @test test_KurlandMcGarvey_ord4_Br_field()
    @test test_KurlandMcGarvey_ord4_Br_temperature()
    @test test_print_composition()
    @test test_group_eigenvalues()
    @test test_print_composition2()
    @test_broken test_print_composition_ErIII()
    @test_broken test_calc_contactshift()
    @test test_KurlandMcGarvey_vs_finitefield_Lebedev_ord4()
    @test test_cubicresponse_spin()
    @test test_STOs()
    @test test_PCS_PDA_finitefield_SH()
    @test test_Wybourne()
    @test test_H_fieldfree_Wyb()
    @test test_JJbeta()
    @test test_JJbeta2()
    @test test_calc_dyadics_Wyb()
    @test test_JJbeta3()
    @test test_Bkq_real()
    @test test_symtensor_trafo()
    @test test_dyadics()
    @test test_SSbeta()
    @test test_SSbeta2()
    @test test_SSbeta3()
    @test test_susceptibility_Ln()
    @test test_susceptibility_fromdyadic()
    @test_broken test_Bleaney()
end
