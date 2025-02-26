using Test, LinearAlgebra

using ParaMag

function test_createSDs()
    ref = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]]
    return ParaMag.create_SDs(3,2) == ref
end

function test_createSDs2()
    nel = 1
    norb = 2
    ref = [[1], [2], [3], [4]]
    return ParaMag.create_SDs(nel, norb) == ref
end

function test_U_complex2real()
    ref = [1/sqrt(2) 0 0 0 -im/sqrt(2);
           0 -1/sqrt(2) 0 im/sqrt(2) 0;
           0 0 1 0 0;
           0 1/sqrt(2) 0 im/sqrt(2) 0;
           1/sqrt(2) 0 0 0 im/sqrt(2)]
    return ParaMag.U_complex2real(2) ≈ ref
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
    lz, lplus, lminus = ParaMag.calc_lops_complex(l)
    return lz ≈ ref_lz && lplus ≈ ref_lplus && lminus ≈ ref_lminus
end

function test_calcERIs_complex()
    F = Dict(0=>3.2, 2=>(1/7/7)*1.72, 4=>(1/21/21)*2.20)
    l = 2
    ERIs = ParaMag.calcERIs_complex(l, F)
    return ERIs[2,1,4,5] ≈ (-2.12/9) && ERIs[3,2,5,5] ≈ 0.0 && ERIs[2,2,5,5] ≈ (195.92/63)
end

function test_calcERIs_real()
    F = Dict(0=>3.2, 2=>(1/7/7)*1.72, 4=>(1/21/21)*2.20)
    l = 2
    ERIs = ParaMag.calcERIs_real(l, F)
    return ERIs[2,2,5,5] ≈ (195.92/63) && ERIs[1,4,2,5] ≈ (-0.64/21)
end

function test_ERIs_symmetries()
    F = Dict(0=>3.2, 2=>1.72, 4=>2.20)
    l = 2
    ERIs = ParaMag.calcERIs_real(l, F)
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
    return ParaMag.spinorb2orbindex(2) == (1, 'β') && ParaMag.spinorb2orbindex(5) == (3, 'α')
end

function test_occ_list()
    SD = [1,3,4,7]
    occ_alpha_list = ParaMag.occ_list(SD, 'α')
    occ_beta_list = ParaMag.occ_list(SD, 'β')
    return occ_alpha_list == [(1, 1),(2, 2),(4, 4)] && occ_beta_list == [(3,2)]
end

function test_orb2spinorb_and_back()
    test1 = ParaMag.spinorb2orbindex(ParaMag.orb2spinorbindex(3,'α')) == (3,'α')
    test2 = ParaMag.orb2spinorbindex(ParaMag.spinorb2orbindex(17)...) == 17
    return test1 && test2
end

function test_unocc_list()
    SD = [1,2,4,5,6,11]
    norb = 7
    unocc_alpha_list = ParaMag.unocc_list(SD, norb, 'α')
    unocc_beta_list = ParaMag.unocc_list(SD, norb, 'β')
    return unocc_alpha_list == [2,4,5,7] && unocc_beta_list == [4,5,6,7]
end

function test_SD2index()
    N = 3
    norb = 7
    M = 2norb
    SDs = ParaMag.create_SDs(N,norb)
    test1 = ParaMag.SD2index(SDs[51], M) == 51
    test2 = ParaMag.SD2index(SDs[32], M) == 32
    return test1 && test2
end

function test_Z_summand()
    M = 14
    N = 3
    i = 2
    P = 3
    return ParaMag.Z_summand(i,P,N,M) == 11
end

function test_calc_exc_equal()
    SD = [1,2]
    norb = 3
    sigma_p = 'α'
    sigma_q = 'α'
    exc = ParaMag.calc_exc_occ2unocc(SD, norb, sigma_p, sigma_q)
    return exc == [(6,2,1,-1), (8,3,1,-1)]
end

function test_calc_exc_minus()
    SD = [3,4]
    norb = 3
    sigma_p = 'β'
    sigma_q = 'α'
    exc = ParaMag.calc_exc_occ2unocc(SD, norb, sigma_p, sigma_q)
    return exc == [(7,1,2,1), (14,3,2,-1)]
end

function test_calc_exc_occ2self()
    SD = [1,2,4]
    norb = 5
    exc_alpha_self = ParaMag.calc_exc_occ2self(SD, norb, 'α')
    exc_beta_self = ParaMag.calc_exc_occ2self(SD, norb,'β')
    test_alpha = exc_alpha_self == [(2,1,1,1)]
    test_beta = exc_beta_self == [(2,1,1,1), (2,2,2,1)]
    return test_alpha && test_beta
end

function test_calc_exclists()
    l = 1
    N = 2
    exc = ParaMag.calc_exclists(l,N)
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
    exc = ParaMag.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    F = Dict(0 => 1.0, 2 => 2.0)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = ParaMag.calc_H_nonrel(param, exc)
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
    exc = ParaMag.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = ParaMag.Racah2F(A,B,C)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = ParaMag.calc_H_nonrel(param, exc)
    E = eigvals(H)
    return E[1]≈-15 && E[22]≈3 && E[27]≈15 && E[36]≈17 && E[45]≈57
end

function test_calc_H_fieldfree()
    l=2
    N=2
    exc = ParaMag.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = ParaMag.Racah2F(A,B,C)
    zeta = 0.5
    param = LFTParam(N, norb, l, hLFT, F, zeta)
    H = ParaMag.calc_H_fieldfree(param, exc)
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
    exc = ParaMag.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = ParaMag.Racah2F(A,B,C)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = ParaMag.calc_H_nonrel(param, exc)
    C = eigvecs(H)
    Lx, Ly, Lz = ParaMag.calc_L(l, exc)
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
    exc = ParaMag.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = ParaMag.Racah2F(A,B,C)
    param = LFTParam(N, norb, l, hLFT, F, 0.0)
    H = ParaMag.calc_H_nonrel(param, exc)
    C = eigvecs(H)
    Sx, Sy, Sz = ParaMag.calc_S(l, exc)
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
    exc = ParaMag.calc_exclists(l,N)
    norb = 2l+1
    hLFT = zeros(norb,norb)
    A = 1
    B = 2
    C = 4
    F = ParaMag.Racah2F(A,B,C)
    zeta = 0.5
    param = LFTParam(N, norb, l, hLFT, F, zeta)
    H = ParaMag.calc_H_fieldfree(param, exc)
    C = eigvecs(H)
    Lx, Ly, Lz = ParaMag.calc_L(l, exc)
    Sx, Sy, Sz = ParaMag.calc_S(l, exc)
    Jx = Lx+Sx
    Jy = Ly+Sy
    Jz = Lz+Sz
    J2 = Jx*Jx + Jy*Jy + Jz*Jz
    J2val = diag(C'*J2*C)
    return J2val[1]≈6 && J2val[6]≈12 && J2val[13]≈20
end

function test_read_AILFT_params_ORCA()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = ParaMag.calc_exclists(param)
    H = ParaMag.calc_H_nonrel(param, exc)
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
    exc = ParaMag.calc_exclists(param)
    H = ParaMag.calc_H_nonrel(param, exc)
    E = eigvals(H)
    E = (E .- E[1])*27.211  # take ground state energy as reference and convert to eV
    test1 = round(real(E[29]), digits=3) == 0.020
    test2 = round(real(E[53]), digits=3) == 2.194
    test3 = round(real(E[75]), digits=3) == 2.238
    return test1 && test2 && test3
end

function test_Ercomplex_SOC()
    param = read_AILFT_params_ORCA("ErCl63-.out", "CASSCF")
    exc = ParaMag.calc_exclists(param)
    H = ParaMag.calc_H_fieldfree(param, exc)
    E = eigvals(H)
    E = (E .- E[1])*219474.63  # take ground state energy as reference and convert to cm-1
    test1 = round(real(E[17]), digits=3) ==   5862.094
    test2 = round(real(E[43]), digits=3) ==   12235.462
    test3 = round(real(E[364]), digits=3) ==  118616.544
    return test1 && test2 && test3
end

function test_calc_free_energy()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = ParaMag.calc_exclists(param)
    H_fieldfree = ParaMag.calc_H_fieldfree(param, exc)
    S = ParaMag.calc_S(param.l, exc)
    L = ParaMag.calc_L(param.l, exc)
    Mel = ParaMag.calc_magneticmoment_operator(L, S)
    B = [0.0, 0.0, 1.0e-5]
    T = 298.0
    F1 = ParaMag.calc_free_energy(H_fieldfree, Mel, B, T)
    return F1 ≈ -0.0013082816934478216
end

function test_average_magnetic_moment()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = ParaMag.calc_exclists(param)
    H_fieldfree = ParaMag.calc_H_fieldfree(param, exc)
    S = ParaMag.calc_S(param.l, exc)
    L = ParaMag.calc_L(param.l, exc)
    Mel = ParaMag.calc_magneticmoment_operator(L, S)
    B0_mol = [0.0, 0.0, 0.0]
    T = 298.0
    energies, states = ParaMag.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Hderiv = -Mel
    Fderiv1 = ParaMag.calc_F_deriv1(energies, states, Hderiv, T)
    Mel_avg = -Fderiv1
    return Mel_avg + [1.0, 1.0, 1.0] ≈ [1.0, 1.0, 1.0]    # magnetization is zero in absence of field
end

# At low temperature, magnetization should be that of the ground state (approximately MS=-3/2)
function test_average_magnetic_moment2()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    exc = ParaMag.calc_exclists(param)
    H_fieldfree = ParaMag.calc_H_fieldfree(param, exc)
    S = ParaMag.calc_S(param.l, exc)
    L = ParaMag.calc_L(param.l, exc)
    Mel = ParaMag.calc_magneticmoment_operator(L, S)
    B0_mol = [0.0, 0.0, 1.0e-4]
    T = 1.0
    energies, states = ParaMag.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Mel = ParaMag.calc_magneticmoment_operator(L,S)
    Mel_avg = ParaMag.calc_Boltzmann_averaged_ops(energies, states, Mel, T)
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
    param = ParaMag.LFTParam(nel, norb, l, hLFT, F, zeta)

    T = 0.0001
    B0_mol = [0, 0, 1.0e-7]

    H_fieldfree, Mel = ParaMag.calc_operators(param)
    energies, states = ParaMag.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Mel_avg_finitefield = ParaMag.calc_Boltzmann_averaged_ops(energies, states, Mel, T)

    return norm(Mel_avg_finitefield - [0.0, 0.0, 0.5]) < 1.0e-4
end

function test_integrate_spherical()
    Y_1_m1(theta,phi) = 0.5*sqrt(3/(2pi)) * exp(-im*phi)*sin(theta)
    Y_1_0(theta,phi) = 0.5*sqrt(3/pi) * cos(theta)
    grid = ParaMag.spherical_product_grid(50,50)
    # Results of integration: inner products [<Y_1_m1|Y_1_m1>, <Y_1_m1|Y_1_0>, <Y_1_0|Y_1_0>]
    f(x,y) = [Y_1_m1(x,y)'*Y_1_m1(x,y), Y_1_m1(x,y)'*Y_1_0(x,y), Y_1_0(x,y)'*Y_1_0(x,y)]
    integrals = ParaMag.integrate_spherical(f, grid)
    integrals = [round(x, digits=3) for x in integrals]
    return integrals ≈ [1.000, 0.000, 1.000]
end

function test_integrate_spherical_Lebedev()
    Y_1_m1(theta,phi) = 0.5*sqrt(3/(2pi)) * exp(-im*phi)*sin(theta)
    Y_1_0(theta,phi) = 0.5*sqrt(3/pi) * cos(theta)
    grid = lebedev_grids[25]
    # Results of integration: inner products [<Y_1_m1|Y_1_m1>, <Y_1_m1|Y_1_0>, <Y_1_0|Y_1_0>]
    f(x,y) = [Y_1_m1(x,y)'*Y_1_m1(x,y), Y_1_m1(x,y)'*Y_1_0(x,y), Y_1_0(x,y)'*Y_1_0(x,y)]
    integrals = ParaMag.integrate_spherical(f, grid)
    integrals = [round(x, digits=10) for x in integrals]
    return integrals ≈ [1.000, 0.000, 1.000]
end

function test_dipole_matrix()
    R = [1,5,-2.0]
    dipmat = ParaMag.calc_dipole_matrix(R)
    ref = [-0.005477225575051661 0.003042903097250923 -0.0012171612389003691;
    0.003042903097250923 0.009128709291752768 -0.006085806194501846;
    -0.0012171612389003691 -0.006085806194501846 -0.0036514837167011074]
    return dipmat ≈ ref
end

function test_dipole_field()
    R = [5,7,9]
    M = [-4.0, 3, 7.5]
    B_dip = ParaMag.dipole_field(M, R)
    ref = [2.9330997775211083e-7, 1.7331548609510163e-7, 1.2230892547235214e-7]
    return B_dip ≈ ref
end

function test_determine_degenerate_sets()
    energies = [0.0, 1e-12, 1e-11, 2.2, 2.2+1e-11, 2.2+1e-9, 7, 9, 9, 9]
    degenerate_sets = ParaMag.determine_degenerate_sets(energies)
    D = ParaMag.DegenerateSet
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
    degenerate_sets = ParaMag.determine_degenerate_sets(energies)
    D = ParaMag.DegenerateSet
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
    lft = ParaMag.LFT(param)
    T = 298
    B0_mol = [0, 0, 1.0e-7]

    # version 1
    chi = ParaMag.calc_susceptibility_vanVleck(lft, T)
    Mel_avg_linear = (1/(4pi*ParaMag.alpha^2))*chi*B0_mol

    # version 2
    H_fieldfree, Mel = ParaMag.calc_operators(param)
    energies, states = ParaMag.calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
    Mel_avg_finitefield = ParaMag.calc_Boltzmann_averaged_ops(energies, states, Mel, T)

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
    lft = ParaMag.LFT(param)
    T = 298   # I actually did not find in the paper at which temperature they recorded it!?
    shifts = ParaMag.calc_shifts_KurlandMcGarvey(lft, R_selected_NiSAL, T)
    ref = [-96.1951957001265, 30.53047905687625, 30.175351679457314, -22.804411834183288, -19.50459482031643]
    # calculated PCS at CASSCF/NEVPT2/QDPT level according to SI of paper:
    # [-61.5, 21.7, 21.3, -15.5, -13.6]
    # My shifts are larger in magnitude by around 50%, but the relative size and sign is correct
    # This could well be an artifact of only using CASSCF for determining the AILFT parameters (and/or of the LFT approximation)
    return shifts ≈ ref
end

function test_KurlandMcGarvey_vs_finitefield()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = ParaMag.LFT(param, R_selected_NiSAL)
    T = 298   # I actually did not find in the paper at which temperature they recorded it!?
    KMcG_shifts = ParaMag.calc_shifts_KurlandMcGarvey(lft, R_selected_NiSAL, T)
    grid = ParaMag.spherical_product_grid(100,100)
    B0 = 1.0e-6
    finitefield_shifts = ParaMag.estimate_shifts_finitefield(lft, B0, T, grid)
    return norm(KMcG_shifts - finitefield_shifts) < 0.1
end

function test_KurlandMcGarvey_vs_finitefield_Lebedev()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = ParaMag.LFT(param, R_selected_NiSAL)
    T = 298   # I actually did not find in the paper at which temperature they recorded it!?
    KMcG_shifts = ParaMag.calc_shifts_KurlandMcGarvey(lft, R_selected_NiSAL, T)
    grid = lebedev_grids[20]
    B0 = 1.0e-7
    finitefield_shifts = ParaMag.estimate_shifts_finitefield(lft, B0, T, grid)
    return norm(KMcG_shifts - finitefield_shifts) < 1.0e-6
end

function test_Fderiv2_numeric_vs_analytic()
    param = read_AILFT_params_ORCA("CrF63-.out", "CASSCF")
    lft = ParaMag.LFT(param)
    B0_mol = [0.0, 0.0, 1.0e-4]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv1_0 = ParaMag.calc_F_deriv1(lft, T, B0_mol)
    Fderiv1_x = ParaMag.calc_F_deriv1(lft, T, B0_mol+[h, 0, 0])
    Fderiv1_y = ParaMag.calc_F_deriv1(lft, T, B0_mol+[0, h, 0])
    Fderiv1_z = ParaMag.calc_F_deriv1(lft, T, B0_mol+[0, 0, h])
    Fderiv2 = ParaMag.calc_F_deriv2(lft, T, B0_mol)
    Fderiv2_numeric = zeros(3,3)
    Fderiv2_numeric[1, :] = (1/h)*(Fderiv1_x - Fderiv1_0)
    Fderiv2_numeric[2, :] = (1/h)*(Fderiv1_y - Fderiv1_0)
    Fderiv2_numeric[3, :] = (1/h)*(Fderiv1_z - Fderiv1_0)
    # numerically, we calculate the derivative with respect to B0, which corresponds to the perturbation operator Mel
    # but analytically, the perturbation operator (base_op of the LFT model) is +Mel
    # => the two numbers will differ by a sign
    return norm(Fderiv2+Fderiv2_numeric) < 0.1
end

function test_Fderiv3_numeric_vs_analytic()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = ParaMag.LFT(param)
    B0_mol = [0.0, 0.0, 1.0e-4]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv2_0 = ParaMag.calc_F_deriv2(lft, T, B0_mol)
    Fderiv2_x = ParaMag.calc_F_deriv2(lft, T, B0_mol+[h, 0, 0])
    Fderiv2_y = ParaMag.calc_F_deriv2(lft, T, B0_mol+[0, h, 0])
    Fderiv2_z = ParaMag.calc_F_deriv2(lft, T, B0_mol+[0, 0, h])
    Fderiv3 = ParaMag.calc_F_deriv3(lft, T, B0_mol)
    Fderiv3_numeric = zeros(3,3,3)
    Fderiv3_numeric[1, :, :] = (1/h)*(Fderiv2_x - Fderiv2_0)
    Fderiv3_numeric[2, :, :] = (1/h)*(Fderiv2_y - Fderiv2_0)
    Fderiv3_numeric[3, :, :] = (1/h)*(Fderiv2_z - Fderiv2_0)
    # numerically, we calculate the derivative with respect to B0, which corresponds to the perturbation operator Mel
    # but analytically, the perturbation operator (base_op of the LFT model) is +Mel
    # => the two numbers will differ by a sign
    # note: elements of the tensor have magnitudes that are all larger than 1e5
    return norm(Fderiv3+Fderiv3_numeric) < 8000
end

function test_Fderiv4_numeric_vs_analytic()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = ParaMag.LFT(param)
    B0_mol = [0.0, 0.0, 1.0e-4]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv3_0 = ParaMag.calc_F_deriv3(lft, T, B0_mol)
    Fderiv3_x = ParaMag.calc_F_deriv3(lft, T, B0_mol+[h, 0, 0])
    Fderiv3_y = ParaMag.calc_F_deriv3(lft, T, B0_mol+[0, h, 0])
    Fderiv3_z = ParaMag.calc_F_deriv3(lft, T, B0_mol+[0, 0, h])
    Fderiv4 = ParaMag.calc_F_deriv4(lft, T, B0_mol)
    Fderiv4_numeric = zeros(3,3,3,3)
    Fderiv4_numeric[1, :, :, :] = (1/h)*(Fderiv3_x - Fderiv3_0)
    Fderiv4_numeric[2, :, :, :] = (1/h)*(Fderiv3_y - Fderiv3_0)
    Fderiv4_numeric[3, :, :, :] = (1/h)*(Fderiv3_z - Fderiv3_0)
    # numerically, we calculate the derivative with respect to B0, which corresponds to the perturbation operator Mel
    # but analytically, the perturbation operator (base_op of the LFT model) is +Mel
    # => the two numbers will differ by a sign
    return rel_diff_norm(Fderiv4, -Fderiv4_numeric) < 1e-4
end

function test_Fderiv4_numeric_vs_analytic_zerofield()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = ParaMag.LFT(param)
    B0_mol = [0.0, 0.0, 0.0]
    h = 1.0e-10   # displacement for numerical derivative
    T = 1.0
    Fderiv3_0 = ParaMag.calc_F_deriv3(lft, T, B0_mol)
    Fderiv3_x = ParaMag.calc_F_deriv3(lft, T, B0_mol+[h, 0, 0])
    Fderiv3_y = ParaMag.calc_F_deriv3(lft, T, B0_mol+[0, h, 0])
    Fderiv3_z = ParaMag.calc_F_deriv3(lft, T, B0_mol+[0, 0, h])
    Fderiv4 = ParaMag.calc_F_deriv4(lft, T, B0_mol)
    Fderiv4_numeric = zeros(3,3,3,3)
    Fderiv4_numeric[1, :, :, :] = (1/h)*(Fderiv3_x - Fderiv3_0)
    Fderiv4_numeric[2, :, :, :] = (1/h)*(Fderiv3_y - Fderiv3_0)
    Fderiv4_numeric[3, :, :, :] = (1/h)*(Fderiv3_z - Fderiv3_0)
    # numerically, we calculate the derivative with respect to B0, which corresponds to the perturbation operator Mel
    # but analytically, the perturbation operator (base_op of the LFT model) is +Mel
    # => the two numbers will differ by a sign
    return rel_diff_norm(Fderiv4, -Fderiv4_numeric) < 1e-4
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
    ParaMag.print_composition(C, U, labels, thresh, buf)
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
    unique_values, indices = ParaMag.group_eigenvalues(values)
    ref_values = [1,2,5,7,10]
    ref_indices = [[1,2], [3], [4], [5,6,7], [8]]
    return (unique_values == ref_values) && (indices == ref_indices)
end

function test_print_composition2()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    exc = ParaMag.calc_exclists(param)
    H_fieldfree = ParaMag.calc_H_fieldfree(param, exc)
    energies_rel, states_rel = ParaMag.calc_solutions(H_fieldfree)
    C_rel_ground = states_rel[:,1]
    C_rel_first = states_rel[:,2]
    C_list, labels_list = ParaMag.adapt_basis_Hnonrel_S2_Sz(param)
    thresh = 0.99
    buf = IOBuffer()
    ParaMag.print_composition(C_rel_ground, C_list, labels_list, thresh, buf)
    printed_string_ground = String(take!(buf))
    ParaMag.print_composition(C_rel_first, C_list, labels_list, thresh, buf)
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
    exc = ParaMag.calc_exclists(param)
    H_fieldfree = ParaMag.calc_H_fieldfree(param, exc)
    energies_rel, states_rel = ParaMag.calc_solutions(H_fieldfree)
    C_rel_ground = states_rel[:,1]
    C_list, labels_list = ParaMag.adapt_basis_L2_S2_J2_Jz(param)
    C_list_alt, labels_list_alt = ParaMag.adapt_basis_L2_S2_J2_Jz(param, "term_symbols")
    thresh = 0.99
    buf = IOBuffer()
    ParaMag.print_composition(C_rel_ground, C_list, labels_list, thresh, buf)
    printed_string_ground = String(take!(buf))
    ParaMag.print_composition(C_rel_ground, C_list_alt, labels_list_alt, thresh, buf)
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
    lft = ParaMag.LFT(param)
    T = 298. 
    B0_MHz = 10. 
    B0 = B0_MHz/42.577478518/2.35051756758e5
    S = 1.0

    shift_ord4 = ParaMag.calc_shifts_KurlandMcGarvey_ord4(lft, R, T, B0, true, true)
    KMcG_shifts = ParaMag.calc_shifts_KurlandMcGarvey(lft, R, T)
    Br_shifts = ParaMag.calc_shifts_KurlandMcGarvey_Br(lft, R, T, B0, S, 2.0, true, true)

    return norm(KMcG_shifts - Br_shifts) < 1.0e-4  && norm(KMcG_shifts - shift_ord4) < 1.0e-4
end

#at very high temperature the difference in field dependent effects simulation should decrease
#due to the population of all the levels of the ground multiplet  
function test_KurlandMcGarvey_ord4_Br_temperature()
    param = read_AILFT_params_ORCA("NiSAL_HDPT.out", "CASSCF")
    lft = ParaMag.LFT(param)
    T = 400. 
    B0_MHz = 400. 
    B0 = B0_MHz/42.577478518/2.35051756758e5
    S = 1.0

    shift_ord4 = ParaMag.calc_shifts_KurlandMcGarvey_ord4(lft, R_selected_NiSAL, T, B0, true, true)
    Br_shifts = ParaMag.calc_shifts_KurlandMcGarvey_Br(lft, R_selected_NiSAL, T, B0, S, 2.0, true, true)

    return norm((shift_ord4 - Br_shifts) .* 100 ./ shift_ord4) < 5.0e-2
end


function test_calc_contactshift()
    mult = 3

    from_au = 27.2113834*8065.54477
    #from casscf-nevpt2 effective hamiltonian
    D = [3.021254 -5.482999 -16.364810; -5.482999 21.938653 -7.107412; -16.364810 -7.107412 -0.931482]
    D *= ParaMag.cmm1_Hartree  #from cm-1 to au (Hartree)

    #from casscf-nevpt2 effective hamiltonian
    g = [2.320356 0.036880 0.109631; 0.036810 2.180382 0.053365; 0.109499 0.053077 2.347002]
    #from DFT calc
    Aiso_values_MHz = [-0.2398,	0.3677	,-0.0953	,0.1169,	5.6311,	0.8206,	2.5466,	0.9669,	2.2237,	-0.2058,	0.3801,	-0.0323,	0.0943,	5.7383,	-0.1162,	-0.1048,	1.3015,	4.0604,	3.9516,	0.929,	-0.1677,	-0.2015,	-3.5469]
    Aiso_values_Hartree = Aiso_values_MHz * 1e6 * 2pi * ParaMag.au_time  # conversion from frequency to energy in atomic units: E = omega = 2pi nu
    Atensors = [Aiso*Matrix(1.0I, 3, 3) for Aiso in Aiso_values_Hartree]
    Nnuc = length(Atensors)
    gamma_1H = 2.6752e8      # proton gyromagnetic ratio in rad/s/T
    gamma_1H *= ParaMag.au_time * ParaMag.au_fluxdensity # proton gyromagnetic ratio in atomic units
    gammas = [gamma_1H for A in 1:Nnuc]

    T = 298.0

    shparam = ParaMag.SHParam(mult, g, D, Atensors, gammas)
    sh = ParaMag.SpinHamiltonian(shparam)
    shiftcon = ParaMag.calc_fieldindep_shifts(sh, T)
    
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

    lft = ParaMag.LFT(param, R)
    # field-independent shifts:
    shift_0 = ParaMag.calc_fieldindep_shifts(lft, T)
    grid = lebedev_grids[20]
    B0_MHz = 50.0  # should not be chosen too small: probably numerical precision issue
    B0 = B0_MHz/42.577478518/2.35051756758e5  # trafo from MHz to Tesla and then to atomic units
    finitefield_shifts = ParaMag.estimate_shifts_finitefield(lft, B0, T, grid)
    shift_ord4 = ParaMag.calc_shifts_2ndorder_total(lft, T, B0)
    diff_ord4 = shift_0 - shift_ord4
    diff_finitefield = shift_0 - finitefield_shifts

    return norm(diff_ord4[1]-diff_finitefield[1])<1e-7

end

function test_cubicresponse_spin()
    T = 298.0
    beta = 1/(ParaMag.kB*T)
    S = 2
    Sp = ParaMag.calc_lplusminus(S, +1)
    Sm = ParaMag.calc_lplusminus(S, -1)
    Sz = ParaMag.calc_lz(S)

    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5im * (Sp-Sm)
    Svec = [Sx, Sy, Sz]

    dim = 2S+1
    energies = zeros(dim)
    states = im*Matrix(1.0I, dim, dim)
    SiSjSkSl = ParaMag.calc_F_deriv4(energies, states, Svec, 298.0)
    ref = zeros(3,3,3,3)
    d = Matrix(1.0I, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        ref[i,j,k,l] = (beta^3/45)*S*(S+1)*(2S^2+2S+1)*(d[i,j]*d[k,l] + d[i,k]*d[j,l] + d[i,l]*d[j,k])
    end
    return norm(SiSjSkSl - ref)/norm(ref) < 1.0e-10
end

function test_STOs()
    TWE = ParaMag.calc_STOs_WE(3.0, ParaMag.RME_Lucas)
    Trec = ParaMag.calc_STOs_recursive(3.0)
    return norm(Trec[(4,-3)]-TWE[(4,-3)]) < 1e-10 && norm(Trec[(5,2)]-TWE[(5,2)]) < 1e-10
end

function test_PCS_PDA_finitefield_SH()
    gamma_1H = 2.6752e8      # proton gyromagnetic ratio in rad/s/T
    gamma_1H *= ParaMag.au_time * ParaMag.au_fluxdensity # proton gyromagnetic ratio in atomic units

    mult = 3   # NiSAL has a triplet ground state

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*ParaMag.cmm1_Hartree   # directly convert from cm-1 to Hartree

    gtensor = [2.1384111    0.0084976    0.0250646;
    0.0074791    2.0934328    0.0112682;
    0.0228213    0.0119502    2.1324169]

    Nnuc = length(R_selected_NiSAL)
    gammas = [gamma_1H for A in 1:Nnuc]
    Atensors = ParaMag.calc_Atensors_PDA(gtensor, gammas, R_selected_NiSAL)
    shparam = ParaMag.SHParam(mult, gtensor, Dtensor, Atensors, gammas)
    sh = ParaMag.SpinHamiltonian(shparam)

    T = 298.0
    grid = lebedev_grids[20]
    B0_MHz = 10.0
    B0 = B0_MHz/42.577478518/2.35051756758e5
    finitefield_shifts = ParaMag.estimate_shifts_finitefield(sh, B0, T, grid)
    KMcG_shifts = ParaMag.calc_shifts_KurlandMcGarvey(sh, R_selected_NiSAL, T)
    return norm(finitefield_shifts-KMcG_shifts) < 1e-5
end

# Test against analytical expressions from Smith and Thornley (1966)
function test_Wybourne()
    J = 5
    Jz, Jp, Jm = ParaMag.calc_lops_complex(J)
    T_Wyb = ParaMag.calc_STOs_WE(J, ParaMag.RME_Wybourne)
    ref_3_0 = 0.5*(5*Jz^3 - (3J*(J+1)-1)*Jz)
    ref_3_3 = -sqrt(5/16)*Jp^3
    return norm(T_Wyb[(3,0)]-ref_3_0)<1e-5 && norm(T_Wyb[(3,3)]-ref_3_3)<1e-5
end

function test_H_fieldfree_Wyb()
    shparam = ParaMag.SHParam_lanthanoid("Bkq_Tb", "Tb", "Wyb_complex")
    sh = ParaMag.SpinHamiltonian(shparam)
    values = eigvals(sh.H_fieldfree)
    sort!(values)
    values = values .- values[1]
    ref = ParaMag.cmm1_Hartree*[0.0, 0.18236297983901295, 72.82195951712006, 74.04535451377458,
    172.9690387781115, 191.64087572361137, 239.46631221328218, 245.49361237777234,
    263.3118254286973, 330.5619011015835, 339.36452977217255, 416.54655566755685, 418.04462983585125]
    return norm(values-ref)<1e-10
end

function test_JJbeta()
    J = 2.5
    Dtensor = zeros(3,3)
    gmatrix = 2*Matrix(1.0I, 3, 3)
    shparam = ParaMag.SHParam(Int64(2J+1), gmatrix, Dtensor)
    JJ_deriv = ParaMag.JJbeta(shparam)
    ref = -(J*(J+1)/3)*Matrix(1.0I, 3, 3)
    return norm(JJ_deriv-ref) <1e-10
end

# tested correctness of JJbeta2 by comparing with exact dyadic minus JJbeta at very high T
function test_JJbeta2()
    T = 298.0       # temperature in Kelvin
    shparam = ParaMag.SHParam_lanthanoid("Bkq_Tb_real", "Tb")
    JJderiv2 = ParaMag.JJbeta2(shparam)
    ref = [-0.008626405358262486 -0.00031406199837191793 -0.000584764690998459;
    -0.00031406199837191793 -0.0037511527284775164 7.13592610231119e-5;
    -0.000584764690998459 7.13592610231119e-5 0.012377558086740003]
    return norm(JJderiv2-ref) < 1e-10
end

# tested correctness of calc_dyadics_Wyb by comparing with sum of JJbeta and JJbeta2 terms (analytical)
# at very high T
function test_calc_dyadics_Wyb()
    T = 298.0       # temperature in Kelvin
    Ln = "Tb"
    shparam = ParaMag.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    exact_dyadic = ParaMag.calc_dyadic(sh, T)
    ref = [-19983.04038430544 -368.8536285344071 131.44423778299375;
    -368.8536285344071 -15564.84416103268 -273.97794073360933;
    131.44423778299375 -273.97794073360933 -8543.79795610935]
    return norm(exact_dyadic - ref) < 1e-7
end

# we choose an unphysically high temperature here in order to
# suppress higher powers in beta in the exact dyadic
function test_JJbeta3()
    T = 100000      # temperature in Kelvin
    Ln = "Tb"
    shparam = ParaMag.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    dyadic_order2 = ParaMag.calc_dyadic_order2(shparam, T)
    dyadic_order3 = ParaMag.calc_dyadic_order3(shparam, T)
    exact_dyadic = ParaMag.calc_dyadic(sh, T)
    residual = exact_dyadic - dyadic_order2
    beta3term = dyadic_order3 - dyadic_order2
    return norm(residual - beta3term)/norm(residual) <1e-2
end

function test_Bkq_real()
    Ln = "Tb"
    shparam = ParaMag.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    values = eigvals(sh.H_fieldfree)
    sort!(values)
    values = values .- values[1]
    ref = ParaMag.cmm1_Hartree*[0.0, 0.18236297983901295, 72.82195951712006, 74.04535451377458,
    172.9690387781115, 191.64087572361137, 239.46631221328218, 245.49361237777234,
    263.3118254286973, 330.5619011015835, 339.36452977217255, 416.54655566755685, 418.04462983585125]
    return norm(values-ref)<1e-10
end

function test_symtensor_trafo()
    tensor = rand(3,3)
    tensor = 0.5*(tensor+tensor')   # symmetrize the random 3x3 tensor
    tensor_0_0, tensor_2 = ParaMag.symtensor_trafo_Cart_sph(tensor)
    tensor_backtransformed = ParaMag.symtensor_trafo_sph_Cart(tensor_0_0, tensor_2)
    return norm(tensor-tensor_backtransformed) < 1e-10
end

function test_dyadics()
    mult = 3   # NiSAL has a triplet ground state
    S =  (mult-1)/2

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*ParaMag.cmm1_Hartree   # directly convert from cm-1 to Hartree

    sh = ParaMag.SpinHamiltonian(mult, 2*Matrix(1.0I, 3, 3), Dtensor)

    T = 298.0

    dyadic = ParaMag.calc_dyadic(sh, T)
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
    shparam = ParaMag.SHParam(mult, gmatrix, Dtensor)
    SS_beta = ParaMag.JJbeta(shparam)
    SS_beta_ref = -S*(S+1)/3 *Matrix(1.0I, 3,3)
    return norm(SS_beta-SS_beta_ref) < 1e-15
end

function test_SSbeta2()
    mult = 3   # NiSAL has a triplet ground state
    S =  (mult-1)/2

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*ParaMag.cmm1_Hartree   # directly convert from cm-1 to Hartree
    Dtensor = Dtensor - (tr(Dtensor)/3)*Matrix(1.0I, 3, 3)  # equation is only valid if D-tensor is traceless

    shparam = ParaMag.SHParam(mult, 2*Matrix(1.0I, 3, 3), Dtensor)
    SS_beta2 = ParaMag.JJbeta2(shparam)
    SS_beta2_ref = S*(S+1)*(2S+3)*(2S-1)/15 * Dtensor
    return norm(SS_beta2-SS_beta2_ref) < 1e-15
end

function test_SSbeta3()
    mult = 5   # For this test, use a fake multiplicity: anisotropic contribution is exactly zero for S=1
    S =  (mult-1)/2

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*ParaMag.cmm1_Hartree   # directly convert from cm-1 to Hartree
    Dtensor = Dtensor - (tr(Dtensor)/3)*Matrix(1.0I, 3, 3)  # equation is only valid if D-tensor is traceless
    D2 = Dtensor * Dtensor
    trD2 = tr(D2)
    D2aniso = D2 - (trD2/3)*Matrix(1.0I, 3, 3)
    shparam = ParaMag.SHParam(mult, 2*Matrix(1.0I, 3, 3), Dtensor)
    SS_beta3 = ParaMag.JJbeta3(shparam)
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

    shparam = ParaMag.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    T = 298.0
    susc = ParaMag.calc_susceptibility_vanVleck(sh, T)
    susc = susc-tr(susc)/3*Matrix(1.0I, 3, 3)   # traceless part
    vals = eigvals(susc)

    au2angstrom3 = (ParaMag.a0/ParaMag.angstrom)^3
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

    shparam = ParaMag.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    T = 298.0
    susc = ParaMag.calc_susceptibility_vanVleck(sh, T)
    dyadic = ParaMag.calc_dyadic(sh, T)
    susc_fromdyadic = ParaMag.calc_susceptibility_fromdyadic(dyadic, shparam.gtensor)
    return norm(susc-susc_fromdyadic) < 1e-10
end

# The results differ by 3/2 for axial part and sqrt(3/2) for rhombic part (no idea why).
# When I insert parameters from Table S6 of Suturina et al. -> can reproduce experimental chi_ax, chi_rh
# When I insert exact dyadic instead of second-order one -> can reproduce CASSCF chi_ax, chi_rh
function test_Bleaney()
    T = 298.0

    # B20 and B22 constants of Tb complex from SI Table S8 of Suturina et al. (2017).
    B20 = -437*ParaMag.cmm1_Hartree
    B22 = -142*ParaMag.cmm1_Hartree
    #B20 = -314*ParaMag.cmm1_Hartree   # Test: reproduce experimental susc
    #B22 = -250*ParaMag.cmm1_Hartree   # Test: reproduce experimental susc
    mu0 = 4pi*ParaMag.alpha^2
    muB = 0.5
    CJ = -157.5   # from Table 1 of Bleaney (1972)

    # Well-known equations for Bleaney's theory:
    chi_ax_Bleaney = -mu0*muB^2*CJ*B20/10/(ParaMag.kB*T)^2
    chi_rh_Bleaney = -mu0*muB^2*CJ*B22/30/(ParaMag.kB*T)^2

    # Susceptibility from dyadic approximated to second order in beta:
    Ln = "Tb"
    shparam = ParaMag.SHParam_lanthanoid("Bkq_$(Ln)_real", Ln)
    sh = ParaMag.SpinHamiltonian(shparam)
    dyadic_order2 = ParaMag.calc_dyadic_order2(shparam, T)
    #dyadic_order2 = ParaMag.calc_dyadic(sh, T)   # For test: exact dyadic!
    susc_fromdyadic = ParaMag.calc_susceptibility_fromdyadic(dyadic_order2, shparam.gtensor)
    susc_fromdyadic_aniso = susc_fromdyadic - tr(susc_fromdyadic)/3*Matrix(1.0I, 3, 3)
    vals = eigvals(susc_fromdyadic_aniso)
    order = sortperm(vals, by=abs)   # convention: |chi_x| <= |chi_y| <= |chi_z|
    chi_x = vals[order[1]]
    chi_y = vals[order[2]]
    chi_z = vals[order[3]]

    chi_ax = chi_z - 0.5*(chi_x+chi_y)
    chi_rh = 0.5*(chi_x - chi_y)

    au2angstrom3 = (ParaMag.a0/ParaMag.angstrom)^3
    #println(chi_ax_Bleaney*au2angstrom3)
    #println(chi_ax*au2angstrom3)
    #println(chi_rh_Bleaney/chi_ax_Bleaney)
    #println(chi_rh/chi_ax)

    return false
end

function test_PDA_SH_general_vs_specific()
    gamma_1H = 2.6752e8      # proton gyromagnetic ratio in rad/s/T
    gamma_1H *= ParaMag.au_time * ParaMag.au_fluxdensity # proton gyromagnetic ratio in atomic units

    mult = 3   # NiSAL has a triplet ground state

    Dtensor = [  3.053926    -5.555174   -16.580693;
    -5.555174    22.210495    -7.191116;
   -16.580693    -7.191116    -0.939858]*ParaMag.cmm1_Hartree   # directly convert from cm-1 to Hartree

    gtensor = [2.1384111    0.0084976    0.0250646;
    0.0074791    2.0934328    0.0112682;
    0.0228213    0.0119502    2.1324169]

    Nnuc = length(R_selected_NiSAL)
    gammas = [gamma_1H for A in 1:Nnuc]
    Atensors = ParaMag.calc_Atensors_PDA(gtensor, gammas, R_selected_NiSAL)
    shparam = ParaMag.SHParam(mult, gtensor, Dtensor, Atensors, gammas)
    sh = ParaMag.SpinHamiltonian(shparam)

    T = 298.0
    KMcG_shifts = ParaMag.calc_shifts_KurlandMcGarvey(sh, R_selected_NiSAL, T)
    shifts_generalfunction = ParaMag.calc_fieldindep_shifts(sh, T)
    return norm(KMcG_shifts-shifts_generalfunction) < 1e-10
end

function test_fielddep_contactshifts()
    gamma_1H = 2.6752e8      # proton gyromagnetic ratio in rad/s/T
    gamma_1H *= ParaMag.au_time * ParaMag.au_fluxdensity # proton gyromagnetic ratio in atomic units

    mult = 3   # NiSAL has a triplet ground state
    S = (mult-1)/2

    Dtensor = [3.021254 -5.482999 -16.364810; -5.482999 21.938653 -7.107412; -16.364810 -7.107412 -0.931482]
    Dtensor *= ParaMag.cmm1_Hartree  #from cm-1 to au (Hartree)

    #from casscf-nevpt2 effective hamiltonian
    gtensor = [2.320356 0.036880 0.109631; 0.036810 2.180382 0.053365; 0.109499 0.053077 2.347002]

    Nnuc = length(R_selected_NiSAL)
    gammas = [gamma_1H for A in 1:Nnuc]
    Aiso_values_MHz = [-0.2398,	0.3677	,-0.0953	,0.1169,	5.6311,	0.8206,	2.5466,	0.9669,	2.2237,	-0.2058,	0.3801,	-0.0323,	0.0943,	5.7383,	-0.1162,	-0.1048,	1.3015,	4.0604,	3.9516,	0.929,	-0.1677,	-0.2015,	-3.5469]
    indices = [1, 23, 5, 14]
    Aiso_values_MHz = Aiso_values_MHz[indices]    # select only some of them
    Aiso_values_Hartree = Aiso_values_MHz * 1e6 * 2pi * ParaMag.au_time  # conversion from frequency to energy in atomic units: E = omega = 2pi nu
    Atensors = [Aiso*Matrix(1.0I, 3, 3) for Aiso in Aiso_values_Hartree]
    shparam = ParaMag.SHParam(mult, gtensor, Dtensor, Atensors, gammas)
    sh = ParaMag.SpinHamiltonian(shparam)

    T = 298.0
    B0_MHz = 1200.0
    B0 = B0_MHz/42.577478518/2.35051756758e5

    shifts_indirect = ParaMag.calc_fielddep_shifts_2ndorder_indirect(sh, T, B0)
    shifts_direct = ParaMag.calc_fielddep_shifts_2ndorder_direct(sh, T, B0)
    ref_indirect = [-0.0008741457857955934, -0.012929556662378742, 0.020527115656353536, 0.0209178930885359]
    ref_direct = [0.06892304308923315, 1.0194459613561342, -1.6184843533769009, -1.6492956553750897]

    return norm(shifts_indirect-ref_indirect)<1e-5 && norm(shifts_direct-ref_direct)<1e-5
end

@testset "ParaMag.jl" begin
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
    @test test_calc_contactshift()
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
    @test test_PDA_SH_general_vs_specific()
    @test test_fielddep_contactshifts()
end
