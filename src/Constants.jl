const kB = 3.166811563e-6    # Boltzmann constant in Eh/K
const alpha = 0.0072973525693  # fine structure constant
const mu0 = 4pi*alpha^2
const cmm1_Hartree = 4.55633525e-06   # cm-1 / Ha

# values for the following two dictionaries are taken from Table 1 of
# Parker et al., Acc. Chem. Res. 2020, 53, 1520-1534.
ground_J = Dict("Eu" => 0,
"Tb" => 6,
"Dy" => 15/2,
"Ho" => 8,
"Er" => 15/2,
"Tm" => 6,
"Yb" => 7/2)

theta_factors = Dict(("Eu", 2) => 0,
("Eu", 4) => 0,
("Eu", 6) => 0,
("Tb", 2) => -1/99,
("Tb", 4) => 2/11/1485,
("Tb", 6) => -1/13/33/2079,
("Dy", 2) => -2/9/35,
("Dy", 4) => -8/11/45/273,
("Dy", 6) => 4/11^2/13^2/3^3/7,
("Ho", 2) => -1/30/15,
("Ho", 4) => -1/11/10/273,
("Ho", 6) => -5/11^2/13^2/3^3/7,
("Er", 2) => 4/45/35,
("Er", 4) => 2/11/15/273,
("Er", 6) => 8/11^2/13^2/3^3/7,
("Tm", 2) => 1/99,
("Tm", 4) => 8/3/11/1485,
("Tm", 6) => -5/13/33/2079,
("Yb", 2) => 2/63,
("Yb", 4) => -2/77/15,
("Yb", 6) => 4/13/33/63,
)