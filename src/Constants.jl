const kB = 3.166811563e-6    # Boltzmann constant in Eh/K
const alpha = 0.0072973525693  # fine structure constant
const mu0 = 4pi*alpha^2
const cmm1_Hartree = 4.55633525e-06   # cm-1 / Ha
const a0 = 5.29177210544e-11  # Bohr radius in m (SI)
const angstrom = 1e-10        # angstrom in m (SI)
const au_time = 2.4188843265864e-17       # atomic unit of time in seconds
const au_fluxdensity = 2.35051757077e5    # atomic unit of magnetic flux density in Tesla

# values for the following two dictionaries are taken from Table 1 of
# Parker et al., Acc. Chem. Res. 2020, 53, 1520-1534.
ground_J = Dict("Tb" => 6,
"Dy" => 15/2,
"Ho" => 8,
"Er" => 15/2,
"Tm" => 6,
"Yb" => 7/2)

# The following is from Table 1 of B. Bleaney, J. Magn. Reson. 8, 91-100 (1972).
Lande_gJ = Dict("Tb" => 3/2,
"Dy" => 4/3,
"Ho" => 5/4,
"Er" => 6/5,
"Tm" => 7/6,
"Yb" => 8/7)

theta_factors = Dict(("Tb", 2) => -1/99,
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

# The following constants are taken from Table 1 in Journal of Computational Chemistry 2014, 35, 1935–1941
conv_Wyb2Stev = Dict((2,0) => 2.0,
(2,1) => 1.0/sqrt(6.0),
(2,2) => 2.0/sqrt(6.0),
(4,0) => 8.0,
(4,1) => 2.0/sqrt(5.0),
(4,2) => 4.0/sqrt(10.0),
(4,3) => 2.0/sqrt(35.0),
(4,4) => 8.0/sqrt(70.0),
(6,0) => 16.0,
(6,1) => 8.0/sqrt(42.0),
(6,2) => 16.0/sqrt(105.0),
(6,3) => 8.0/sqrt(105.0),
(6,4) => 16.0/(3.0*sqrt(14.0)),
(6,5) => 8.0/(3.0*sqrt(77.0)),
(6,6) => 16.0/sqrt(231.0)
)
