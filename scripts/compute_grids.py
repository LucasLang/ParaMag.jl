import os

os.system("mkdir -p lebedevgrids")

gridsizes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]

for N in gridsizes:

    Nstr = str(N)
    digits = len(Nstr)
    zeros = 4-digits
    subroutine = "LD" + zeros*"0" + Nstr


    fortranfile = fr"""program write_grids

implicit none

integer, parameter :: dp = selected_real_kind(15,307)
integer :: i

real(dp), dimension({N}) :: X
real(dp), dimension({N}) :: Y
real(dp), dimension({N}) :: Z
real(dp), dimension({N}) :: W
integer :: N

call {subroutine}(X,Y,Z,W,N)

open(1, file = 'grid_{N}', status = 'new')
do i = 1, N
    write(1,*) X(i), Y(i), Z(i), W(i)
enddo
close(1)

end program
"""

    with open("lebedevgrids/write_grid.f90", "w") as f:
        f.writelines(fortranfile)

    os.system("gfortran Lebedev-Laikov.F lebedevgrids/write_grid.f90 -o lebedevgrids/write_grid")
    os.chdir("lebedevgrids")
    os.system("./write_grid")
    os.system("rm write_grid")
    os.system("rm write_grid.f90")
    os.chdir("..")
