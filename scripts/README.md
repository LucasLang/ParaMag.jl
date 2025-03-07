The Lebedev grids hard-coded in `ParaMag.jl` for performing numerical integration over the unit sphere for finite field shifts were generated using the script `compute_grids.py`.

In order to run the script, you need to download the Fortran source file `Lebedev-Laikov.F` for generating Lebedev quadrature grids from

http://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/

In order to generate all grid files (columns are X, Y, Z, weight), run:

```
python compute_grids.py
```

For each grid, the script generates a Fortran source file, compiles it (if you want to use a different compiler than `gfortran`, you need to manually change that in `compute_grids.py`), and runs it to generate the grid.
