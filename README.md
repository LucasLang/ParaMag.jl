`ParaMag.jl` is a Julia package for the simulation of magnetic properties of open-shell paramagnetic molecules.
Currently, two different computational models are supported: Spin Hamiltonians and ligand field theory.

A specialty of the code are higher-order derivatives of the Helmholtz free energy, which arise in the form of hypersusceptibilities and in the context of describing field-dependent chemical shieldings in NMR spectra of paramagnetic molecules.

Some functionality (e.g. reading ligand field parameters from an ORCA AILFT output file) requires the [OutputParser.jl](https://github.com/LucasLang/OutputParser.jl) package.

If you use `ParaMag.jl` in your research, please cite

[Lang et al.: Theory of Field-Dependent NMR shifts in Paramagnetic Molecules](https://doi.org/10.26434/chemrxiv-2025-1z8v9)
