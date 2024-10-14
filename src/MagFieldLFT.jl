module MagFieldLFT

using LinearAlgebra, Permutations, OutputParser, DelimitedFiles, Printf, TensorOperations

export read_AILFT_params_ORCA, LFTParam, lebedev_grids

include("Basics.jl")
include("Read.jl")
include("LFT.jl")
include("SpinHamiltonians.jl")
include("MagProp.jl")
include("PrintComposition.jl")

function F_deriv_param2states(calc_F_derivx::Function)
    function calc_F_deriv_param(param::LFTParam, T::Real, B0_mol::Vector{Float64})
        H_fieldfree, Mel = calc_operators_SDbasis(param)
        Hderiv = -Mel
        energies, states = calc_solutions_magfield(H_fieldfree, Mel, B0_mol)
        return calc_F_derivx(energies, states, Hderiv, T)
    end
    return calc_F_deriv_param
end

calc_F_deriv1(param::LFTParam, T::Real, B0_mol::Vector{Float64}) = F_deriv_param2states(calc_F_deriv1)(param, T, B0_mol)
calc_F_deriv2(param::LFTParam, T::Real, B0_mol::Vector{Float64}) = F_deriv_param2states(calc_F_deriv2)(param, T, B0_mol)
calc_F_deriv3(param::LFTParam, T::Real, B0_mol::Vector{Float64}) = F_deriv_param2states(calc_F_deriv3)(param, T, B0_mol)
calc_F_deriv4(param::LFTParam, T::Real, B0_mol::Vector{Float64}) = F_deriv_param2states(calc_F_deriv4)(param, T, B0_mol)

function F_chi1_chi3_fromparam(F_calc_shift::Function)
    function calc_chi1_chi3_fromparam(param::LFTParam, R::Vector{Vector{Float64}}, T::Real, B0::Real, direct::Bool=false, indirect::Bool=false)
        chi1 = calc_susceptibility_vanVleck(param, T)
        chi3 = zeros(Float64, 3, 3, 3, 3)
        if direct
            chi3 = -4pi*alpha^2*calc_F_deriv4(param, T, [0.0,0.0,0.0]) 
        end
        return F_calc_shift(chi1, chi3, R, T, B0, direct, indirect)
    end
end

calc_shifts_KurlandMcGarvey_ord4(param::LFTParam, R::Vector{Vector{Float64}}, T::Real, B0::Real, direct::Bool=false, indirect::Bool=false) = F_chi1_chi3_fromparam(calc_shifts_KurlandMcGarvey_ord4)(param, R, T, B0, direct, indirect)

function F_chi1_fromparam(F_calc_shift::Function)
    function calc_chi1_fromparam(param::LFTParam, R::Vector{Vector{Float64}}, T::Real, B0::Float64, S::Float64, gfactor::Float64, direct::Bool=false, indirect::Bool=false)
        chi1 = calc_susceptibility_vanVleck(param, T)
        return F_calc_shift(chi1, R, T, B0, S, gfactor, direct, indirect)
    end
end

calc_shifts_KurlandMcGarvey_Br(param::LFTParam, R::Vector{Vector{Float64}}, T::Real, B0::Float64, S::Float64, gfactor::Float64, direct::Bool=false, indirect::Bool=false) = F_chi1_fromparam(calc_shifts_KurlandMcGarvey_Br)(param, R, T, B0, S, gfactor, direct, indirect)

end