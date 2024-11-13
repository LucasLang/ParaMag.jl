module MagFieldLFT

using LinearAlgebra, Permutations, OutputParser, DelimitedFiles, Printf, TensorOperations, WignerSymbols

export read_AILFT_params_ORCA, LFTParam, lebedev_grids
export calc_dyadics

include("Basics.jl")
include("MagProp.jl")
include("Read.jl")
include("LFT.jl")
include("SpinHamiltonians.jl")
include("PrintComposition.jl")

# XXXLucasXXX: Something in the following function is strange -> try to understand how it works
function F_chi1_chi3_fromparam(F_calc_shift::Function)
    function calc_chi1_chi3_fromparam(model::CompModel, R::Vector{Vector{Float64}}, T::Real, B0::Real, direct::Bool=false, indirect::Bool=false)
        chi1 = calc_susceptibility_vanVleck(model, T)
        chi3 = zeros(Float64, 3, 3, 3, 3)
        if direct
            chi3 = -4pi*alpha^2*calc_F_deriv4(model, T, [0.0,0.0,0.0]) 
        end
        return F_calc_shift(chi1, chi3, R, T, B0, direct, indirect)
    end
end

calc_shifts_KurlandMcGarvey_ord4(model::CompModel, R::Vector{Vector{Float64}}, T::Real, B0::Real, direct::Bool=false, indirect::Bool=false) = F_chi1_chi3_fromparam(calc_shifts_KurlandMcGarvey_ord4)(model, R, T, B0, direct, indirect)

# XXXLucasXXX: Something in the following function is strange -> try to understand how it works
function F_chi1_fromparam(F_calc_shift::Function)
    function calc_chi1_fromparam(model::CompModel, R::Vector{Vector{Float64}}, T::Real, B0::Float64, S::Float64, gfactor::Float64, direct::Bool=false, indirect::Bool=false)
        chi1 = calc_susceptibility_vanVleck(model, T)
        return F_calc_shift(chi1, R, T, B0, S, gfactor, direct, indirect)
    end
end

calc_shifts_KurlandMcGarvey_Br(model::CompModel, R::Vector{Vector{Float64}}, T::Real, B0::Float64, S::Float64, gfactor::Float64, direct::Bool=false, indirect::Bool=false) = F_chi1_fromparam(calc_shifts_KurlandMcGarvey_Br)(model, R, T, B0, S, gfactor, direct, indirect)

end