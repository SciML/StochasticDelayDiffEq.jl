module StochasticDelayDiffEq
#########################################################
#########################################################
#########################################################
#########################################################
# using DifferentialEquations
using Reexport
@reexport using StochasticDiffEq
# import StochasticDiffEq: calc_J, calc_J!, calc_tderivative!, calc_tderivative
using LinearAlgebra, StaticArrays
using Parameters, DataStructures
using Logging
using RecursiveArrayTools
using DiffEqBase: AbstractSDDEProblem, AbstractSDDEAlgorithm, AbstractRODESolution, AbstractRODEFunction, AbstractSDEIntegrator, AbstractSDDEIntegrator, DEIntegrator, DEAlgorithm, AbstractRODEAlgorithm, AbstractSDEAlgorithm

import DelayDiffEq: constant_extrapolant!, constant_extrapolant, AbstractMethodOfStepsAlgorithm, iscomposite, MethodOfSteps
using DiffEqNoiseProcess

import RandomNumbers: Xorshifts
using Random
import Base: convert

# #########################################################
# #########################################################
# #########################################################
# #########################################################

include("discontinuity_type.jl")
include("integrators/type.jl")
include("integrators/interface.jl")
include("integrators/utils.jl")
include("functionwrapper.jl")
include("history_function.jl")
include("utils.jl")
include("solve.jl")

end # module
