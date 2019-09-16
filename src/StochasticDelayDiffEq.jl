module StochasticDelayDiffEq
#########################################################
#########################################################
#########################################################
#########################################################
# using DifferentialEquations
using Reexport
@reexport using StochasticDiffEq
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
# # include("algorithms.jl")

end # module
