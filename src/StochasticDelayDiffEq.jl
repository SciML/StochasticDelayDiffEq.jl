# __precompile__(false)
module StochasticDelayDiffEq
#########################################################
#########################################################
#########################################################
#########################################################
# using DifferentialEquations
using LinearAlgebra, StaticArrays
using Parameters, DataStructures
using Logging
using RecursiveArrayTools
using DiffEqBase
using DiffEqBase: @add_kwonly, add_kwonly
import DiffEqBase:
DEProblem, DEAlgorithm, DEIntegrator, AbstractRODEProblem,
AbstractODESolution, AbstractRODESolution, AbstractHistoryFunction, AbstractDiffEqFunction, AbstractDiffEqLinearOperator, AbstractSDEIntegrator,
isinplace, __has_analytic, has_analytic, __has_jac, __has_tgrad, __has_Wfact, __has_Wfact_t, __has_paramjac, __has_syms, __has_colorvec, RECOMPILE_BY_DEFAULT, promote_tspan, isadaptive, __solve, __init, postamble!, savevalues,
    AbstractSDEAlgorithm, AbstractRODEAlgorithm, 
    ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK
import DelayDiffEq: constant_extrapolant!, constant_extrapolant, AbstractMethodOfStepsAlgorithm, iscomposite, Discontinuity, MethodOfSteps
import OrdinaryDiffEq:  alg_maximum_order, alg_extrapolates, uses_uprev
using DiffEqNoiseProcess
using StochasticDiffEq
import StochasticDiffEq: alg_order, alg_mass_matrix_compatible, alg_compatible, alg_needs_extra_process, is_diagonal_noise, initialize_callbacks!, handle_dt!, modify_dt_for_tstops!, tstop_saveat_disc_handling, 
    EM, LambaEM, RSWM, StochasticDiffEqAlgorithm, AbstractSDEAlgorithm, StochasticDiffEqAdaptiveAlgorithm, StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm, SDEIntegrator

import RandomNumbers: Xorshifts

import Base: convert
# #########################################################
# #########################################################
# #########################################################
# #########################################################

include("ToDiffEqBase/DiffEqBase.jl")
include("ToDiffEqBase/diffeqfunction.jl")
include("ToDiffEqBase/problems/sdde_problems.jl")
include("ToDiffEqBase/problems/problem_traits.jl")
include("ToDiffEqBase/solutions/rode_solutions.jl")

include("integrators/type.jl")
include("integrators/integrator_utils.jl")
include("functionwrapper.jl")
include("history_function.jl")
include("utils.jl")
include("solve.jl")
# # include("algorithms.jl")
# # include("history_function.jl")
# # include("solve.jl")

export __init, __solve, SDDEProblem

end # module
