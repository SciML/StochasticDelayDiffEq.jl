module StochasticDelayDiffEq

#########################################################
#########################################################
#########################################################
#########################################################
using DifferentialEquations
using LinearAlgebra
using DifferentialEquations.DiffEqBase: @add_kwonly, add_kwonly
import DifferentialEquations.DiffEqBase:
DEProblem, DEAlgorithm, DEIntegrator,
AbstractODESolution, AbstractRODESolution, AbstractHistoryFunction, AbstractDiffEqFunction, AbstractDiffEqLinearOperator,
 __has_analytic, __has_jac, __has_tgrad, __has_Wfact, __has_Wfact_t, __has_paramjac, __has_syms, __has_colorvec, RECOMPILE_BY_DEFAULT, promote_tspan
import Base: convert
#########################################################
#########################################################
#########################################################
#########################################################

include("ToDiffEqBase/DiffEqBase.jl")
include("ToDiffEqBase/diffeqfunction.jl")
include("ToDiffEqBase/problems/sdde_problems.jl")

export SDDEProblem

end # module
