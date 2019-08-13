# Create SDDEProblem and solution

abstract type AbstractSDDEProblem{uType,tType,lType,isinplace,ND} <: DEProblem end
abstract type AbstractConstantLagSDDEProblem{uType,tType,lType,isinplace,ND} <:
                      AbstractSDDEProblem{uType,tType,lType,isinplace,ND} end

abstract type AbstractSDDEAlgorithm <: DEAlgorithm end

abstract type AbstractSDDEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

abstract type AbstractSDDESolution{T,N} <: AbstractRODESolution{T,N} end


abstract type AbstractSDDEFunction{iip} <: AbstractDiffEqFunction{iip} end

struct SDDEFunction{iip,F,G,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,GG,TCV} <: AbstractSDDEFunction{iip}
  f::F
  g::G
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  ggprime::GG
  syms::S
  colorvec::TCV
end
