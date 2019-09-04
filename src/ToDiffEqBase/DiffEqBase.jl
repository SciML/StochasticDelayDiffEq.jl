# Create SDDEProblem and solution

abstract type AbstractSDDEProblem{uType,tType,lType,isinplace,ND} <: DEProblem end
abstract type AbstractConstantLagSDDEProblem{uType,tType,lType,isinplace,ND} <:
                      AbstractSDDEProblem{uType,tType,lType,isinplace,ND} end

abstract type AbstractSDDEAlgorithm <: DEAlgorithm end

abstract type AbstractSDDEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

abstract type AbstractSDDESolution{T,N} <: AbstractRODESolution{T,N} end
