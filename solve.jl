function DiffEqBase.__solve(prob::AbstractSDDEProblem, # TODO: DiffEqBase.AbstractSDDEProblem
                            alg::AbstractMethodOfStepsAlgorithm, args...;
                            kwargs...)
  integrator = DiffEqBase.__init(prob, alg, args...; kwargs...)
  DiffEqBase.solve!(integrator)
  integrator.sol
end
