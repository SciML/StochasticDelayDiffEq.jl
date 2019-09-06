function build_solution(
        prob::AbstractSDDEProblem,
        # prob::AbstractSDDEProblem,
        alg,t,u;W=nothing,timeseries_errors = length(u) > 2,
        dense = false,dense_errors=dense,calculate_error=true,
        interp = LinearInterpolation(t,u),
        retcode = :Default,
        seed = UInt64(0), destats=nothing, kwargs...)

  T = eltype(eltype(u))
  N = length((size(prob.u0)..., length(u)))

  if typeof(prob.f) <: Tuple
    f = prob.f[1]
  else
    f = prob.f
  end

  if has_analytic(f)
    u_analytic = Vector{typeof(prob.u0)}()
    errors = Dict{Symbol,real(eltype(prob.u0))}()
    sol = RODESolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(W),
                      typeof(prob),typeof(alg),typeof(interp),typeof(destats)}(
                      u,u_analytic,errors,t,W,prob,alg,interp,dense,0,destats,retcode,seed)

    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol  
  else
    return RODESolution{T,N,typeof(u),Nothing,Nothing,typeof(t),
                        typeof(W),typeof(prob),typeof(alg),typeof(interp),typeof(destats)}(
                        u,nothing,nothing,t,W,prob,alg,interp,dense,0,destats,retcode,seed)
  end
end