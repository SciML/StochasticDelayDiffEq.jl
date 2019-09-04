"""
    has_constant_lags(integrator::DDEIntegrator)

Return if the DDE problem of the `integrator` contains constant delays.
"""
has_constant_lags(integrator::SDDEIntegrator) = has_constant_lags(integrator.sol.prob)

"""
    has_dependent_lags(integrator::DDEIntegrator)

Return if the DDE problem of the `integrator` contains dependent delays.
"""
has_dependent_lags(integrator::SDDEIntegrator) = has_dependent_lags(integrator.sol.prob)

"""
    has_constant_lags(prob::DDEProblem)

Return if the DDE problem `prob` contains constant delays.
"""
has_constant_lags(prob::SDDEProblem) =
  prob.constant_lags !== nothing && !isempty(prob.constant_lags)

"""
    has_dependent_lags(prob::DDEProblem)

Return if the DDE problem `prob` contains dependent delays.
"""
has_dependent_lags(prob::SDDEProblem) =
  prob.dependent_lags !== nothing && !isempty(prob.dependent_lags)


"""
    solution_arrays(u, tspan, rate_prototype; kwargs...)

Return arrays of saved time points, states, and rates, initialized with the solution at the
first time point if `save_start = true` (the default).
"""
function solution_arrays(u::uType, tspan, rate_prototype;
  timeseries_init = typeof(u)[],
  ts_init = eltype(tspan)[],
  save_idxs = nothing,
  save_start = true) where {uType}
# determine types of time and state
  # uType = typeof(u)
  tType = eltype(tspan)

  # initialize vector of saved time points
  ts = convert(Vector{tType}, ts_init)
  alg_choice = Int[]

  # initialize vector of saved states
  if save_idxs === nothing
    timeseries = convert(Vector{uType}, timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)}, timeseries_init)
  end

  # save solution at initial time point
  if save_start
    saveiter=1
    copyat_or_push!(ts, 1, first(tspan))
    if save_idxs === nothing
      copyat_or_push!(timeseries, 1, u)
    else
      copyat_or_push!(timeseries, 1, u_initial, Val{false})
    end
  else
    saveiter=0;
  end

  ts, timeseries, saveiter
end
