getalg(alg0::DEAlgorithm) = alg0
getalg(alg0::MethodOfSteps) = alg0.alg

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

function u_uprev(u0; alias_u0 = false)
  if  typeof(u0) <: Tuple
    u = ArrayPartition(prob.u0,Val{true})
  else
    if alias_u0
      u = u0
    else
      u = recursivecopy(u0)
    end
  end
  uprev = recursivecopy(u)
  u, uprev
end

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
"""
    sizehint!(sol::DESolution, n)

Suggest that solution `sol` reserves capacity for at least `n` elements.
"""
function Base.sizehint!(sol::RODESolution, n)
  sizehint!(sol.u, n)
  sizehint!(sol.t, n)

  nothing
end

"""
    sizehint!(sol::DESolution, alg, tspan, tstops, saveat; kwargs...)

Suggest that solution `sol` reserves capacity for a number of elements that
depends on the parameter settings of the numerical solver.
"""
function Base.sizehint!(sol::RODESolution, alg, tspan, tstops, saveat;
                        save_everystep = isempty(saveat),
                        adaptive = StochasticDiffEq.isadaptive(getalg(alg)), 
                        internalnorm = DiffEqBase.ODE_DEFAULT_NORM, 
                        dt = zero(eltype(tspan)))
  # obtain integration time
  t0 = first(tspan)
  integrationtime = last(tspan) - t0

  if !adaptive && save_everystep && tspan[2]-tspan[1] != Inf
    # determine number of steps if known a priori
    if iszero(dt)
      steps = length(tstops)
    else
      steps = ceil(Int, internalnorm(integrationtime / dt, t0))
    end

    sizehint!(sol,steps+1)
  elseif save_everystep
    sizehint!(sol,50)
  elseif !isempty(saveat)
    sizehint!(sol,length(saveat)+1)
  else
    sizehint!(sol,2)
  end

  nothing
end

function build_history_function(prob, alg, reltol, rate_prototype, noise_rate_prototype, W, _seed, dense;
  dt = zero(eltype(prob.tspan)),
  adaptive = StochasticDiffEq.isadaptive(getalg(alg)),
  calck = false,
  internalnorm = DiffEqBase.ODE_DEFAULT_NORM)
  @unpack f, g, h, u0, tspan, p = prob

  t0 = first(tspan)
  tType = eltype(tspan)
  tTypeNoUnits = typeof(one(tType))
  tdir = sign(last(tspan) - t0)

  uEltypeNoUnits = recursive_unitless_eltype(u0)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u0)

  # bootstrap an SDE integrator
  # - whose solution captures the dense history of the simulation
  # - that is used for extrapolation of the history for time points past the
  #   already fixed history
  # - that is used for interpolation of the history for time points in the
  #   current integration step (so the interpolation is fixed while updating the stages)
  # we wrap the user-provided history function such that function calls during the setup
  # of the integrator do not fail
  # sde_f = SDEFunctionWrapper(f, prob.h)
  # sde_g = SDEDiffusionTermWrapper(g,prob.h)

  sde_f,sde_g = wrap_functions_and_history(f, g, h)
  sde_prob = SDEProblem{isinplace(prob)}(sde_f, sde_g, u0, tspan, p)

  # get states of ODE integrator (do not alias uprev)
  sde_u, sde_uprev = u_uprev(u0; alias_u0 = false)

  # # initialize output arrays
  sde_ts, sde_timeseries, sde_saveiter = solution_arrays(sde_u, tspan, rate_prototype,
                                      save_idxs = nothing,
                                      save_start = true)

  # # obtain cache (we alias uprev2 and uprev)
  sde_cache = StochasticDiffEq.alg_cache(getalg(alg),prob, sde_u, W.dW, W.dZ, p, rate_prototype,noise_rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, sde_uprev, f, t0, dt, Val{isinplace(prob)})

  # build dense interpolation of history
  id = StochasticDiffEq.LinearInterpolationData(sde_timeseries,sde_ts)
  if typeof(getalg(alg)) <: StochasticDiffEq.StochasticDiffEqCompositeAlgorithm
    alg_choice = Int[]
    sde_sol = DiffEqBase.build_solution(prob,alg,sde_ts,sde_timeseries,W=W,
                                    destats = DiffEqBase.DEStats(0),
                                    calculate_error = false, alg_choice=alg_choice,
                                    interp = id, dense = dense, seed = _seed)
  else
    sde_sol = DiffEqBase.build_solution(prob,alg,sde_ts,sde_timeseries,W=W,
                                    destats = DiffEqBase.DEStats(0),
                                    calculate_error = false,
                                    interp = id, dense = dense, seed = _seed)
  end

  # # reserve capacity
  sizehint!(sde_sol, getalg(alg), tspan, (), ();
    save_everystep = true, adaptive = adaptive, internalnorm = internalnorm, dt = tType(dt))

  # # create simple integrator
  sde_integrator = HistorySDEIntegrator{typeof(getalg(alg)),isinplace(prob),typeof(prob.u0),
            tType, typeof(sde_sol),typeof(sde_cache)}(
              sde_sol, sde_u, t0, zero(tType), sde_uprev,
              t0, getalg(alg), zero(tType), tdir, 1, sde_cache)

  # # combine the user-provided history function and the ODE integrator with dense solution
  # # to a joint dense history of the DDE
  # # we use this history information to create a problem function of the DDE with all
  # # available history information that is of the form f(du,u,p,t) or f(u,p,t) such that
  # # SDE algorithms can be applied
  HistoryFunction(prob.h, sde_integrator)
end