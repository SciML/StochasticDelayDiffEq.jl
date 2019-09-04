function solution_arrays(u, tspan, rate_prototype;
                         timeseries_init = typeof(u)[],
                         ts_init = eltype(tspan)[],
                         ks_init = [],
                         save_idxs = nothing,
                         save_start = true)
  # determine types of time and state
  uType = typeof(u)
  tType = eltype(tspan)

  # initialize vector of saved time points
  ts = convert(Vector{tType}, ts_init)

  # initialize vector of saved states
  if save_idxs === nothing
    timeseries = convert(Vector{uType}, timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)}, timeseries_init)
  end

  # initialize vector of saved rates
  if save_idxs === nothing
    ksEltype = Vector{typeof(rate_prototype)}
  else
    ks_prototype = rate_prototype[save_idxs]
    ksEltype = Vector{typeof(ks_prototype)}
  end
  ks = convert(Vector{ksEltype}, ks_init)

  # save solution at initial time point
  if save_start
    copyat_or_push!(ts, 1, first(tspan))
    if save_idxs === nothing
      copyat_or_push!(timeseries, 1, u)
      copyat_or_push!(ks, 1, [rate_prototype])
    else
      u_initial = u[save_idxs]
      copyat_or_push!(timeseries, 1, u_initial, Val{false})
      copyat_or_push!(ks, 1, [ks_prototype])
    end
  end

  ts, timeseries, ks
end

function build_history_function(prob, alg, rate_prototype, reltol;
                                dt = zero(eltype(prob.tspan)),
                                adaptive = DiffEqBase.isadaptive(alg.alg),
                                calck = false,
                                internalnorm = DiffEqBase.ODE_DEFAULT_NORM)
  @unpack f ,g, u0, tspan, p = prob

  t0 = first(tspan)
  tType = eltype(tspan)
  tTypeNoUnits = typeof(one(tType))
  tdir = sign(last(tspan) - t0)

  uEltypeNoUnits = recursive_unitless_eltype(u0)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u0)

  # bootstrap an ODE integrator
  # - whose solution captures the dense history of the simulation
  # - that is used for extrapolation of the history for time points past the
  #   already fixed history
  # - that is used for interpolation of the history for time points in the
  #   current integration step (so the interpolation is fixed while updating the stages)
  # we wrap the user-provided history function such that function calls during the setup
  # of the integrator do not fail
  sde_f = SDEFunctionWrapper(f, prob.h)
  sde_g = SDEDiffusionFunctionWrapper(g,prob.h)
  sde_prob = SDEProblem{isinplace(prob)}(sde_f,sde_g, u0, tspan, p)

  # get states of ODE integrator (do not alias uprev)
  sde_u = recursivecopy(u0)
  sde_uprev = recursivecopy(u)

  # initialize output arrays
  sde_k = typeof(rate_prototype)[]
  sde_ts, sde_timeseries, sde_ks = DelayDiffEq.solution_arrays(sde_u, tspan, rate_prototype, save_idxs = nothing, save_start = true)

  #####################################
  ########## TODO from here ##########
  #####################################
  # obtain cache (we alias uprev2 and uprev)
  ode_cache = StochasticDiffEq.alg_cache(alg.alg, sde_u, rate_prototype, uEltypeNoUnits,
                                       uBottomEltypeNoUnits, tTypeNoUnits, sde_uprev,
                                       sde_uprev, sde_f, sde_g, t0, zero(tType), reltol, p, calck,
                                       Val{isinplace(prob)})
  (alg::IIF2M,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  # build dense interpolation of history
  if iscomposite(alg)
    ode_alg_choice = Int[]
    ode_id = OrdinaryDiffEq.CompositeInterpolationData(ode_f, ode_timeseries, ode_ts, ode_ks,
                                                       ode_alg_choice, true, ode_cache) # dense = true
    ode_sol = DiffEqBase.build_solution(ode_prob, alg.alg, ode_ts, ode_timeseries;
                                        dense = true, k = ode_ks, interp = ode_id,
                                        alg_choice = ode_alg_choice,
                                        calculate_error = false, destats = DiffEqBase.DEStats(0))
  else
    ode_id = OrdinaryDiffEq.InterpolationData(ode_f, ode_timeseries, ode_ts, ode_ks, true, ode_cache) # dense = true
    ode_sol = DiffEqBase.build_solution(ode_prob, alg.alg, ode_ts, ode_timeseries;
                                        dense = true, k = ode_ks, interp = ode_id,
                                        calculate_error = false, destats = DiffEqBase.DEStats(0))
  end

  # reserve capacity
  sizehint!(ode_sol, alg.alg, tspan, (), ();
            save_everystep = true, adaptive = adaptive, internalnorm = internalnorm, dt = tType(dt))

  # create simple integrator
  tdirType = typeof(sign(zero(tType)))
  ode_integrator = HistoryODEIntegrator{typeof(alg.alg),isinplace(prob),typeof(prob.u0),
                                        tType,tdirType,typeof(ode_k),
                                        typeof(ode_sol),typeof(ode_cache)}(
                                          ode_sol, ode_u, ode_k, t0, zero(tType), ode_uprev,
                                          t0, alg.alg, zero(tType), tdir, 1, 1, ode_cache)

  # combine the user-provided history function and the ODE integrator with dense solution
  # to a joint dense history of the DDE
  # we use this history information to create a problem function of the DDE with all
  # available history information that is of the form f(du,u,p,t) or f(u,p,t) such that
  # ODE algorithms can be applied
  HistoryFunction(prob.h, ode_integrator)
end
