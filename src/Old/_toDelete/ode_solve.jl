function DiffEqBase.__init(prob::DiffEqBase.AbstractODEProblem,
                           alg::OrdinaryDiffEqAlgorithm,
                           recompile::Type{Val{recompile_flag}} = Val{true};
                           d_discontinuities= eltype(prob.tspan)[],
                           dense = save_everystep && !(typeof(alg) <: FunctionMap) && isempty(saveat),
                           dt = typeof(alg) <: FunctionMap && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
                           gamma = gamma_default(alg),
                           qmin = qmin_default(alg),
                           qmax = qmax_default(alg),
                           qsteady_min = qsteady_min_default(alg),
                           qsteady_max = qsteady_max_default(alg),
                           progress_name = "ODE",
                           alias_u0 = false) where recompile_flag
  if typeof(prob.f)<:DynamicalODEFunction && typeof(prob.f.mass_matrix)<:Tuple
    if any(mm != I for mm in prob.f.mass_matrix)
      error("This solver is not able to use mass matrices.")
    end
  elseif !(typeof(prob)<:DiscreteProblem) &&
         !is_mass_matrix_alg(alg) &&
         prob.f.mass_matrix != I
    error("This solver is not able to use mass matrices.")
  end



  if (((!(typeof(alg) <: OrdinaryDiffEqAdaptiveAlgorithm) && !(typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm)) || !adaptive) && dt == tType(0) && isempty(tstops)) && !(typeof(alg) <: Union{FunctionMap,LinearExponential})
      error("Fixed timestep methods require a choice of dt or choosing the tstops")
  end

  # Get the control variables

  if alias_u0
    u = prob.u0
  else
    u = recursivecopy(prob.u0)
  end

  uType = typeof(u)

  ks = Vector{uType}(undef, 0)

  uEltypeNoUnits = recursive_unitless_eltype(u)
  tTypeNoUnits   = typeof(one(tType))

  # dtmin is all abs => does not care about sign already.

  if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number && uBottomEltypeNoUnits == uBottomEltype # Could this be more efficient for other arrays?
    if !(typeof(u) <: ArrayPartition)
      rate_prototype = recursivecopy(u)
    else
      rate_prototype = similar(u, typeof.(oneunit.(recursive_bottom_eltype.(u.x))./oneunit(tType))...)
    end
  else
    if uBottomEltypeNoUnits == uBottomEltype
      rate_prototype = u
    else # has units!
      rate_prototype = u/oneunit(tType)
    end
  end
  rateType = typeof(rate_prototype) ## Can be different if united

  callbacks_internal = CallbackSet(callback,prob.callback)

  max_len_cb = DiffEqBase.max_vector_callback_length(callbacks_internal)
  if max_len_cb isa VectorContinuousCallback
    callback_cache = DiffEqBase.CallbackCache(max_len_cb.len,uBottomEltype,uBottomEltype)
  else
    callback_cache = nothing
  end

  ### Algorithm-specific defaults ###
  if save_idxs === nothing
    ksEltype = Vector{rateType}
  else
    ks_prototype = rate_prototype[save_idxs]
    ksEltype = Vector{typeof(ks_prototype)}
  end

  # Have to convert incase passed in wrong.
  if save_idxs === nothing
    timeseries = convert(Vector{uType},timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)},timeseries_init)
  end
  ts = convert(Vector{tType},ts_init)
  ks = convert(Vector{ksEltype},ks_init)
  alg_choice = Int[]

  if !adaptive && save_everystep && tspan[2]-tspan[1] != Inf
    dt == 0 ? steps = length(tstops) :
              steps = ceil(Int,internalnorm((tspan[2]-tspan[1])/dt,tspan[1]))
    sizehint!(timeseries,steps+1)
    sizehint!(ts,steps+1)
    sizehint!(ks,steps+1)
  elseif save_everystep
    sizehint!(timeseries,50)
    sizehint!(ts,50)
    sizehint!(ks,50)
  elseif !isempty(saveat_internal)
    sizehint!(timeseries,length(saveat_internal)+1)
    sizehint!(ts,length(saveat_internal)+1)
    sizehint!(ks,length(saveat_internal)+1)
  else
    sizehint!(timeseries,2)
    sizehint!(ts,2)
    sizehint!(ks,2)
  end

  if save_start
    saveiter = 1 # Starts at 1 so first save is at 2
    saveiter_dense = 1
    copyat_or_push!(ts,1,t)
    if save_idxs === nothing
      copyat_or_push!(timeseries,1,u)
      copyat_or_push!(ks,1,[rate_prototype])
    else
      copyat_or_push!(timeseries,1,u_initial,Val{false})
      copyat_or_push!(ks,1,[ks_prototype])
    end

    if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
      copyat_or_push!(alg_choice,1,1)
    end
  else
    saveiter = 0 # Starts at 0 so first save is at 1
    saveiter_dense = 0
  end


  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    id = CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,dense,cache)
    beta2 === nothing && ( beta2=beta2_default(alg.algs[cache.current]) )
    beta1 === nothing && ( beta1=beta1_default(alg.algs[cache.current],beta2) )
  else
    id = InterpolationData(f,timeseries,ts,ks,dense,cache)
    beta2 === nothing && ( beta2=beta2_default(alg) )
    beta1 === nothing && ( beta1=beta1_default(alg,beta2) )
  end

  opts = DEOptions{typeof(abstol_internal),typeof(reltol_internal),QT,tType,
                   typeof(internalnorm),typeof(internalopnorm),typeof(callbacks_internal),typeof(isoutofdomain),
                   typeof(progress_message),typeof(unstable_check),typeof(tstops_internal),
                   typeof(d_discontinuities_internal),typeof(userdata),typeof(save_idxs),
                   typeof(maxiters),typeof(tstops),typeof(saveat),
                   typeof(d_discontinuities)}(
                       maxiters,save_everystep,adaptive,abstol_internal,
                       reltol_internal,QT(gamma),QT(qmax),
                       QT(qmin),QT(qsteady_max),
                       QT(qsteady_min),QT(failfactor),tType(dtmax),
                       tType(dtmin),internalnorm,internalopnorm,save_idxs,tstops_internal,saveat_internal,
                       d_discontinuities_internal,
                       tstops,saveat,d_discontinuities,
                       userdata,progress,progress_steps,
                       progress_name,progress_message,timeseries_errors,dense_errors,
                       QT(beta1),QT(beta2),QT(qoldinit),dense,
                       save_on,save_start,save_end,callbacks_internal,isoutofdomain,
                       unstable_check,verbose,
                       calck,force_dtmin,advance_to_tstop,stop_at_next_tstop)

  destats = DiffEqBase.DEStats(0)

  if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
    sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=id,
                      alg_choice=alg_choice,
                      calculate_error = false, destats=destats)
  else
    sol = DiffEqBase.build_solution(prob,alg,ts,timeseries,
                      dense=dense,k=ks,interp=id,
                      calculate_error = false, destats=destats)
  end

  if recompile_flag == true
    FType = typeof(f)
    SolType = typeof(sol)
    cacheType = typeof(cache)
  else
    FType = Function
    SolType = DiffEqBase.AbstractODESolution
    cacheType =  OrdinaryDiffEqCache
  end

  # rate/state = (state/time)/state = 1/t units, internalnorm drops units
  eigen_est = one(uBottomEltypeNoUnits)/one(tType)
  tprev = t
  dtcache = tType(dt)
  dtpropose = tType(dt)
  iter = 0
  kshortsize = 0
  reeval_fsal = false
  u_modified = false
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  isout = false
  accept_step = false
  force_stepfail = false
  last_stepfail = false
  event_last_time = 0
  vector_event_last_time = 1
  last_event_error = zero(uBottomEltypeNoUnits)
  dtchangeable = isdtchangeable(alg)
  q11 = tTypeNoUnits(1)
  success_iter = 0
  erracc = tTypeNoUnits(1)
  dtacc = tType(1)

  integrator = ODEIntegrator{typeof(alg),isinplace(prob),uType,tType,typeof(p),typeof(eigen_est),
                             QT,typeof(tdir),typeof(k),SolType,
                             FType,cacheType,
                             typeof(opts),fsal_typeof(alg,rate_prototype),
                             typeof(last_event_error),typeof(callback_cache)}(
                             sol,u,k,t,tType(dt),f,p,uprev,uprev2,tprev,
                             alg,dtcache,dtchangeable,
                             dtpropose,tdir,eigen_est,EEst,QT(qoldinit),q11,
                             erracc,dtacc,success_iter,
                             iter,saveiter,saveiter_dense,cache,callback_cache,
                             kshortsize,force_stepfail,last_stepfail,
                             just_hit_tstop,event_last_time,vector_event_last_time,last_event_error,
                             accept_step,
                             isout,reeval_fsal,
                             u_modified,opts,destats)
  if initialize_integrator
    initialize_callbacks!(integrator, initialize_save)
    initialize!(integrator,integrator.cache)
    save_start && typeof(alg) <: CompositeAlgorithm && copyat_or_push!(alg_choice,1,integrator.cache.current)
  end

  handle_dt!(integrator)

  integrator
end

function DiffEqBase.solve!(integrator::ODEIntegrator)
  @inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir * integrator.t < top(integrator.opts.tstops)
      loopheader!(integrator)
      if check_error!(integrator) != :Success
        return integrator.sol
      end
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
      if isempty(integrator.opts.tstops)
        break
      end
    end
    handle_tstop!(integrator)
  end
  postamble!(integrator)

  f = integrator.sol.prob.f

  if DiffEqBase.has_analytic(f)
    DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors=integrator.opts.timeseries_errors,dense_errors=integrator.opts.dense_errors)
  end
  if integrator.sol.retcode != :Default
    return integrator.sol
  end
  integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol,:Success)
end


function initialize_callbacks!(integrator, initialize_save = true)
  t = integrator.t
  u = integrator.u
  callbacks = integrator.opts.callback
  integrator.u_modified = true

  u_modified = initialize!(callbacks,u,t,integrator)

  # if the user modifies u, we need to fix previous values before initializing
  # FSAL in order for the starting derivatives to be correct
  if u_modified

    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev,integrator.u)
    else
      integrator.uprev = integrator.u
    end

    if alg_extrapolates(integrator.alg)
      if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.uprev2,integrator.uprev)
      else
        integrator.uprev2 = integrator.uprev
      end
    end

    if initialize_save &&
      (any((c)->c.save_positions[2],callbacks.discrete_callbacks) ||
      any((c)->c.save_positions[2],callbacks.continuous_callbacks))
      savevalues!(integrator,true)
    end
  end

  # reset this as it is now handled so the integrators should proceed as normal
  integrator.u_modified = false
end

function tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan)
  t0, tf = tspan
  tType = eltype(tspan)
  tdir = sign(tf - t0)

  tdir_t0 = tdir * t0
  tdir_tf = tdir * tf

  # time stops
  tstops_internal = BinaryMinHeap{tType}()
  if isempty(d_discontinuities) && isempty(tstops) # TODO: Specialize more
    push!(tstops_internal, tdir_tf)
  else
    for t in tstops
      tdir_t = tdir * t
      tdir_t0 < tdir_t ≤ tdir_tf && push!(tstops_internal, tdir_t)
    end

    for t in d_discontinuities
      tdir_t = tdir * t
      tdir_t0 < tdir_t ≤ tdir_tf && push!(tstops_internal, tdir_t)
    end

    push!(tstops_internal, tdir_tf)
  end

  # saving time points
  saveat_internal = BinaryMinHeap{tType}()
  if typeof(saveat) <: Number
    if (t0:saveat:tf)[end] == tf
      for t in (t0 + saveat):saveat:tf
        push!(saveat_internal, tdir * t)
      end
    else
      for t in (t0 + saveat):saveat:(tf - saveat)
        push!(saveat_internal, tdir * t)
      end
    end
  elseif !isempty(saveat)
    for t in saveat
      tdir_t = tdir * t
      tdir_t0 < tdir_t < tdir_tf && push!(saveat_internal, tdir_t)
    end
  end

  # discontinuities
  d_discontinuities_internal = BinaryMinHeap{tType}()
  sizehint!(d_discontinuities_internal.valtree, length(d_discontinuities))
  for t in d_discontinuities
    push!(d_discontinuities_internal, tdir * t)
  end

  tstops_internal, saveat_internal, d_discontinuities_internal
end
