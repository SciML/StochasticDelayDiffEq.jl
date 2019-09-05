@inline function loopheader!(integrator::SDDEIntegrator)
    # Apply right after iterators / callbacks
  
    # Accept or reject the step
    if integrator.iter > 0
      if ((integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive) && !integrator.force_stepfail
        integrator.success_iter += 1
        StochasticDiffEq.apply_step!(integrator)
      elseif integrator.opts.adaptive && !integrator.accept_step
        if integrator.isout
          integrator.dtnew = integrator.dt*integrator.opts.qmin
        elseif !integrator.force_stepfail
          integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
        end
        StochasticDiffEq.choose_algorithm!(integrator,integrator.cache)
        StochasticDiffEq.fix_dtnew_at_bounds!(integrator)
        StochasticDiffEq.modify_dtnew_for_tstops!(integrator)
        reject_step!(integrator.W,integrator.dtnew)
        integrator.dt = integrator.dtnew
        integrator.sqdt = sqrt(abs(integrator.dt))
      end
    end
  
    integrator.iter += 1
    integrator.force_stepfail = false
  end

@inline function loopfooter!(integrator::SDDEIntegrator)
  ttmp = integrator.t + integrator.dt
  if integrator.force_stepfail
    if integrator.opts.adaptive
      integrator.dtnew = integrator.dt/integrator.opts.failfactor
    elseif integrator.last_stepfail
      return
    end
    integrator.last_stepfail = true
    integrator.accept_step = false
  elseif integrator.opts.adaptive
    @fastmath integrator.q11 = integrator.EEst^integrator.opts.beta1
    @fastmath integrator.q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    @fastmath integrator.q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma))
    @fastmath integrator.dtnew = integrator.dt/integrator.q
    integrator.isout = integrator.opts.isoutofdomain(integrator.u,integrator.p,ttmp)
    integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0) || (integrator.opts.force_dtmin && integrator.dt <= integrator.opts.dtmin)
    if integrator.accept_step # Accepted
      integrator.last_stepfail = false
      integrator.tprev = integrator.t
      if typeof(integrator.t)<:AbstractFloat && !isempty(integrator.opts.tstops)
        tstop = integrator.tdir * top(integrator.opts.tstops)
        @fastmath abs(ttmp - tstop) < 10eps(integrator.t) ? (integrator.t = tstop) : (integrator.t = ttmp)
      else
        integrator.t = ttmp
      end
      StochasticDiffEq.calc_dt_propose!(integrator)
      StochasticDiffEq.handle_callbacks!(integrator)
    end
  else # Non adaptive
    integrator.tprev = integrator.t
    if typeof(integrator.t)<:AbstractFloat && !isempty(integrator.opts.tstops)
      tstop = integrator.tdir * top(integrator.opts.tstops)
      # For some reason 10eps(integrator.t) is slow here
      # TODO: Allow higher precision but profile
      @fastmath abs(ttmp - tstop) < 10eps(max(integrator.t,tstop)) ? (integrator.t = tstop) : (integrator.t = ttmp)
    else
      integrator.t = ttmp
    end
    integrator.last_stepfail = false
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    StochasticDiffEq.handle_callbacks!(integrator)
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    @logmsg(-1,
    integrator.opts.progress_name,
    _id = :StochasticDelayDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress=integrator.t/integrator.sol.prob.tspan[2])
  end
end


@inline function postamble!(integrator::SDDEIntegrator)
  StochasticDiffEq.solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  if integrator.opts.progress
    @logmsg(-1,
    integrator.opts.progress_name,
    _id = :StochasticDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress="done")
  end
  return nothing
end

@inline function DiffEqBase.savevalues!(integrator::SDDEIntegrator,force_save=false)::Tuple{Bool,Bool}
  saved, savedexactly = false, false
  !integrator.opts.save_on && return saved, savedexactly
  tdir_t = integrator.tdir * integrator.t
  while !isempty(integrator.opts.saveat) && top(integrator.opts.saveat) <= tdir_t # Perform saveat
    integrator.saveiter += 1; saved = true
    curt = integrator.tdir * pop!(integrator.opts.saveat)
    if curt!=integrator.t # If <t, interpolate
      Θ = (curt - integrator.tprev)/integrator.dt
      val = StochasticDiffEq.sde_interpolant(Θ,integrator,integrator.opts.save_idxs,Val{0}) # out of place, but force copy later
      save_val = val
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,save_val,Val{false})
      if typeof(integrator.alg) <: StochasticDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      savedexactly = true
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      if integrator.opts.save_idxs === nothing
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      else
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
      end
      if typeof(integrator.alg) <: Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    end
  end
  if force_save || integrator.opts.save_everystep
    integrator.saveiter += 1; saved, savedexactly = true, true
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if typeof(integrator.alg) <: Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
  return saved, savedexactly
end