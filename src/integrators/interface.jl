@inline function loopheader!(integrator::SDDEIntegrator)
    # Apply right after iterators / callbacks

    # Accept or reject the step
    if integrator.iter > 0
        if ((integrator.opts.adaptive && integrator.accept_step) ||
            !integrator.opts.adaptive) && !integrator.force_stepfail
            integrator.success_iter += 1
            # StochasticDiffEq.apply_step!(integrator)
            StochasticDiffEq.apply_step!(integrator)
        elseif integrator.opts.adaptive && !integrator.accept_step
            if integrator.isout
                integrator.dtnew = integrator.dt * integrator.opts.qmin
            elseif !integrator.force_stepfail
                step_reject_controller!(integrator, getalg(integrator.alg))
            end
            StochasticDiffEq.choose_algorithm!(integrator, integrator.cache)
            StochasticDiffEq.fix_dtnew_at_bounds!(integrator)
            StochasticDiffEq.modify_dtnew_for_tstops!(integrator)
            reject_step!(integrator, integrator.dtnew)
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
            integrator.dtnew = integrator.dt / integrator.opts.failfactor
        elseif integrator.last_stepfail
            return
        end
        integrator.last_stepfail = true
        integrator.accept_step = false
    elseif integrator.opts.adaptive
        stepsize_controller!(integrator, getalg(integrator.alg))
        integrator.isout = integrator.opts.isoutofdomain(integrator.u, integrator.p, ttmp)
        integrator.accept_step = (!integrator.isout && accept_step_controller(integrator,
                                                         integrator.opts.controller)) ||
                                 (integrator.opts.force_dtmin &&
                                  integrator.dt <= integrator.opts.dtmin)
        if integrator.accept_step # Accepted
            step_accept_controller!(integrator, getalg(integrator.alg))
            integrator.last_stepfail = false
            integrator.tprev = integrator.t
            if integrator.t isa AbstractFloat && !isempty(integrator.opts.tstops)
                tstop = integrator.tdir * first(integrator.opts.tstops)
                @fastmath abs(ttmp - tstop) < 10eps(integrator.t) ? (integrator.t = tstop) :
                          (integrator.t = ttmp)
            else
                integrator.t = ttmp
            end
            StochasticDiffEq.calc_dt_propose!(integrator)
            StochasticDiffEq.handle_callbacks!(integrator)
        end
    else # Non adaptive
        integrator.tprev = integrator.t
        if integrator.t isa AbstractFloat && !isempty(integrator.opts.tstops)
            tstop = integrator.tdir * first(integrator.opts.tstops)
            # For some reason 10eps(integrator.t) is slow here
            # TODO: Allow higher precision but profile
            @fastmath abs(ttmp - tstop) < 10eps(max(integrator.t, tstop)) ?
                      (integrator.t = tstop) : (integrator.t = ttmp)
        else
            integrator.t = ttmp
        end
        integrator.last_stepfail = false
        integrator.accept_step = true
        integrator.dtpropose = integrator.dt
        StochasticDiffEq.handle_callbacks!(integrator)
    end
    if integrator.opts.progress && integrator.iter % integrator.opts.progress_steps == 0
        @logmsg(-1,
                integrator.opts.progress_name,
                _id=:StochasticDelayDiffEq,
                message=integrator.opts.progress_message(integrator.dt, integrator.u,
                                                         integrator.p, integrator.t),
                progress=integrator.t / integrator.sol.prob.tspan[2])
    end
end

@inline function postamble!(integrator::SDDEIntegrator)
    StochasticDiffEq.solution_endpoint_match_cur_integrator!(integrator)
    resize!(integrator.sol.t, integrator.saveiter)
    resize!(integrator.sol.u, integrator.saveiter)
    if integrator.opts.progress
        @logmsg(-1,
                integrator.opts.progress_name,
                _id=:StochasticDiffEq,
                message=integrator.opts.progress_message(integrator.dt, integrator.u,
                                                         integrator.p, integrator.t),
                progress="done")
    end
    return nothing
end

function DiffEqBase.auto_dt_reset!(integrator::SDDEIntegrator)
    # @unpack f,g, u, t, tdir, opts, sol, stats = integrator
    @unpack f, g, u, t, tdir, opts, sol = integrator
    @unpack prob = sol
    @unpack abstol, reltol, internalnorm = opts

    # determine maximal time step
    if has_constant_lags(prob)
        dtmax = tdir * min(abs(opts.dtmax), minimum(abs, prob.constant_lags))
    else
        dtmax = opts.dtmax
    end

    # determine initial time step
    sde_prob = SDEProblem(f, g, prob.u0, prob.tspan, prob.p;
                          noise_rate_prototype = prob.noise_rate_prototype,
                          noise = prob.noise,
                          seed = prob.seed,
                          prob.kwargs...)
    integrator.dt = StochasticDiffEq.sde_determine_initdt(integrator.u, integrator.t,
                                                          integrator.tdir,
                                                          integrator.opts.dtmax,
                                                          integrator.opts.abstol,
                                                          integrator.opts.reltol,
                                                          integrator.opts.internalnorm,
                                                          sde_prob,
                                                          StochasticDiffEq.get_current_alg_order(getalg(integrator.alg),
                                                                                                 integrator.cache),
                                                          integrator)
    # update statistics
    # destats.nf += 2
    # destats.nf2 += 2

    nothing
end

DiffEqBase.has_reinit(integrator::SDDEIntegrator) = false # TODO true - syncronize with DDE!
function DiffEqBase.reinit!(integrator::SDDEIntegrator, u0 = integrator.sol.prob.u0;
                            t0 = integrator.sol.prob.tspan[1],
                            tf = integrator.sol.prob.tspan[2],
                            erase_sol = true,
                            tstops = integrator.opts.tstops_cache,
                            saveat = integrator.opts.saveat_cache,
                            d_discontinuities = integrator.opts.d_discontinuities_cache,
                            order_discontinuity_t0 = t0 == integrator.sol.prob.tspan[1] &&
                                u0 == integrator.sol.prob.u0 ?
                                                     integrator.order_discontinuity_t0 : 0,
                            reinit_cache = true, reinit_callbacks = true,
                            initialize_save = true,
                            reset_dt = (integrator.dtcache == zero(integrator.dt)) &&
                                integrator.opts.adaptive)
    if DiffEqBase.isinplace(integrator.sol.prob)
        recursivecopy!(integrator.u, u0)
        recursivecopy!(integrator.uprev, integrator.u)
    else
        integrator.u = u0
        integrator.uprev = integrator.u
    end

    integrator.t = t0
    integrator.tprev = t0

    tType = typeof(integrator.t)
    maximum_order = StochasticDiffEq.alg_order(getalg(integrator.alg))
    tstops_internal, saveat_internal, d_discontinuities_internal = tstop_saveat_disc_handling(tstops,
                                                                                              saveat,
                                                                                              d_discontinuities,
                                                                                              (tType(t0),
                                                                                               tType(tf)),
                                                                                              order_discontinuity_t0,
                                                                                              maximum_order,
                                                                                              integrator.sol.prob.constant_lags,
                                                                                              integrator.sol.prob.neutral)

    integrator.opts.tstops = tstops_internal
    integrator.opts.saveat = saveat_internal
    integrator.opts.d_discontinuities = d_discontinuities_internal
    integrator.order_discontinuity_t0 = order_discontinuity_t0

    if erase_sol
        if integrator.opts.save_start
            resize_start = 1
        else
            resize_start = 0
        end
        resize!(integrator.sol.u, resize_start)
        resize!(integrator.sol.t, resize_start)
        if integrator.sol.u_analytic !== nothing
            resize!(integrator.sol.u_analytic, 0)
        end
        if typeof(getalg(integrator.alg)) <:
           StochasticDiffEq.StochasticDiffEqCompositeAlgorithm
            resize!(integrator.sol.alg_choice, resize_start)
        end
        integrator.saveiter = resize_start

        if order_discontinuity_t0 ≤ maximum_order
            resize!(integrator.tracked_discontinuities, 1)
            integrator.tracked_discontinuities[1] = Discontinuity(integrator.tdir *
                                                                  integrator.t,
                                                                  Rational{Int}(order_discontinuity_t0))
        else
            resize!(integrator.tracked_discontinuities, 0)
        end
    end

    integrator.iter = 0
    integrator.success_iter = 0

    # full re-initialize the PI in timestepping
    integrator.qold = integrator.opts.qoldinit
    integrator.q11 = typeof(integrator.t)(1)
    integrator.u_modified = false

    if reset_dt
        DiffEqBase.auto_dt_reset!(integrator)
    end

    if reinit_callbacks
        StochasticDiffEq.initialize_callbacks!(integrator, initialize_save)
    end

    if reinit_cache
        StochasticDiffEq.initialize!(integrator.integrator, integrator.cache)
        # StochasticDiffEq.initialize!(integrator,integrator.cache)
    end

    reinit!(integrator.W, integrator.dt)
    nothing
end

@inline function DiffEqBase.get_du(integrator::SDDEIntegrator)
    (integrator.u - integrator.uprev) / integrator.dt
end

@inline function DiffEqBase.get_du!(out, integrator::SDDEIntegrator)
    DiffEqBase.@.. out = (integrator.u - integrator.uprev) / integrator.dt
end

@inline function DiffEqBase.add_tstop!(integrator::SDDEIntegrator, t)
    t < integrator.t &&
        error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
    push!(integrator.opts.tstops, integrator.tdir * t)
end

function DiffEqBase.change_t_via_interpolation!(integrator::SDDEIntegrator, t,
                                                modify_save_endpoint::Type{Val{T}} = Val{
                                                                                         false
                                                                                         }) where {
                                                                                                   T
                                                                                                   }
    StochasticDiffEq.change_t_via_interpolation!(integrator, t, modify_save_endpoint)
end

# update integrator when u is modified by callbacks
function StochasticDiffEq.handle_callback_modifiers!(integrator::SDDEIntegrator)
    # copied from StochasticDiffEq
    if integrator.P !== nothing && integrator.opts.adaptive
        if integrator.cache isa StochasticDiffEqMutableCache
            oldrate = integrator.P.cache.currate
            integrator.P.cache.rate(oldrate, integrator.u, integrator.p, integrator.t)
        else
            integrator.P.cache.currate = integrator.P.cache.rate(integrator.u, integrator.p,
                                                                 integrator.t)
        end
    end

    # update heap of discontinuities
    # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
    push!(integrator.opts.d_discontinuities,
          Discontinuity(integrator.tdir * integrator.t, 0 // 1))
end

function DiffEqBase.u_modified!(integrator::SDDEIntegrator, bool::Bool)
    integrator.u_modified = bool
    nothing
end

DiffEqBase.get_proposed_dt(integrator::SDDEIntegrator) = integrator.dtpropose
function DiffEqBase.set_proposed_dt!(integrator::SDDEIntegrator, dt::Number)
    (integrator.dtpropose = dt; integrator.dtcache = dt)
end

function DiffEqBase.set_proposed_dt!(integrator::SDDEIntegrator,
                                     integrator2::SDDEIntegrator)
    integrator.dtpropose = integrator2.dtpropose
    integrator.qold = integrator2.qold
end

@inline function postamble!(integrator::HistorySDEIntegrator)
    if integrator.saveiter == 0 || integrator.sol.t[integrator.saveiter] != integrator.t
        integrator.saveiter += 1
        copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
        if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
        else
            copyat_or_push!(integrator.sol.u, integrator.saveiter,
                            integrator.u[integrator.opts.save_idxs], Val{false})
        end
    end

    resize!(integrator.sol.t, integrator.saveiter)
    resize!(integrator.sol.u, integrator.saveiter)
end

function DiffEqBase.postamble!(integrator::SDDEIntegrator)
    # clean up solution of the SDE integrator
    postamble!(integrator.integrator)

    # # clean solution of the SDDE integrator
    # StochasticDiffEq._postamble!(integrator)
end

@inline function DiffEqBase.savevalues!(integrator::SDDEIntegrator,
                                        force_save = false)::Tuple{Bool, Bool}
    saved, savedexactly = false, false
    !integrator.opts.save_on && return saved, savedexactly
    tdir_t = integrator.tdir * integrator.t
    while !isempty(integrator.opts.saveat) && first(integrator.opts.saveat) <= tdir_t # Perform saveat
        integrator.saveiter += 1
        saved = true
        curt = integrator.tdir * pop!(integrator.opts.saveat)
        if curt != integrator.t # If <t, interpolate
            Θ = (curt - integrator.tprev) / integrator.dt
            val = StochasticDiffEq.sde_interpolant(Θ, integrator, integrator.opts.save_idxs,
                                                   Val{0}) # out of place, but force copy later
            save_val = val
            copyat_or_push!(integrator.sol.t, integrator.saveiter, curt)
            copyat_or_push!(integrator.sol.u, integrator.saveiter, save_val, Val{false})
            if integrator.alg isa StochasticDiffEq.StochasticDiffEqCompositeAlgorithm
                copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter,
                                integrator.cache.current)
            end
        else # ==t, just save
            savedexactly = true
            copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
            if integrator.opts.save_idxs === nothing
                copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
            else
                copyat_or_push!(integrator.sol.u, integrator.saveiter,
                                integrator.u[integrator.opts.save_idxs], Val{false})
            end
            if integrator.alg isa
               Union{StochasticDiffEq.StochasticDiffEqCompositeAlgorithm,
                     StochasticDiffEq.StochasticDiffEqRODECompositeAlgorithm}
                copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter,
                                integrator.cache.current)
            end
        end
    end
    if force_save || integrator.opts.save_everystep
        integrator.saveiter += 1
        saved, savedexactly = true, true
        if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
        else
            copyat_or_push!(integrator.sol.u, integrator.saveiter,
                            integrator.u[integrator.opts.save_idxs], Val{false})
        end
        copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
        if integrator.alg isa
           Union{StochasticDiffEq.StochasticDiffEqCompositeAlgorithm,
                 StochasticDiffEq.StochasticDiffEqRODECompositeAlgorithm}
            copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter,
                            integrator.cache.current)
        end
    end
    return saved, savedexactly
end

@inline function DiffEqNoiseProcess.setup_next_step!(integrator::SDDEIntegrator)
    !isnothing(integrator.W) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.W, integrator.u, integrator.p)
    !isnothing(integrator.P) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.P, integrator.u, integrator.p)
end

@inline function DiffEqNoiseProcess.reject_step!(integrator::SDDEIntegrator,
                                                 dtnew = integrator.dtnew)
    !isnothing(integrator.W) &&
        reject_step!(integrator.W, dtnew, integrator.u, integrator.p)
    !isnothing(integrator.P) &&
        reject_step!(integrator.P, dtnew, integrator.u, integrator.p)
end

@inline function DiffEqNoiseProcess.accept_step!(integrator::SDDEIntegrator, setup)
    !isnothing(integrator.W) &&
        accept_step!(integrator.W, integrator.dt, integrator.u, integrator.p, setup)
    !isnothing(integrator.P) &&
        accept_step!(integrator.P, integrator.dt, integrator.u, integrator.p, setup)
end

@inline function DiffEqNoiseProcess.save_noise!(integrator::SDDEIntegrator)
    !isnothing(integrator.W) && DiffEqNoiseProcess.save_noise!(integrator.W)
    !isnothing(integrator.P) && DiffEqNoiseProcess.save_noise!(integrator.P)
end

@inline function DiffEqBase.get_tmp_cache(integrator::SDDEIntegrator)
    get_tmp_cache(integrator, integrator.alg, integrator.cache)
end
# avoid method ambiguity
for typ in (StochasticDiffEq.StochasticDiffEqAlgorithm,
            StochasticDiffEq.StochasticDiffEqNewtonAdaptiveAlgorithm)
    @eval @inline function DiffEqBase.get_tmp_cache(integrator::SDDEIntegrator, alg::$typ,
                                                    cache::StochasticDiffEq.StochasticDiffEqConstantCache)
        nothing
    end
end
@inline DiffEqBase.get_tmp_cache(integrator::SDDEIntegrator, alg, cache) = (cache.tmp,)
@inline function DiffEqBase.get_tmp_cache(integrator::SDDEIntegrator,
                                          alg::StochasticDiffEq.StochasticDiffEqNewtonAdaptiveAlgorithm,
                                          cache)
    (cache.nlsolver.tmp, cache.nlsolver.ztmp)
end
@inline function DiffEqBase.get_tmp_cache(integrator::SDDEIntegrator,
                                          alg::StochasticDiffEq.StochasticCompositeAlgorithm,
                                          cache)
    get_tmp_cache(integrator, alg.algs[1], cache.caches[1])
end
