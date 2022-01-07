function DiffEqBase.__solve(prob::AbstractSDDEProblem, # TODO: DiffEqBase.AbstractSDDEProblem
    alg::AbstractSDEAlgorithm, args...; # TODO: Method of steps???
    kwargs...)
    integrator = DiffEqBase.__init(prob, alg, args...; kwargs...)
    DiffEqBase.solve!(integrator)
    integrator.sol
end

function DiffEqBase.__init(prob::AbstractSDDEProblem,# TODO DiffEqBasee.AbstractSDDEProblem
    alg::Union{AbstractRODEAlgorithm,AbstractSDEAlgorithm},
    timeseries_init = typeof(prob.u0)[],
    ts_init = eltype(prob.tspan)[],
    ks_init = nothing,
    recompile::Type{Val{recompile_flag}} = Val{true};
    saveat = eltype(prob.tspan)[],
    tstops = eltype(prob.tspan)[],
    d_discontinuities = Discontinuity{eltype(prob.tspan),Rational{Int}}[],
    save_idxs = nothing,
    save_everystep = isempty(saveat),
    save_noise = save_everystep && (typeof(prob.f) <: Tuple ?
                 DiffEqBase.has_analytic(prob.f[1]) : DiffEqBase.has_analytic(prob.f)),
    save_on = true,
    save_start = save_everystep || isempty(saveat) || typeof(saveat) <: Number ? true : prob.tspan[1] in saveat,
    save_end = nothing,
    callback = nothing,
    dense = save_everystep && isempty(saveat),
    calck = (!isempty(setdiff(saveat, tstops)) || dense),
    dt = eltype(prob.tspan)(0),
    adaptive = StochasticDiffEq.isadaptive(getalg(alg)),
    gamma = 9 // 10, # TODO gamma_default(alg.alg) ?
    abstol = nothing,
    reltol = nothing,
    qmax = StochasticDiffEq.qmax_default(getalg(alg)),
    qmin = StochasticDiffEq.qmin_default(getalg(alg)),
    qsteady_min = StochasticDiffEq.qsteady_min_default(alg),
    qsteady_max = StochasticDiffEq.qsteady_max_default(alg),
    qoldinit = 1 // 10^4, fullnormalize = true,
    controller = nothing,
    failfactor = 2,
    beta2 = nothing,
    beta1 = nothing,
    delta = StochasticDiffEq.delta_default(getalg(alg)),
    maxiters = adaptive ? 1000000 : typemax(Int),
    dtmax = eltype(prob.tspan)((prob.tspan[end] - prob.tspan[1])),
    dtmin = typeof(one(eltype(prob.tspan))) <: AbstractFloat ? eps(eltype(prob.tspan)) :
            typeof(one(eltype(prob.tspan))) <: Integer ? 0 :
            eltype(prob.tspan)(1 // 10^(10)),
    internalnorm = DiffEqBase.ODE_DEFAULT_NORM,
    isoutofdomain = DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN,
    unstable_check = DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK,
    verbose = true,force_dtmin = false,
    timeseries_errors = true, dense_errors = false,
    advance_to_tstop = false,stop_at_next_tstop = false,
    initialize_save = true,
    progress = false,progress_steps = 1000,progress_name = "SDDE",
    progress_message = DiffEqBase.ODE_DEFAULT_PROG_MESSAGE,
    userdata = nothing,
    initialize_integrator = true,
    seed = UInt64(0), alias_u0 = false,
    #  Keywords for Delay problems (from DDE)
    discontinuity_interp_points::Int = 10,
    discontinuity_abstol = eltype(prob.tspan)(1 // Int64(10)^12),
    discontinuity_reltol = 0, kwargs...) where recompile_flag

    # alg = getalg(alg0);
    if typeof(prob.f) <: Tuple
        if any(mm != I for mm in prob.f.mass_matrix)
            error("This solver is not able to use mass matrices.")
        end
    elseif prob.f.mass_matrix != I && !alg_mass_matrix_compatible(getalg(alg))
        error("This solver is not able to use mass matrices.")
    end

    if !isempty(saveat) && dense
        @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
    end

    if typeof(prob.noise) <: NoiseProcess && prob.noise.bridge === nothing && adaptive
        error("Bridge function must be given for adaptivity. Either declare this function in noise process or set adaptive=false")
    end

    #if !StochasticDiffEq.alg_compatible(prob, getalg(alg))
    #    error("The algorithm is not compatible with the chosen noise type. Please see the documentation on the solver methods")
    #end

    if haskey(kwargs, :initial_order)
        @warn "initial_order has been deprecated. Please specify order_discontinuity_t0 in the DDEProblem instead."
        order_discontinuity_t0 = kwargs[:initial_order]
    else
        order_discontinuity_t0 = prob.order_discontinuity_t0
    end

    if haskey(kwargs, :minimal_solution)
        @warn "minimal_solution is ignored"
    end
    progress && @logmsg(-1,progress_name,_id = _id = :StochasticDiffEq,progress = 0)

    tType = eltype(prob.tspan)
    noise = prob.noise
    tspan = prob.tspan
    tdir = sign(tspan[end] - tspan[1])

    t = tspan[1]

    if !adaptive && iszero(dt) && isempty(tstops)
        error("Fixed timestep methods require a choice of dt or choosing the tstops")
    end
    if has_dependent_lags(prob)
        @warn "Discontinuities are not handled for dependent lags. Please consider it when evaluating the soution!"
    end

    @unpack f, g, u0, h, tspan, p, constant_lags, dependent_lags, neutral = prob
    # g = prob isa AbstractSDDEProblem ? prob.g : nothing # TODO DiffEqBase.AbstractSDDEProblem

    # no fixed-point iterations for stochastic delay problems,
    # and thus `dtmax` should match minimal lag
    if  has_constant_lags(prob) # && isconstrained(alg)
        dtmax = tdir * min(abs(dtmax), minimum(abs, constant_lags))
    end

    u, uprev = u_uprev(prob.u0, alias_u0 = alias_u0)

    uType = typeof(u)
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    ks = Vector{uType}(undef, 0)

    uEltypeNoUnits = recursive_unitless_eltype(u)
    tTypeNoUnits   = typeof(one(tType))

    if abstol === nothing
        if uBottomEltypeNoUnits == uBottomEltype
            abstol_internal = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
        else
            abstol_internal = real.(oneunit.(u) .* 1 // 10^2)
        end
    else
        abstol_internal = real.(abstol)
    end

    if reltol === nothing
        if uBottomEltypeNoUnits == uBottomEltype
            reltol_internal = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
        else
            reltol_internal = real.(oneunit.(u) .* 1 // 10^2)
        end
    else
        reltol_internal = real.(reltol)
    end

    if isinplace(prob) && typeof(u) <: AbstractArray && eltype(u) <: Number && uBottomEltypeNoUnits == uBottomEltype # Could this be more efficient for other arrays?
        if !(typeof(u) <: ArrayPartition)
            rate_prototype = recursivecopy(u)
        else
            rate_prototype = similar(u, typeof.(oneunit.(recursive_bottom_eltype.(u.x)) ./ oneunit(tType))...)
        end
    else
        if uBottomEltypeNoUnits == uBottomEltype
            rate_prototype = u
        else # has units!
            rate_prototype = u / oneunit(tType)
        end
    end
    rateType = typeof(rate_prototype) ## Can be different if united

    if StochasticDiffEq.is_diagonal_noise(prob)
        noise_rate_prototype = rate_prototype
    else
        if prob isa AbstractSDDEProblem # TODO DiffEqBase.AbstractSDDEProblem
            noise_rate_prototype = copy(prob.noise_rate_prototype)
        else
            noise_rate_prototype = copy(prob.rand_prototype)
        end
    end

    #= TODO: Jump handling
    if typeof(_prob) <: JumpProblem && _prob.regular_jump !== nothing

      if !isnothing(_prob.regular_jump.mark_dist) == nothing # https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/250
        error("Mark distributions are currently not supported in SimpleTauLeaping")
      end

      jump_prototype = zeros(_prob.regular_jump.numjumps)
      c = _prob.regular_jump.c

      if isinplace(_prob.regular_jump)
        rate_constants = zeros(_prob.regular_jump.numjumps)
        _prob.regular_jump.rate(rate_constants,u./u,prob.p,tspan[1])
        P = CompoundPoissonProcess!(_prob.regular_jump.rate,t,jump_prototype,
                                    computerates = !alg_control_rate(alg) || !adaptive,
                                    save_everystep=save_noise,
                                    rng = Xorshifts.Xoroshiro128Plus(_seed))
        alg_control_rate(alg) && adaptive && P.cache.rate(P.cache.currate,u,p,tspan[1])
      else
        rate_constants = _prob.regular_jump.rate(u./u,prob.p,tspan[1])
        P = CompoundPoissonProcess(_prob.regular_jump.rate,t,jump_prototype,
                                   save_everystep=save_noise,
                                   computerates = !alg_control_rate(alg) || !adaptive,
                                   rng = Xorshifts.Xoroshiro128Plus(_seed))
        alg_control_rate(alg) && adaptive && (P.cache.currate = P.cache.rate(u,p,tspan[1]))
      end

    else
    =#
      jump_prototype = nothing
      c = nothing
      P = nothing
      rate_constants = nothing
    #end

    # tstops_internal, saveat_internal, d_discontinuities_internal =
    #   tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan) # TODO add delays to discontinuities
    # retrieve time stops, time points at which solutions is saved, and discontinuities
    maximum_order = StochasticDiffEq.alg_order(getalg(alg))
    tstops_internal, saveat_internal, d_discontinuities_internal =
      tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan,
                                                order_discontinuity_t0, maximum_order,
                                                constant_lags, neutral)

    tracked_discontinuities = Discontinuity{tType,Rational{Int}}[]
    if order_discontinuity_t0 ≤ maximum_order
        push!(tracked_discontinuities, Discontinuity(tdir * t, Rational{Int}(order_discontinuity_t0)))
    end

    callbacks_internal = CallbackSet(callback)

    max_len_cb = DiffEqBase.max_vector_callback_length(callbacks_internal)
    if max_len_cb isa DiffEqBase.VectorContinuousCallback
        callback_cache = DiffEqBase.CallbackCache(max_len_cb.len, uBottomEltype, uBottomEltype)
    else
        callback_cache = nothing
    end

    QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits

    if !(uType <: AbstractArray)
        rand_prototype = zero(u / u) # Strip units and type info
        randType = typeof(rand_prototype)
    else
        randElType = uBottomEltypeNoUnits # Strip units and type info
        if StochasticDiffEq.is_diagonal_noise(prob)
            if typeof(u) <: SArray
                rand_prototype = zero(u) # TODO: Array{randElType} for units
            else
                rand_prototype = (u .- u) ./ sqrt(oneunit(t))
            end
        elseif prob isa AbstractSDDEProblem # TODO DiffEqBase.AbstractSDDEProblem
            rand_prototype = false .* noise_rate_prototype[1,:]
        else
            rand_prototype = copy(prob.rand_prototype)
        end
        randType = typeof(rand_prototype) # Strip units and type info
    end

    _seed = iszero(seed) ? (iszero(prob.seed) ? rand(UInt64) : prob.seed) : seed

    if prob.noise === nothing
        rswm = StochasticDiffEq.isadaptive(getalg(alg)) ? RSWM(adaptivealg = :RSwM3) : RSWM(adaptivealg = :RSwM1)
        if isinplace(prob)
            if StochasticDiffEq.alg_needs_extra_process(getalg(alg))
                W = WienerProcess!(t, rand_prototype, rand_prototype,
                             save_everystep = save_noise,
                             rswm = rswm,
                             rng = Xorshifts.Xoroshiro128Plus(_seed))
            else
                W = WienerProcess!(t, rand_prototype,
                             save_everystep = save_noise,
                             rswm = rswm,
                             rng = Xorshifts.Xoroshiro128Plus(_seed))
            end
        else
            if StochasticDiffEq.alg_needs_extra_process(getalg(alg))
                W = WienerProcess(t, rand_prototype, rand_prototype,
                             save_everystep = save_noise,
                             rswm = rswm,
                             rng = Xorshifts.Xoroshiro128Plus(_seed))
            else
                W = WienerProcess(t, rand_prototype,
                             save_everystep = save_noise,
                             rswm = rswm,
                             rng = Xorshifts.Xoroshiro128Plus(_seed))
            end
        end
    else
        W = prob.noise
        if W.reset
            if W.t[end] != t
                reinit!(W, t)
            end
        # Reseed
            if typeof(W) <: NoiseProcess && W.reseed
                Random.seed!(W.rng, _seed)
            end
        elseif W.t[end] != t
            error("Starting time in the noise process is not the starting time of the simulation. The noise process should be re-initialized for repeated use")
        end
    end

    ts, timeseries, saveiter = solution_arrays(u, tspan, rate_prototype,
                    timeseries_init = timeseries_init, ts_init = ts_init, save_idxs = save_idxs, save_start = save_start)

    if !adaptive && save_everystep && tspan[2] - tspan[1] != Inf
        iszero(dt) ? steps = length(tstops) : steps = ceil(Int, internalnorm((tspan[2] - tspan[1]) / dt, tspan[1]))
        sizehint!(timeseries, steps + 1)
        sizehint!(ts, steps + 1)
    elseif save_everystep
        sizehint!(timeseries, 50)
        sizehint!(ts, 50)
    elseif !isempty(saveat_internal)
        sizehint!(timeseries, length(saveat_internal) + 1)
        sizehint!(ts, length(saveat_internal) + 1)
    else
        sizehint!(timeseries, 2)
        sizehint!(ts, 2)
    end

    alg_choice = Int[]
    if save_start && typeof(getalg(alg)) <: StochasticDiffEq.StochasticDiffEqCompositeAlgorithm
        copyat_or_push!(alg_choice, 1, 1)
    end

    # create a history function
    history = build_history_function(prob, alg, reltol_internal,
                                    rate_prototype, noise_rate_prototype, jump_prototype, W, _seed, dense;
                                    dt = dt, adaptive = adaptive,
                                    internalnorm = internalnorm)
    f_with_history, g_with_history = wrap_functions_and_history(f, g, history)

    sde_integrator = history.integrator;

    cache = StochasticDiffEq.alg_cache(getalg(alg), prob, u, W.dW, W.dZ, p, rate_prototype, noise_rate_prototype, jump_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, f_with_history, t, dt, Val{isinplace(prob)})

    # id = StochasticDiffEq.LinearInterpolationData(timeseries,ts)
    id = StochasticDiffEq.LinearInterpolationData(sde_integrator.sol.u, sde_integrator.sol.t)

    save_end_user = save_end
    save_end = save_end === nothing ? save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[2] in saveat : save_end

    # Setting up the step size controller
    if (beta1 !== nothing || beta2 !== nothing) && controller !== nothing
      throw(ArgumentError(
        "Setting both the legacy PID parameters `beta1, beta2 = $((beta1, beta2))` and the `controller = $controller` is not allowed."))
    end

    if (beta1 !== nothing || beta2 !== nothing)
      message = "Providing the legacy PID parameters `beta1, beta2` is deprecated. Use the keyword argument `controller` instead."
      Base.depwarn(message, :init)
      Base.depwarn(message, :solve)
    end

  if controller === nothing
    controller = StochasticDiffEq.default_controller(_alg, cache, convert(QT,qoldinit),
                                    beta1 === nothing ? nothing : convert(QT,beta1),
                                    beta2 === nothing ? nothing : convert(QT,beta2))
  end

    opts = StochasticDiffEq.SDEOptions(maxiters, save_everystep,
                      adaptive, abstol_internal,
                      reltol_internal,
                      QT(gamma),
                      QT(qmax), QT(qmin),
                      QT(qsteady_max),QT(qsteady_min),
                      QT(qoldinit),
                      QT(failfactor),
                      tType(dtmax), tType(dtmin),
                      controller,
                      internalnorm, save_idxs,
                      tstops_internal, saveat_internal,
                      d_discontinuities_internal,
                      tstops, saveat, d_discontinuities,
                      userdata,
                      progress, progress_steps,
                      progress_name, progress_message,
                      timeseries_errors, dense_errors,
                      convert.(uBottomEltypeNoUnits, delta),
                      dense, save_on, save_start, save_end, save_end_user, save_noise,
                      callbacks_internal, isoutofdomain, unstable_check,
                      verbose, calck, force_dtmin,
                      advance_to_tstop, stop_at_next_tstop)

    destats = DiffEqBase.DEStats(0)
    if typeof(getalg(alg)) <: StochasticDiffEq.StochasticDiffEqCompositeAlgorithm
      # TODO: DISCONNECT!!!!
        sol =  DiffEqBase.build_solution(prob, alg, sde_integrator.sol.t, sde_integrator.sol.u, W = W,
                                      destats = destats,
                                      calculate_error = false, alg_choice = alg_choice,
                                      interp = id, dense = dense, seed = _seed)
    # separate statistics of the integrator and the history
    else
      # TODO: DISCONNECT!!!!
        sol = DiffEqBase.build_solution(prob, alg, sde_integrator.sol.t, sde_integrator.sol.u, W = W,
                                      destats = destats,
                                      calculate_error = false,
                                      interp = id, dense = dense, seed = _seed)
    # separate statistics of the integrator and the history
    end

    if recompile_flag == true
        FType = typeof(f_with_history)
        GType = typeof(g_with_history)
        SolType = typeof(sol)
        cacheType = typeof(cache)
    else
        FType = Function
        GType = Function
        SolType = DiffEqBase.AbstractRODESolution
        cacheType = StochasticDiffEq.StochasticDiffEqCache
    end

    tprev = t
    dtcache = tType(dt)
    iter = 0
    u_modified = false
    eigen_est = one(uBottomEltypeNoUnits) / oneunit(tType) # rate/state = (state/time)/state = 1/t units
    EEst = tTypeNoUnits(1)
    just_hit_tstop = false
    do_error_check = true
    isout = false
    accept_step = false
    force_stepfail = false
    last_stepfail = false
    event_last_time = 0
    vector_event_last_time = 1
    last_event_error = zero(uBottomEltypeNoUnits)
    dtchangeable = true
    q11 = tTypeNoUnits(1)
    success_iter = 0
    q = tTypeNoUnits(1)

    integrator = SDDEIntegrator{typeof(getalg(alg)),isinplace(prob),uType,
                                uBottomEltype,tType,typeof(p),typeof(eigen_est),
                                QT,uEltypeNoUnits,typeof(W),typeof(P),
                                rateType,typeof(sol),typeof(cache),FType,GType,typeof(c),
                                typeof(opts),typeof(noise),typeof(last_event_error),
                                typeof(callback_cache),typeof(history),
                                typeof(sde_integrator)}(f_with_history,
                                g_with_history, c, noise, uprev, tprev,
                      order_discontinuity_t0, tracked_discontinuities,
                      t, u, p, tType(dt), tType(dt), tType(dt), dtcache, tspan[2], tdir,
                      just_hit_tstop, do_error_check,
                      isout, event_last_time, vector_event_last_time, last_event_error, accept_step,
                      last_stepfail, force_stepfail, dtchangeable, u_modified, saveiter, getalg(alg), sol,
                      cache, callback_cache, tType(dt), W, P,
                      opts, iter, success_iter, eigen_est, EEst, q, QT(qoldinit), q11, history, destats, sde_integrator)

    if initialize_integrator
        StochasticDiffEq.initialize_callbacks!(integrator, initialize_save)
        initialize!(integrator, integrator.cache)

        save_start && typeof(alg) <: StochasticDiffEq.StochasticDiffEqCompositeAlgorithm && copyat_or_push!(alg_choice, 1, integrator.cache.current)
    end

    StochasticDiffEq.handle_dt!(integrator)

    ## Modify the first dt for tstops
    StochasticDiffEq.modify_dt_for_tstops!(integrator)
    ### Needs to be done before first rand
    integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))

    integrator.W.dt = integrator.dt
    DiffEqNoiseProcess.setup_next_step!(integrator)

    integrator
end


function DiffEqBase.solve!(integrator::SDDEIntegrator)
    @inbounds while !isempty(integrator.opts.tstops)
        while integrator.tdir * integrator.t < first(integrator.opts.tstops)
            loopheader!(integrator)
            if DiffEqBase.check_error!(integrator) != :Success
                return integrator.sol
            end
            StochasticDiffEq.perform_step!(integrator, integrator.cache)
            loopfooter!(integrator)
            if isempty(integrator.opts.tstops)
                break
            end
        end
        @inbounds handle_tstop!(integrator)
    end
    postamble!(integrator)

    f = typeof(integrator.sol.prob.f) <: Tuple ? integrator.sol.prob.f[1] : integrator.sol.prob.f

    if DiffEqBase.has_analytic(f)
        DiffEqBase.calculate_solution_errors!(integrator.sol;timeseries_errors = integrator.opts.timeseries_errors,dense_errors = integrator.opts.dense_errors)
    end
    if integrator.sol.retcode != :Default
        return integrator.sol
    end
    integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, :Success)
end


# function StochasticDiffEq.tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan)
function tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan, order_discontinuity_t0, alg_maximum_order, constant_lags, neutral)

    tType = eltype(tspan)
    t0, tf = tspan
    tdir = sign(tf - t0)
    tdir_t0 = tdir * t0
    tdir_tf = tdir * tf

  # add discontinuities propagated from initial discontinuity
    add_propagated_constant_lags = order_discontinuity_t0 ≤ alg_maximum_order && constant_lags !== nothing && !isempty(constant_lags)

    if add_propagated_constant_lags
        maxlag = tdir_tf - tdir_t0
        next_order = neutral ? order_discontinuity_t0 : order_discontinuity_t0 + 1//2
    end

  # time stops
    tstops_internal = BinaryMinHeap{tType}()
    if isempty(d_discontinuities) && !add_propagated_constant_lags && isempty(tstops) # TODO: Specialize more
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

  # add propagated discontinuities
        if add_propagated_constant_lags
            for lag in constant_lags
                if tdir * lag < maxlag
                    push!(tstops_internal, tdir * (t0 + lag))
                end
            end
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
    d_discontinuities_internal = BinaryMinHeap{Discontinuity{tType,Rational{Int}}}()
    if add_propagated_constant_lags
        sizehint!(d_discontinuities_internal.valtree, length(d_discontinuities) + length(constant_lags))
    else
        sizehint!(d_discontinuities_internal.valtree, length(d_discontinuities))
    end

    for d in d_discontinuities
        tdir_t = tdir * d.t

        if tdir_t0 < tdir_t < tdir_tf && d.order ≤ alg_maximum_order + 1
            push!(d_discontinuities_internal, Discontinuity{tType,Rational{Int}}(tdir_t, d.order))
        end
    end

    if add_propagated_constant_lags
        for lag in constant_lags
            if tdir * lag < maxlag
                push!(d_discontinuities_internal, Discontinuity{tType,Rational{Int}}(tdir * (t0 + lag), next_order))
            end
        end
    end

    tstops_internal, saveat_internal, d_discontinuities_internal
end
