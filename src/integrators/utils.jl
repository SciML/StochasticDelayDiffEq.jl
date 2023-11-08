@inline function handle_tstop!(integrator)
    tstops = integrator.opts.tstops
    if !isempty(tstops)
        tdir_t = integrator.tdir * integrator.t
        tdir_ts_top = first(tstops)
        if tdir_t == tdir_ts_top
            pop!(tstops)
            integrator.just_hit_tstop = true
        elseif tdir_t > tdir_ts_top
            if !integrator.dtchangeable
                change_t_via_interpolation!(integrator, integrator.tdir * pop!(tstops),
                                            Val{true})
                integrator.just_hit_tstop = true
            else
                error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
            end
        end
        isempty(integrator.opts.d_discontinuities) || handle_discontinuities!(integrator)
    end
end
#=
Dealing with discontinuities

If we hit a discontinuity (this is checked in `apply_step!`), then we remove the
discontinuity, additional discontinuities at the current time point (if present), and
maybe also discontinuities and time stops coming shortly after the current time point
in `handle_discontinuities!`. The order of the discontinuity at the current time point is
defined as the lowest order of all these discontinuities.

If the problem is not neutral, we will only add additional discontinuities if
this order is less or equal to the order of the algorithm in
`add_next_discontinuities!`. If we add discontinuities, we add discontinuities
of the next order caused by constant lags (these we can calculate explicitly and
just add them to `d_discontinuities` and `tstops`) and we add the current
discontinuity to `tracked_discontinuities` which is the array of old
discontinuities that are checked by a `DiscontinuityCallback` (if existent).
=#

# handle discontinuities at the current time point of the `integrator`
function handle_discontinuities!(integrator::SDDEIntegrator)
    # remove all discontinuities at current time point and calculate minimal order
    # of these discontinuities
    d = pop!(integrator.opts.d_discontinuities)
    order = d.order
    while !isempty(integrator.opts.d_discontinuities) &&
        first(integrator.opts.d_discontinuities) == integrator.tdir * integrator.t
        d2 = pop!(integrator.opts.d_discontinuities)
        order = min(order, d2.order)
    end

    # remove all discontinuities close to the current time point as well and
    # calculate minimal order of these discontinuities
    # integrator.EEst has unitless type of integrator.t
    if integrator.EEst isa AbstractFloat
        maxΔt = 10eps(integrator.t)

        while !isempty(integrator.opts.d_discontinuities) &&
            abs(first(integrator.opts.d_discontinuities).t - integrator.tdir * integrator.t) <
            maxΔt
            d2 = pop!(integrator.opts.d_discontinuities)
            order = min(order, d2.order)
        end

        # also remove all corresponding time stops
        while !isempty(integrator.opts.tstops) &&
            abs(first(integrator.opts.tstops) - integrator.tdir * integrator.t) < maxΔt
            pop!(integrator.opts.tstops)
        end
    end

    # add discontinuities of next order to integrator
    add_next_discontinuities!(integrator, order)
    nothing
end

"""
    add_next_discontinuities!(integrator::DDEIntegrator, order[, t=integrator.t])

Add discontinuities of next order that are propagated from discontinuity of
order `order` at time `t` in `integrator`, but only if `order` is less or equal
than the order of the applied method or the problem is neutral.

Discontinuities caused by constant delays are immediately calculated, and
discontinuities caused by dependent delays are tracked by a callback.
"""
function add_next_discontinuities!(integrator, _order, t = integrator.t)
    neutral = integrator.sol.prob.neutral
    order = Rational{Int}(_order)
    next_order = neutral ? order : order + 1 // 2

    # only track discontinuities up to order of the applied method
    alg_order = StochasticDiffEq.alg_order(getalg(integrator.alg))
    next_order <= alg_order + 1 // 2 || return

    # discontinuities caused by constant lags
    if has_constant_lags(integrator)
        constant_lags = integrator.sol.prob.constant_lags
        maxlag = integrator.tdir * (integrator.sol.prob.tspan[end] - t)

        for lag in constant_lags
            if integrator.tdir * lag < maxlag
                # calculate discontinuity and add it to heap of discontinuities and time stops
                d = Discontinuity(integrator.tdir * (t + lag), next_order)
                push!(integrator.opts.d_discontinuities, d)
                push!(integrator.opts.tstops, d.t)
            end
        end
    end

    # track propagated discontinuities with callback
    push!(integrator.tracked_discontinuities, Discontinuity(integrator.tdir * t, order))
    nothing
end
