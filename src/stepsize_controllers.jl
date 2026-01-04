function stepsize_controller!(integrator::SDDEIntegrator, controller::PIController, alg)
    integrator.q11 = DiffEqBase.value(FastPower.fastpower(integrator.EEst, controller.beta1))
    integrator.q = DiffEqBase.value(
        integrator.q11 /
            FastPower.fastpower(integrator.qold, controller.beta2)
    )
    return @fastmath integrator.q = DiffEqBase.value(
        max(
            inv(integrator.opts.qmax),
            min(
                inv(integrator.opts.qmin),
                integrator.q / integrator.opts.gamma
            )
        )
    )
end

@inline function step_accept_controller!(integrator::SDDEIntegrator, alg)
    return step_accept_controller!(integrator, integrator.opts.controller, alg)
end

function step_accept_controller!(integrator::SDDEIntegrator, controller::PIController, alg)
    return integrator.dtnew = DiffEqBase.value(integrator.dt / integrator.q) *
        oneunit(integrator.dt)
end

function step_reject_controller!(integrator::SDDEIntegrator, controller::PIController, alg)
    return integrator.dtnew = integrator.dt / min(
        inv(integrator.opts.qmin),
        integrator.q11 / integrator.opts.gamma
    )
end
