using StochasticDelayDiffEq, DiffEqDevTools
using Test

#Hayes Equation
begin
    function hayes_modelf(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        du .= a .* u .+ b .* h(p, t - τ) .+ c
    end
    function hayes_modelg(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        du .= α .* u .+ γ
    end
    h(p, t) = (ones(1) .+ t)
    tspan = (0.0, 10.0)
end

pmul = [1.0, -4.0, -2.0, 10.0, -1.3, -1.2, 1.1]
padd = [1.0, -4.0, -2.0, 10.0, -0.0, -0.0, 0.1]

prob = @test_nowarn SDDEProblem(hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
                                constant_lags = (pmul[1],));
@test_nowarn SDDEProblem(hayes_modelf, hayes_modelg, h, tspan, padd;
                         constant_lags = (padd[1],));

dts = (1 / 2) .^ (8:-1:4)
test_dt = 1 / 2^9
sim2 = analyticless_test_convergence(dts, prob, EM(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob, LambaEM(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob, EulerHeun(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, LambaEulerHeun(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, RKMil(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
                                     test_dt, trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, WangLi3SMil_A(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, WangLi3SMil_B(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, WangLi3SMil_C(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, WangLi3SMil_D(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, WangLi3SMil_E(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, WangLi3SMil_F(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3

# Test SROCK methods
println("SROCK methods")
prob.p .= pmul;
sim2 = analyticless_test_convergence(dts, prob, SROCK1(), test_dt, trajectories = 100,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, SROCK1(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
                                     test_dt, trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, SROCKEM(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, SROCKEM(strong_order_1 = false), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob, SKSROCK(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob, SKSROCK(; post_processing = true), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3

# Test Implicit methods
println("implicit methods")
sim2 = analyticless_test_convergence(dts, prob, ImplicitEM(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob,
                                     ImplicitEM(symplectic = true, theta = 1 / 2), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob, ImplicitEulerHeun(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob,
                                     ImplicitEulerHeun(symplectic = true, theta = 1 / 2),
                                     test_dt, trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, ISSEM(), test_dt, trajectories = 1000,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.35
sim2 = analyticless_test_convergence(dts, prob, ISSEM(symplectic = true, theta = 1 / 2),
                                     test_dt, trajectories = 500, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 0.5) < 0.3
sim2 = analyticless_test_convergence(dts, prob, ImplicitRKMil(), test_dt,
                                     trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob,
                                     ImplicitRKMil(symplectic = true, theta = 1 / 2),
                                     test_dt, trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob,
                                     ImplicitRKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich,
                                                   symplectic = true, theta = 1 / 2),
                                     test_dt, trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob, ISSEulerHeun(), test_dt, trajectories = 300,
                                     use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
sim2 = analyticless_test_convergence(dts, prob,
                                     ISSEulerHeun(symplectic = true, theta = 1 / 2),
                                     test_dt, trajectories = 300, use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
