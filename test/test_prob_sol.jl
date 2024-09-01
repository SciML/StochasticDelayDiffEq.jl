using StochasticDelayDiffEq
using Test

#Hayes Equation
begin
    function hayes_modelf(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        du .= a .* u .+ b .* h(p, t - τ) .+ c
    end
    function hayes_modelg(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        du .= α .* u .+ β .* h(p, t - τ) .+ γ
    end
    h(p, t) = (ones(1) .+ t)
    tspan = (0.0, 0.1)
end

pmul = [1.0, -4.0, -2.0, 10.0, -1.3, -1.2, 1.1]
padd = [1.0, -4.0, -2.0, 10.0, -0.0, -0.0, 0.1]

prob = @test_nowarn SDDEProblem(hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
                                constant_lags = (pmul[1],));
@test_nowarn SDDEProblem(hayes_modelf, hayes_modelg, h, tspan, padd;
                         constant_lags = (padd[1],));

sol = @test_nowarn solve(prob, EM(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, LambaEM(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, EulerHeun(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, LambaEulerHeun(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, RKMil(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, RKMilCommute(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, RKMilCommute(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, WangLi3SMil_A(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, WangLi3SMil_B(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, WangLi3SMil_C(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, WangLi3SMil_D(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, WangLi3SMil_E(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, WangLi3SMil_F(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SRI(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SRIW1(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SRIW2(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SOSRI(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SOSRI2(), dt = 0.01)
@test sol.u[end] != zeros(1)

# Additive problems
println("additive problems")
prob.p .= padd;
sol = @test_nowarn solve(prob, SRA(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SRA1(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SRA2(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SRA3(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SOSRA(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SOSRA2(), dt = 0.01)
@test sol.u[end] != zeros(1)
# @test_nowarn solve(prob,SKenCarp(),dt=0.01) # Not working

# Test SROCK methods
println("SROCK methods")
prob.p .= pmul;
sol = @test_nowarn solve(prob, SROCK1(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SROCK1(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SROCKEM(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SROCKEM(strong_order_1 = false), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SROCK2(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SKSROCK(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SKSROCK(; post_processing = true), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, TangXiaoSROCK2(), dt = 0.01)
@test sol.u[end] != zeros(1)
for i in 1:5
    sol = @test_nowarn solve(prob, TangXiaoSROCK2(version_num = i), dt = 0.01)
    @test sol.u[end] != zeros(1)
end

# Test Implicit methods
println("implicit methods")
sol = @test_nowarn solve(prob, ImplicitEM(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ImplicitEM(symplectic = true, theta = 1 / 2), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ImplicitEulerHeun(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ImplicitEulerHeun(symplectic = true, theta = 1 / 2),
                         dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ImplicitRKMil(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ImplicitRKMil(symplectic = true, theta = 1 / 2), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob,
                         ImplicitRKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich, symplectic = true,
                                       theta = 1 / 2), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ISSEM(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ISSEM(symplectic = true, theta = 1 / 2), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ISSEulerHeun(), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, ISSEulerHeun(symplectic = true, theta = 1 / 2), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = @test_nowarn solve(prob, SKenCarp(), dt = 0.01)
@test sol.u[end] != zeros(1)
