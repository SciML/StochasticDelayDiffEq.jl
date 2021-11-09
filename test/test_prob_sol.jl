using StochasticDelayDiffEq
using Test

#Hayes Equation
begin
    function hayes_modelf(du,u,h,p,t)
        τ,a,b,c,α,β,γ = p
        du .= a.*u .+ b .* h(p,t-τ) .+ c
    end
    function hayes_modelg(du,u,h,p,t)
        τ,a,b,c,α,β,γ = p
        du .= α.*u .+ β.*h(p,t-τ) .+ γ
    end
    h(p,t) = (ones(1) .+ t);
    tspan = (0.,10.)
end


pmul = [1.0,-4.,-2.,10.,-1.3,-1.2, 1.1]
padd = [1.0,-4.,-2.,10.,-0.0,-0.0, 0.1]

prob = SDDEProblem(hayes_modelf, hayes_modelg, [1.], h, tspan, pmul; constant_lags = (pmul[1],));
SDDEProblem(hayes_modelf, hayes_modelg, h, tspan, padd; constant_lags = (padd[1],));

sol = solve(prob,EM(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,LambaEM(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,EulerHeun(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,LambaEulerHeun(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,RKMil(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,RKMil(interpretation=:Stratonovich),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,RKMilCommute(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,RKMilCommute(interpretation=:Stratonovich),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,WangLi3SMil_A(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,WangLi3SMil_B(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,WangLi3SMil_C(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,WangLi3SMil_D(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,WangLi3SMil_E(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,WangLi3SMil_F(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SRI(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SRIW1(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SRIW2(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SOSRI(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SOSRI2(),dt=0.01)
@test sol.u[end] != zeros(1)

# Additive problems
println("additive problems")
prob.p .= padd;
sol = solve(prob,SRA(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SRA1(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SRA2(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SRA3(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SOSRA(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SOSRA2(),dt=0.01)
@test sol.u[end] != zeros(1)
# solve(prob,SKenCarp(),dt=0.01) # Not working

# Test SROCK methods
println("SROCK methods")
prob.p .= pmul;
sol = solve(prob,SROCK1(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SROCK1(interpretation=:Stratonovich),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SROCKEM(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SROCKEM(strong_order_1=false),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SROCK2(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SKSROCK(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SKSROCK(;post_processing=true),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,TangXiaoSROCK2(),dt=0.01)
@test sol.u[end] != zeros(1)
for i in 1:5
    sol = solve(prob,TangXiaoSROCK2(version_num=i),dt=0.01)
    @test sol.u[end] != zeros(1)
end

# Test Implicit methods
println("implicit methods")
sol = solve(prob,ImplicitEM(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ImplicitEM(symplectic=true, theta = 1/2),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ImplicitEulerHeun(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ImplicitEulerHeun(symplectic=true, theta = 1/2),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ImplicitRKMil(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ImplicitRKMil(symplectic=true, theta = 1/2),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ImplicitRKMil(interpretation=:Stratonovich, symplectic = true, theta = 1/2),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ISSEM(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ISSEM(symplectic=true, theta = 1/2),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ISSEulerHeun(),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,ISSEulerHeun(symplectic=true, theta = 1/2),dt=0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob,SKenCarp(),dt=0.01)
@test sol.u[end] != zeros(1)
