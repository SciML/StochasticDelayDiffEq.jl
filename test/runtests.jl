import Pkg; Pkg.activate("."); #Pkg.instantiate()

using Revise
using StochasticDelayDiffEq
using Plots

begin
    const out = zeros(1);
    function hayes_modelf(du,u,h,p,t)
        τ,a,b,c,α,β,γ = p
        du .= a.*u .+ b .* h(p,t-τ) .+ c
    end
    function hayes_modelg(du,u,h,p,t)
        τ,a,b,c,α,β,γ = p
        du .= α.*u .+ β.*h(p,t-τ) .+ γ
    end
    h(p,t) = (ones(1) .+ t);
    tspan = (0.,100.)
end

p = [1.0,-4.,-2.,10.,-0.3,-0.2, 0.1]
# p = [1.0,-3.,-2.,0.,-0.0,-0.0, 0.2]

prob = SDDEProblem(hayes_modelf, hayes_modelg, h, tspan, p; constant_lags = (p[1],));

solEM=StochasticDelayDiffEq.solve(prob,StochasticDelayDiffEq.EM(),dt=0.001)
solLEM=StochasticDelayDiffEq.solve(prob,StochasticDelayDiffEq.LambaEM(),dt=0.001)

plot(sol)