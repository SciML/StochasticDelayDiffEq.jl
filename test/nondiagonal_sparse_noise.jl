using StochasticDelayDiffEq
using Random
using SparseArrays

function sir_dde!(du, u, h, p, t)
    (S, I, R) = u
    (β, c, τ) = p
    N = S + I + R
    infection = β * c * I / N * S
    (Sd, Id, Rd) = h(p, t - τ) # Time delayed variables
    Nd = Sd + Id + Rd
    recovery = β * c * Id / Nd * Sd
    @inbounds begin
        du[1] = -infection
        du[2] = infection - recovery
        du[3] = recovery
    end
    return nothing
end;

# Define a sparse matrix by making a dense matrix and setting some values as not zero
A = zeros(3, 2)
A[1, 1] = 1
A[2, 1] = 1
A[2, 2] = 1
A[3, 2] = 1
A = SparseArrays.sparse(A);

# Make `g` write the sparse matrix values
# Use max(0, x) to prevent DomainError from sqrt when stochastic dynamics cause negative values
function sir_delayed_noise!(du, u, h, p, t)
    (S, I, R) = u
    (β, c, τ) = p
    N = S + I + R
    infection = β * c * I / N * S
    (Sd, Id, Rd) = h(p, t - τ) # Time delayed variables
    Nd = Sd + Id + Rd
    recovery = β * c * Id / Nd * Sd
    du[1, 1] = -sqrt(max(0.0, infection))
    du[2, 1] = sqrt(max(0.0, infection))
    du[2, 2] = -sqrt(max(0.0, recovery))
    return du[3, 2] = sqrt(max(0.0, recovery))
end;

function condition(u, t, integrator) # Event when event_f(u,t) == 0
    return u[2]
end;
function affect!(integrator)
    return integrator.u[2] = 0.0
end;
cb = ContinuousCallback(condition, affect!);

δt = 0.1
tmax = 40.0
tspan = (0.0, tmax)
t = 0.0:δt:tmax;
u0 = [990.0, 10.0, 0.0]; # S,I,R

function sir_history(p, t)
    return [1000.0, 0.0, 0.0]
end;

p = [0.05, 10.0, 4.0]; # β,c,τ
Random.seed!(1234);

prob_sdde = SDDEProblem(
    sir_dde!, sir_delayed_noise!, u0, sir_history, tspan, p;
    noise_rate_prototype = A
);
sol_sdde = solve(prob_sdde, LambaEM(), callback = cb);
