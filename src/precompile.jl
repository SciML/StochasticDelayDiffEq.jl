using PrecompileTools

@setup_workload begin
    # Minimal imports needed for precompilation workload
    # Define a simple Hayes equation model (standard test case)
    function _precompile_f(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        du .= a .* u .+ b .* h(p, t - τ) .+ c
        return nothing
    end

    function _precompile_g(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        du .= α .* u .+ β .* h(p, t - τ) .+ γ
        return nothing
    end

    _precompile_h(p, t) = ones(1) .+ t

    _precompile_tspan = (0.0, 0.01)  # Very short timespan for fast precompilation
    _precompile_p = [0.01, -4.0, -2.0, 10.0, -1.3, -1.2, 1.1]  # Small lag for short integration
    _precompile_u0 = [1.0]

    @compile_workload begin
        # Precompile SDDEProblem creation - this is always the first step users take
        prob = SDDEProblem(
            _precompile_f,
            _precompile_g,
            _precompile_u0,
            _precompile_h,
            _precompile_tspan,
            _precompile_p;
            constant_lags = (_precompile_p[1],)
        )

        # Precompile the most common solver: EM (Euler-Maruyama)
        # This is the simplest and most commonly used method
        solve(prob, EM(); dt = 0.001)

        # Precompile SRIW1 - a common adaptive method
        solve(prob, SRIW1(); dt = 0.001)
    end
end
