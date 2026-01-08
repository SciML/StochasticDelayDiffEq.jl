using PrecompileTools

@setup_workload begin
    # Simple SDDE problem for precompilation
    function _precompile_f!(du, u, h, p, t)
        τ = p[1]
        a = p[2]
        b = p[3]
        du[1] = a * u[1] + b * h(p, t - τ)[1]
        return nothing
    end
    function _precompile_g!(du, u, h, p, t)
        α = p[4]
        du[1] = α * u[1]
        return nothing
    end
    _precompile_h(p, t) = ones(1)
    _precompile_tspan = (0.0, 0.01)
    _precompile_p = [0.1, -1.0, 0.5, 0.1]
    _precompile_u0 = [1.0]

    @compile_workload begin
        # Precompile SDDEProblem construction
        prob = SDDEProblem(
            _precompile_f!,
            _precompile_g!,
            _precompile_u0,
            _precompile_h,
            _precompile_tspan,
            _precompile_p;
            constant_lags = (_precompile_p[1],)
        )

        # Precompile solve with EM (most common fixed-step algorithm)
        sol = solve(prob, EM(); dt = 0.001, save_everystep = false)
    end
end
