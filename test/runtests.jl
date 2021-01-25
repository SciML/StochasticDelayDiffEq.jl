using SafeTestsets

@safetestset "SDDEProblem, solve" begin include("test_prob_sol.jl") end
@safetestset "Analyticless Convergence Tests" begin include("analyticless_convergence_tests.jl") end
@safetestset "Event handling" begin include("events.jl") end
