# StochasticDelayDiffEq.jl

[![Build Status](https://github.com/SciML/StochasticDelayDiffEq.jl/workflows/CI/badge.svg)](https://github.com/SciML/StochasticDelayDiffEq.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/SciML/StochasticDelayDiffEq.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/StochasticDelayDiffEq.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/StochasticDelayDiffEq.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/StochasticDelayDiffEq.jl?branch=master)

StochasticDelayDiffEq.jl is a component package in the DifferentialEquations ecosystem.
It holds the stochastic delay differential equation solvers and utilities.
It is built on top of StochasticDiffEq to extend those solvers for stochastic delay differential equations. While completely independent and usable on its own, users interested in using this functionality should check out [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) (documentation coming soon).

## API

StochasticDelayDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes a StochasticDiffEq.jl algorithm to `solve`.

Both constant and state-dependent lags are supported.
Interfacing with StochasticDiffEq.jl for implicit methods for stiff equations is not yet supported, but it is coming soon.

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl SDE Solvers page](https://diffeq.sciml.ai/stable/solvers/sde_solve/). For options for the `solve` command, see the [common solver options page](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
