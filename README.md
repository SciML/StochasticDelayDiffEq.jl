# StochasticDelayDiffEq.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/StochasticDelayDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/StochasticDelayDiffEq.jl)

StochasticDelayDiffEq.jl is a component package in the DifferentialEquatinos ecosystem.
It holds the stochastic delay differential equation solvers and utilities.
It is built on top of StochasticDiffEq to extend those solvers for stochastic delay differential equations. While completely independent and usable on its own, users interested in using this functionality should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) (documentation coming soon).

## API

StochasticDelayDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes a StochasticDiffEq.jl algorithm to `solve`.

Both constant and state-dependent lags are supported.
Interfacing with StochasticDiffEq.jl for implicit methods for stiff equations is not yet supported, but it is coming soon.

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl SDE Solvers page](http://docs.juliadiffeq.org/latest/solvers/sde_solve.html). For options for the `solve` command, see the [common solver options page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html).
