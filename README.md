# StochasticDelayDiffEq.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/StochasticDelayDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/StochasticDelayDiffEq.jl)

StochasticDelayDiffEq.jl is a component package in the DifferentialEquatinos ecosystem.
It holds the delay differential equation solvers and utilities.
s built on top of StochasticDiffEq to extend those solvers for delay differential equations. While completely independent and usable on its own, users interested in using this functionality should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## API

StochasticDelayDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes a DelayDiffEq.jl algorithm to `solve`. For example, we can solve the [DDE tutorial from the documentation](http://docs.juliadiffeq.org/latest/tutorials/dde_example.html) using the `MethodOfSteps(Tsit5())` algorithm:


Both constant and state-dependent lags are supported. Interfacing with OrdinaryDiffEq.jl for implicit methods for stiff equations is also supported.

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl DDE Solvers page](http://docs.juliadiffeq.org/latest/solvers/dde_solve.html). For options for the `solve` command, see the [common solver options page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html).
