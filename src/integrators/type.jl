mutable struct HistorySDEIntegrator{algType, IIP, uType, tType, SolType, CacheType} <:
               AbstractSDEIntegrator{algType, IIP, uType, tType}
    sol::SolType
    u::uType
    t::tType
    dt::tType
    uprev::uType
    tprev::tType
    alg::algType
    dtcache::tType
    tdir::Int
    saveiter::Int
    cache::CacheType
end

function (integrator::HistorySDEIntegrator)(t, deriv::Type = Val{0}; idxs = nothing)
    StochasticDiffEq.current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::HistorySDEIntegrator)(val::AbstractArray,
                                            t::Union{Number, AbstractArray},
                                            deriv::Type = Val{0}; idxs = nothing)
    StochasticDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
end

mutable struct SDDEIntegrator{algType, IIP, uType, uEltype, tType, P, eigenType,
                              tTypeNoUnits, uEltypeNoUnits, randType, randType2, rateType,
                              solType, cacheType, F, G, F6, OType, noiseType,
                              EventErrorType, CallbackCacheType, H, IType, IA} <:
               AbstractSDDEIntegrator{algType, IIP, uType, tType}
    f::F
    g::G
    c::F6
    noise::noiseType
    uprev::uType
    tprev::tType
    # prev_idx::Int - from DDEIntegrator TODO
    # prev2_idx::Int - from DDEIntegrator TODO
    # fpsolver::FP - from DDEIntegrator TODO
    order_discontinuity_t0::Int  #- from DDEIntegrator TODO
    tracked_discontinuities::Vector{Discontinuity{tType, Rational{Int}}}
    # discontinuity_interp_points::Int #- from DDEIntegrator TODO
    # discontinuity_abstol::dAbsType #- from DDEIntegrator TODO
    # discontinuity_reltol::dRelType #- from DDEIntegrator TODO
    t::tType
    u::uType
    p::P
    dt::tType
    dtnew::tType
    dtpropose::tType
    dtcache::tType
    T::tType
    tdir::Int
    just_hit_tstop::Bool
    do_error_check::Bool
    isout::Bool
    event_last_time::Int
    vector_event_last_time::Int
    last_event_error::EventErrorType
    accept_step::Bool
    last_stepfail::Bool
    force_stepfail::Bool
    dtchangeable::Bool
    u_modified::Bool
    saveiter::Int
    alg::algType
    sol::solType
    cache::cacheType
    callback_cache::CallbackCacheType
    sqdt::tType
    W::randType
    P::randType2
    opts::OType
    iter::Int
    success_iter::Int
    eigen_est::eigenType
    EEst::tTypeNoUnits
    q::tTypeNoUnits
    qold::tTypeNoUnits
    q11::tTypeNoUnits
    history::H
    stats::DiffEqBase.Stats
    integrator::IType # history integrator
    initializealg::IA
end

function (integrator::SDDEIntegrator)(t, deriv::Type = Val{0}; idxs = nothing)
    StochasticDiffEq.current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::SDDEIntegrator)(val::AbstractArray, t::Union{Number, AbstractArray},
                                      deriv::Type = Val{0}; idxs = nothing)
    StochasticDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
end
