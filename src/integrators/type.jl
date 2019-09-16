mutable struct HistorySDEIntegrator{algType,IIP,uType,tType,SolType,CacheType} <: AbstractSDEIntegrator{algType,IIP,uType,tType}
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

mutable struct
    SDDEIntegrator{algType,IIP,uType,uEltype,tType,P,eigenType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,cacheType,F,G,OType,noiseType,EventErrorType,CallbackCacheType,H,IType} <: AbstractSDDEIntegrator{algType,IIP,uType,tType}
    f::F
    g::G
    noise::noiseType
    uprev::uType
    tprev::tType
    # prev_idx::Int - from DDEIntegrator TODO
    # prev2_idx::Int - from DDEIntegrator TODO
    # fpsolver::FP - from DDEIntegrator TODO
    order_discontinuity_t0::Int  #- from DDEIntegrator TODO
    tracked_discontinuities::Vector{Discontinuity{tType}}
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
    opts::OType
    iter::Int
    success_iter::Int
    eigen_est::eigenType
    EEst::tTypeNoUnits
    q::tTypeNoUnits
    qold::tTypeNoUnits
    q11::tTypeNoUnits
    history::H
    integrator::IType # history integrator
  end
  
  function (integrator::SDDEIntegrator)(t, deriv::Type=Val{0}; idxs=nothing)
      StochasticDiffEq.current_interpolant(t, integrator, idxs, deriv)
  end
  
  function (integrator::SDDEIntegrator)(val::AbstractArray, t::Union{Number,AbstractArray},deriv::Type=Val{0}; idxs=nothing)
      StochasticDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
  end