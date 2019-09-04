abstract type AbstractSDDEFunction{iip} <: AbstractDiffEqFunction{iip} end

struct SDDEFunction{iip,F,G,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,GG,TCV} <: AbstractSDDEFunction{iip}
  f::F
  g::G
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  ggprime::GG
  syms::S
  colorvec::TCV
end

function SDDEFunction{iip,true}(f,g;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 syms = nothing,
                 colorvec = nothing)  where iip
                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 SDDEFunction{iip,typeof(f),typeof(g),
                 typeof(mass_matrix),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),
                 typeof(ggprime),typeof(_colorvec)}(
                 f,g,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,ggprime,syms,_colorvec)
end

function SDDEFunction{iip,false}(f,g;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 SDDEFunction{iip,Any,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),Any,typeof(_colorvec)}(
                 f,g,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,ggprime,syms,_colorvec)
end
SDDEFunction{iip}(f,g; kwargs...) where iip = SDDEFunction{iip,RECOMPILE_BY_DEFAULT}(f,g; kwargs...)
SDDEFunction{iip}(f::SDDEFunction,g; kwargs...) where iip = f
SDDEFunction(f,g; kwargs...) = SDDEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f,g; kwargs...)
SDDEFunction(f::SDDEFunction; kwargs...) = f


function Base.convert(::Type{SDDEFunction},f,g)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  SDDEFunction(f,g;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end

function Base.convert(::Type{SDDEFunction{iip}},f,g) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  SDDEFunction{iip,RECOMPILE_BY_DEFAULT}(f,g;analytic=analytic,
              tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end
