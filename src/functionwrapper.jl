struct SDEDiffusionTermWrapper{iip, G, H}
    g::G
    h::H
end
function SDEDiffusionTermWrapper(g::G, h::H) where {G, H}
    SDEDiffusionTermWrapper{isinplace(g, 5), G, H}(g, h)
end
(g::SDEDiffusionTermWrapper{true})(du, u, p, t) = g.g(du, u, g.h, p, t)
(g::SDEDiffusionTermWrapper{false})(u, p, t) = g.g(u, g.h, p, t)

struct SDEFunctionWrapper{iip, F, G, H, TMM, Ta, Tt, TJ, JVP, VJP, JP, SP, TW, TWt, TPJ, GG,
                          TCV, ID, S} <: DiffEqBase.AbstractRODEFunction{iip}
    f::F
    g::G
    h::H
    mass_matrix::TMM
    analytic::Ta
    tgrad::Tt
    jac::TJ
    jvp::JVP
    vjp::VJP
    jac_prototype::JP
    sparsity::SP
    Wfact::TW
    Wfact_t::TWt
    paramjac::TPJ
    ggprime::GG
    colorvec::TCV
    initialization_data::ID
    sys::S
end

(f::SDEFunctionWrapper{true})(du, u, p, t) = f.f(du, u, f.h, p, t)
(f::SDEFunctionWrapper{false})(u, p, t) = f.f(u, f.h, p, t)

function wrap_functions_and_history(f::SDDEFunction, g, h)
    gwh = SDEDiffusionTermWrapper{isinplace(g, 5), typeof(g), typeof(h)}(g, h)

    if f.jac === nothing
        jac = nothing
    else
        if isinplace(f)
            jac = let f_jac = f.jac, h = h
                (J, u, p, t) -> f_jac(J, u, h, p, t)
            end
        else
            jac = let f_jac = f.jac, h = h
                (u, p, t) -> f_jac(u, h, p, t)
            end
        end
    end

    SDEFunctionWrapper{isinplace(f), typeof(f.f), typeof(gwh), typeof(h),
                       typeof(f.mass_matrix),
                       typeof(f.analytic), typeof(f.tgrad), typeof(jac), typeof(f.jvp),
                       typeof(f.vjp), typeof(f.jac_prototype), typeof(f.sparsity),
                       typeof(f.Wfact), typeof(f.Wfact_t), typeof(f.paramjac),
                       typeof(f.ggprime), typeof(f.colorvec), typeof(f.initialization_data),
                       typeof(f.sys)}(f.f, gwh, h, f.mass_matrix, f.analytic, f.tgrad, jac,
                       f.jvp, f.vjp, f.jac_prototype, f.sparsity, f.Wfact, f.Wfact_t, f.paramjac,
                       f.ggprime, f.colorvec, f.initialization_data, f.sys),
    gwh
end
