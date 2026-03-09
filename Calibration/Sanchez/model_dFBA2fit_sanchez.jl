using JuMP
using Gurobi

# Estructura para mantener el estado de warm start
def mutable struct WarmStartState
    vbasis::Union{Nothing, Vector{Int}}
    cbasis::Union{Nothing, Vector{Int}}
    status::Union{Nothing, String}
    sol::Union{Nothing, Vector{Float64}}
    WarmStartState() = new(nothing, nothing, nothing, nothing)
end

global warm_state = WarmStartState()

function model_dFBA2fit_sanchez!(dxdt, x, p, t)
    # x: vector de estados [X, G, E]
    # p: (theta, model, index)
    theta, model, index = p

    X = x[1] # Biomasa
    G = x[2] # Glucosa
    E = x[3] # Etanol

    index_mu, index_glc, index_o2, index_etoh = index

    MW_glc  = 180.156/1e3
    MW_etoh = 46.07/1e3
    MW_lac  = 90.08/1e3
    MW_gly  = 92.09382/1e3
    MW_cit  = 192.124/1e3

    # Parámetros estimados
    v_glc, K_G, K_iE, lb_o2, f_gly, f_cit, f_lac = theta

    # Limitar consumo si se agota
    rS = G > 1e-3 ? v_glc*G/((K_G+G)*(1+E/K_iE)) : 0.0

    # Actualizar restricciones en el modelo
    set_lower_bound(model, index_glc, -rS)
    set_lower_bound(model, index_o2, -lb_o2)

    # Warm start
    if warm_state.vbasis !== nothing && warm_state.cbasis !== nothing && warm_state.status == "OPTIMAL"
        MOI.set(model, Gurobi.VariableBasis(), warm_state.vbasis)
        MOI.set(model, Gurobi.ConstraintBasis(), warm_state.cbasis)
    end

    optimize!(model)
    status = termination_status(model)

    if status == MOI.OPTIMAL
        mu   = value(model[:x][index_mu])
        rS_  = value(model[:x][index_glc])
        rE   = value(model[:x][index_etoh])
        # Guardar basis para warm start
        warm_state.vbasis = MOI.get(model, Gurobi.VariableBasis())
        warm_state.cbasis = MOI.get(model, Gurobi.ConstraintBasis())
        warm_state.status = "OPTIMAL"
        # Ecuaciones diferenciales
        dxdt[1] = mu * X
        dxdt[2] = rS_ * X * MW_glc
        dxdt[3] = rE * X * MW_etoh
        dxdt[4] = f_gly * mu * X * MW_gly
        dxdt[5] = f_cit * mu * X * MW_cit
        dxdt[6] = f_lac * mu * X * MW_lac
    else
        warm_state.status = "INFEASIBLE"
        dxdt .= 0.0
    end
end
