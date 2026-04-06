function model = Reduce_S_model(model, idx_act_rxn)
%--------------------------------------------------------------------------
% Reduces the stoichiometric matrix of the metabolic model by removing rows that correspond to
% metabolites that are not involved in the active reactions. This is done to reduce the
% dimensionality of the optimization problem.
%
% Inputs:
%   model       - COBRA metabolic model structure.
%   idx_act_rxn - Indices of active reactions in the model.
%
% Outputs:
%   model - Reduced metabolic model structure with only active reactions and involved metabolites.
%--------------------------------------------------------------------------

non_zero_row   = find(~all(model.S(:, idx_act_rxn) == 0, 2));

model.S    = model.S(non_zero_row, idx_act_rxn);
model.b    = model.b(non_zero_row);
model.lb   = model.lb(idx_act_rxn);
model.ub   = model.ub(idx_act_rxn);
model.rxns = model.rxns(idx_act_rxn);
model.c    = model.c(idx_act_rxn);
model.mets = model.mets(non_zero_row);

end
