function gmodel = model2gurobi(model)
%--------------------------------------------------------------------------
% Converts a metabolic model in COBRA format to a Gurobi optimization model.
% Inputs:
%   model - COBRA metabolic model structure.
%
% Outputs:
%   gmodel - Gurobi optimization model structure.
%--------------------------------------------------------------------------
% Gurobi params
gmodel.A     = sparse(model.S);
gmodel.lb    = model.lb;
gmodel.ub    = model.ub;
gmodel.obj   = model.c;
gmodel.rhs   = model.b;
gmodel.sense = [repmat('=', size(model.S,1), 1)];
gmodel.modelsense = 'Max';

% Extra attributes of metabolic model
gmodel.rxns = model.rxns;
gmodel.rxnNames = model.rxnNames;
gmodel.mets = model.mets;
gmodel.metNames = model.metNames;


end