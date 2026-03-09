function gmodel = model2gurobi(model)

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