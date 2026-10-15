function model = make_aa_irrev(model, index_PAN)

% New components
ub_mod = model.ub;
ub_mod(index_PAN) = ub_mod(index_PAN)*0;

lb_mod = model.lb;
lb_mod(index_PAN) = lb_mod(index_PAN)*0;

%Rxn names
rxn_code  = model.rxns(index_PAN);
rxn_names = model.rxnNames(index_PAN);

% Modify model

model.S  = [model.S, model.S(:,index_PAN)*-1];
model.lb = [lb_mod,; ub_mod(index_PAN)];
model.ub = [model.ub; model.ub(index_PAN)];
model.rxns = [model.rxns; rxn_code];
model.rxnNames = [model.rxnNames; rxn_names];
model.c        = [model.c; zeros(size(index_PAN))];


end