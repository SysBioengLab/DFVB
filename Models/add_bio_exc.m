%% e coli core model
model = load("e_coli_core.mat").e_coli_core;
e_coli_core = addMetabolite(model,'biomass_e');

idx_growth_rxn = model.c==1;
S_mod = [model.S, zeros(size(model.S,1),1);
         zeros(1, size(model.S,2)+1)];
S_mod(end, idx_growth_rxn) = 1;
S_mod(end, end) = -1;

e_coli_core.S               = S_mod;
e_coli_core.rxns{end+1}     = 'EX_bio_e';
e_coli_core.rxnNames{end+1} = 'Biomass Exchange';
e_coli_core.lb(end+1)       = -1000;
e_coli_core.ub(end+1)       =  1000;
e_coli_core.c               = [zeros(size(e_coli_core.c,1),1); 1];

save('Models/e_coli_core_mod.mat', 'e_coli_core')

%% iML1515 model
model   = load("iML1515.mat").iML1515;
IML1515 = addMetabolite(model,'biomass_e');

idx_growth_rxn = model.c==1;
S_mod = [model.S, zeros(size(model.S,1),1);
         zeros(1, size(model.S,2)+1)];
S_mod(end, idx_growth_rxn) = 1;
S_mod(end, end) = -1;

IML1515.S               = S_mod;
IML1515.rxns{end+1}     = 'EX_bio_e';
IML1515.rxnNames{end+1} = 'Biomass Exchange';
IML1515.lb(end+1)       = -1000;
IML1515.ub(end+1)       =  1000;
IML1515.c               = [zeros(size(IML1515.c,1),1); 1];



save('Models/iML1515_mod.mat', 'IML1515')

