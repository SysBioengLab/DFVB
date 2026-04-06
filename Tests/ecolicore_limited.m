%------------------------------------------------------%
%------------ File to test DFVB ecolicore -------------%
%------------------------------------------------------%

changeCobraSolverParams('LP','optTol',1e-8);
changeCobraSolverParams('LP','feasTol',1e-8);

%% Load the model

% initCobraToolbox(false)
% clear, clc
model    = load('e_coli_core_mod.mat').e_coli_core;
N_sparse = sparseNull(model.S); 

%% Blocking non-essential reactions

exchanges  = model.rxns(findExcRxns(model));

FuentesC   = {'EX_acald_e', 'EX_akg_e', 'EX_etoh_e', 'EX_for_e', 'EX_fru_e',...
              'EX_fum_e', 'EX_lac__D_e', 'EX_mal__L_e','EX_pyr_e', 'EX_succ_e','EX_co2_e'};

acetato =  'EX_ac_e';

esenciales = {'EX_co2_e', 'EX_h_e', 'EX_h2o_e', 'EX_o2_e', 'EX_pi_e'};

FuentesN ={'EX_gln__L_e','EX_glu__L_e'};

%% Lower and upper bounds 
 
model = changeRxnBounds(model,'EX_o2_e',-20,'l');
model = changeRxnBounds(model,'EX_ac_e',-1000,'l');
model = changeRxnBounds(model, 'EX_glc__D_e', -1000, 'l');

%% Calculate the model with conversions for DFVB and DFBA with gurobi

modelo_DFVB_N_sparse   = build_DFVB(model,N_sparse);
modelo_DFBA            = model2gurobi(model);

%% Reaction indices

bio_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_bio_e'));
bio_met_DFBA      = findMetIDs(model, 'biomass_e');
bio_flux_DFBA     = findRxnIDs(model, 'EX_bio_e');

glc_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_glc__D_e'));
glc_met_DFBA      = findMetIDs(model, 'glc__D_e');
glc_flux_DFBA     = findRxnIDs(model, 'EX_glc__D_e');
glc_ex_fba        = find(strcmp(model.rxns(findExcRxns(model)), 'EX_glc__D_e'));

ac_met_DFVB_N_sparse     = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_ac_e'));
ac_met_DFBA       = findMetIDs(model, 'ac_e');
ac_flux_DFBA      = findRxnIDs(model, 'EX_ac_e');
ac_ex_fba         = find(strcmp(model.rxns(findExcRxns(model)), 'EX_ac_e'));

o2_met_DFVB_N_sparse     = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_o2_e'));
o2_met_DFBA       = find(startsWith(model.mets, 'o2_e'));
o2_flux_DFBA      = find(startsWith(model.rxns, 'EX_o2_e'));
o2_ex_fba        = find(strcmp(model.rxns(findExcRxns(model)), 'EX_o2_e'));

co2_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_co2_e'));
co2_met_DFBA      = findMetIDs(model, 'co2_e');
co2_flux_DFBA     = findRxnIDs(model, 'EX_co2_e');
co2_ex_fba        = find(strcmp(model.rxns(findExcRxns(model)), 'EX_co2_e'));

nh4_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_nh4_e'));
nh4_met_DFBA      = findMetIDs(model, 'nh4_e');
nh4_flux_DFBA     = findRxnIDs(model, 'EX_nh4_e');
nh4_ex_fba        = find(strcmp(model.rxns(findExcRxns(model)), 'EX_nh4_e'));

etoh_met_DFVB_N_sparse       = find(startsWith(modelo_DFVB_N_sparse.rxns, 'EX_etoh_e'));
etoh_met_DFBA     = findMetIDs(model, 'etoh_e');
etoh_flux_DFBA    = findRxnIDs(model, 'EX_etoh_e');
etoh_ex_fba       = find(strcmp(model.rxns(findExcRxns(model)), 'EX_etoh_e'));

tags = {'Time', 'Biomass', 'Glucose', 'Oxygen', 'Acetate', 'Ammonium', 'Ethanol', 'CO2'};
%% Execution of DFVB Real conversions

% ODE configuration
t         = [0, 3];                  % Time interval
y0        = zeros(7,1);              % Initial conditions
options   = odeset('AbsTol',1e-9, 'RelTol', 1e-8, 'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 0.5;                    % Initial biomass condition
y0(2)      = 10;                     % Initial glucose condition
y0(3)      = 0.3125;                 % Initial oxygen condition
y0(4)      = 0;                      % Initial acetate condition
y0(5)      = 27.7;                   % Initial ammonium condition
y0(6)      = 0;                      % Initial ethanol condition
y0(7)      = 0;                      % Initial CO2 condition
params     = [10 0.015 20 0.003 50];                  % Glucose consumption parameters (V_max, K_S)

index_N_sparse   = [bio_met_DFVB_N_sparse, glc_met_DFVB_N_sparse, o2_met_DFVB_N_sparse, ac_met_DFVB_N_sparse, nh4_met_DFVB_N_sparse, etoh_met_DFVB_N_sparse, co2_met_DFVB_N_sparse];
alpha_dim_N_sparse = size(modelo_DFVB_N_sparse.nullS,2);

% nullSparse execution
tic; %Start the timer
[TN_sparse,YN_sparse] = ode113(@DFVB_o2_limited,t,y0,options,modelo_DFVB_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);
tiempo_transcurrido_N_sparse = toc; %Save total execution time
disp(['DFVB execution time: ' num2str(tiempo_transcurrido_N_sparse) ' seconds']);
Results_N_sparse = table(TN_sparse, YN_sparse(:,1), YN_sparse(:,2), YN_sparse(:,3), YN_sparse(:,4), YN_sparse(:,5), YN_sparse(:,6), YN_sparse(:,7), VariableNames=tags);

%% Average DFVB time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
%sparseNull
tic;
[~,~] = ode113(@DFVB_o2_limited,t,y0,options,modelo_DFVB_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);
t_N_sparse(i) = toc;
end

fprintf('\nMean CPU time of DFVB = %.4f s\n',mean(t_N_sparse))
fprintf('Std CPU time of DFVB = %.4f s\n',std(t_N_sparse)*2)

%% Recover non-forced fluxes and alphas

[alpha_N_sparse  , fluxes_N_sparse]   = Variable_recovery(TN_sparse,YN_sparse,modelo_DFVB_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);
Results_fluxes_N_sparse = table(TN_sparse, fluxes_N_sparse(:,bio_met_DFVB_N_sparse), fluxes_N_sparse(:,glc_met_DFVB_N_sparse), fluxes_N_sparse(:,o2_met_DFVB_N_sparse), fluxes_N_sparse(:,ac_met_DFVB_N_sparse), fluxes_N_sparse(:,nh4_met_DFVB_N_sparse), fluxes_N_sparse(:,etoh_met_DFVB_N_sparse), fluxes_N_sparse(:,co2_met_DFVB_N_sparse), VariableNames=tags);

% sparseNull
alpha_no_nulo_N_sparse = [];
for i=1:size(alpha_N_sparse,2)
    if any(alpha_N_sparse(:,i))
        alpha_no_nulo_N_sparse = [alpha_no_nulo_N_sparse i];   
    end
end
inactivos_N_sparse = setdiff(1:size(alpha_N_sparse,2), alpha_no_nulo_N_sparse);

fprintf('\nTotal conversions of sparseNull = %d', size(alpha_N_sparse,2))
fprintf('\nConversions used with sparseNull = %d\n', size(alpha_N_sparse,2) - length(inactivos_N_sparse))


%% DFBA execution

% ODE configuration
t         = [0, 3];                  % Time interval
y0        = ones(73,1)*0;            % Initial conditions
options   = odeset('AbsTol',1e-9, 'RelTol', 1e-8, 'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(bio_met_DFBA) = 0.5;                 % Initial glucose condition
y0(glc_met_DFBA) = 10;                  % Initial biomass condition
y0(o2_met_DFBA)  = 0.3125;              % Initial oxygen condition
y0(ac_met_DFBA)  = 0;                   % Initial acetate condition
y0(co2_met_DFBA) = 0;                   % Initial co2 condition
y0(nh4_met_DFBA) = 27.7;                % Initial nh4 condition
params     = [10 0.015 20 0.003 50];                  % Glucose consumption parameters (V_max, K_S)
idx_met    = [bio_met_DFBA, glc_met_DFBA, o2_met_DFBA, ac_met_DFBA, co2_met_DFBA, nh4_met_DFBA, etoh_met_DFBA];
idx_flux   = [bio_flux_DFBA, glc_flux_DFBA, o2_flux_DFBA, ac_flux_DFBA, co2_flux_DFBA, nh4_flux_DFBA, etoh_flux_DFBA];

tic; %Start the timer
[T2,Y2] = ode113(@DFBA_o2_limited_gurobi,t,y0,options,modelo_DFBA,idx_met,idx_flux,params);

tiempo_transcurrido = toc; %Save total execution time
disp(['DFBA execution time: ' num2str(tiempo_transcurrido) ' seconds']);

Results_DFBA = table(T2, Y2(:,bio_met_DFBA), Y2(:,glc_met_DFBA), Y2(:,o2_met_DFBA), Y2(:,ac_met_DFBA), Y2(:,nh4_met_DFBA), Y2(:,etoh_met_DFBA), Y2(:,co2_met_DFBA), VariableNames=tags);


%% Average DFBA time calculation with gurobi
n = 5;
t_DFBA = zeros(n,1);

for i = 1:n
tic;
[~,~] = ode113(@DFBA_o2_limited_gurobi,t,y0,options,modelo_DFBA,idx_met,idx_flux,params);
t_DFBA(i) = toc;
end
fprintf('\nMean CPU time of DFBA = %.4f s\n',mean(t_DFBA))
fprintf('Std CPU time of DFBA = %.4f s\n',std(t_DFBA))

%% Recuperar fluxes FBA

[fluxes_fba, fva_min, fva_max] = Variable_recovery_FBA(T2,Y2,model,idx_met,idx_flux,params, true);
fva_range.fva_min = fva_min;
fva_range.fva_max = fva_max;
Results_flux_DFBA = table(T2, fluxes_fba(bio_flux_DFBA,:)', fluxes_fba(glc_flux_DFBA,:)', fluxes_fba(o2_flux_DFBA,:)', fluxes_fba(ac_flux_DFBA,:)', fluxes_fba(nh4_flux_DFBA,:)', fluxes_fba(etoh_flux_DFBA,:)', fluxes_fba(co2_flux_DFBA,:)', VariableNames=tags);


%% Comparison plot

% sparseNull
figure();
subplot(1,2,1);
plot(TN_sparse,YN_sparse(:,1),'m','LineWidth',1.5)
hold on
plot(T2,Y2(:,bio_met_DFBA),'b--','LineWidth',1.5)
plot(TN_sparse,YN_sparse(:,2),'r','LineWidth',1.5)
plot(T2,Y2(:,glc_met_DFBA),'g--','LineWidth',1.5)
grid on
legend('Biomass DFVB','Biomass DFBA','Glucose DFVB','Glucose DFVB','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Biomass and Glucose')  

subplot(1,2,2)

yyaxis left
plot(TN_sparse,YN_sparse(:,3),'m',T2,Y2(:,o2_met_DFBA),'b--','LineWidth',1.5)


yyaxis right
hold on
plot(T2, Y2(:,ac_met_DFBA), 'g', 'LineWidth', 1.5)
plot(TN_sparse, YN_sparse(:,4), 'r--', 'LineWidth', 1.5)
xlabel('Time')
legend('DFVB', 'DFBA')
grid on
legend('Oxygen DFVB','Oxygen DFBA','Acetate DFVB','Acetate DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Oxygen and Acetate') 

%% Identify inactive

inactivos_N_sparse = setdiff(1:size(alpha_N_sparse,2), alpha_no_nulo_N_sparse);
activos_N_sparse   = setdiff(1:size(modelo_DFVB_N_sparse.nullS,2),inactivos_N_sparse);

%% Identify inactive with FVA

[min_f, max_f] = fluxVariability(model,100);
fva_inactivadas = find(min_f == 0 & max_f == 0);
exch_inactivas_fva = intersect(fva_inactivadas, findRxnIDs(model,exchanges));
int_inactivas_fva  = setdiff(fva_inactivadas, exch_inactivas_fva);

%% Identify removed metabolites
modelo_N_sparse_reducido   = build_DFVB(model,N_sparse(:,activos_N_sparse));

Inactive_rxn_N_sparse = [modelo_N_sparse_reducido.inac_ex; modelo_N_sparse_reducido.inac_int];
Inactive_rxn_FVA      = model.rxns(fva_inactivadas);

Idx_act_rxn_N_sparse = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_N_sparse));
Idx_act_rxn_FVA      = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_FVA));

Idx_inact_mets_N_sparse   = met_reduced(model, Idx_act_rxn_N_sparse);
Idx_inact_mets_FVA        = met_reduced(model, Idx_act_rxn_FVA);
modelo_FBA_reducido       = Reduce_S_model(model,Idx_act_rxn_FVA);
modelo_FBA_reducido_gurobi= model2gurobi(modelo_FBA_reducido);

fprintf('\nMetabolites removed with FVA = %d\n', size(Idx_inact_mets_FVA,1))
fprintf('\nMetabolites removed with DFVB - SparseNull = %d\n', size(Idx_inact_mets_N_sparse,1))


fprintf('\nReactions removed with FVA = %d\n', size(Inactive_rxn_FVA,1))
fprintf('\nReactions removed with DFVB - SparseNull = %d\n', size(modelo_N_sparse_reducido.inac_ex,1) + size(modelo_N_sparse_reducido.inac_int,1))


%% Run with reduced model

% Índices de reacciones

bio_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_bio'));
glc_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_glc__D_e'));
ac_met_DFVB_N_sparse     = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_ac_e'));
o2_met_DFVB_N_sparse     = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_o2_e'));
co2_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_co2_e'));
nh4_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_nh4_e'));
etoh_met_DFVB_N_sparse   = find(startsWith(modelo_N_sparse_reducido.rxns, 'EX_etoh_e'));

bio_met_DFBA      = findMetIDs(modelo_FBA_reducido, 'biomass_e');
bio_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, 'EX_bio_e');

glc_met_DFBA      = findMetIDs(modelo_FBA_reducido, 'glc__D_e');
glc_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, 'EX_glc__D_e');
glc_ex_fba        = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), 'EX_glc__D_e'));

ac_met_DFBA       = findMetIDs(modelo_FBA_reducido, 'ac_e');
ac_flux_DFBA      = findRxnIDs(modelo_FBA_reducido, 'EX_ac_e');
ac_ex_fba         = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), 'EX_ac_e'));

o2_met_DFBA       = find(startsWith(modelo_FBA_reducido.mets, 'o2_e'));
o2_flux_DFBA      = find(startsWith(modelo_FBA_reducido.rxns, 'EX_o2_e'));
o2_ex_fba         = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), 'EX_o2_e'));

co2_met_DFBA      = findMetIDs(modelo_FBA_reducido, 'co2_e');
co2_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, 'EX_co2_e');
co2_ex_fba        = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), 'EX_co2_e'));

nh4_met_DFBA      = findMetIDs(modelo_FBA_reducido, 'nh4_e');
nh4_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, 'EX_nh4_e');
nh4_ex_fba        = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), 'EX_nh4_e'));

etoh_met_DFBA     = findMetIDs(modelo_FBA_reducido, 'etoh_e');
etoh_flux_DFBA    = findRxnIDs(modelo_FBA_reducido, 'EX_etoh_e');
etoh_ex_fba       = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), 'EX_etoh_e'));

%% Execution of DFVB Real conversions


% ODE configuration
t         = [0, 3];                  % Time interval
y0        = zeros(7,1);              % Initial conditions
options   = odeset('AbsTol',1e-9, 'RelTol',1e-8,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 0.5;                      % Initial biomass condition
y0(2)      = 10;                     % Initial glucose condition
y0(3)      = 0.3125;                     % Initial oxygen condition
y0(4)      = 0;                      % Initial acetate condition
y0(5)      = 27.7;                     % Initial ammonium condition
y0(6)      = 0;                      % Initial ethanol condition
y0(7)      = 0;                      % Initial CO2 condition
params     = [10 0.015 20 0.003 50];                  % Glucose consumption parameters (V_max, K_S)

index_N_sparse   = [bio_met_DFVB_N_sparse, glc_met_DFVB_N_sparse, o2_met_DFVB_N_sparse, ac_met_DFVB_N_sparse, nh4_met_DFVB_N_sparse, etoh_met_DFVB_N_sparse, co2_met_DFVB_N_sparse];
alpha_dim_N_sparse = size(modelo_N_sparse_reducido.nullS,2);

% nullSparse execution
tic; %Start the timer
[TN_sparse,YN_sparse] = ode113(@DFVB_o2_limited,t,y0,options,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);
tiempo_transcurrido_N_sparse = toc; %Save total execution time
disp(['DFVB execution time: ' num2str(tiempo_transcurrido_N_sparse) ' seconds']);

Results_N_sparse_reducido = table(TN_sparse, YN_sparse(:,1), YN_sparse(:,2), YN_sparse(:,3), YN_sparse(:,4), YN_sparse(:,5), YN_sparse(:,6), YN_sparse(:,7), VariableNames=tags);

[alpha_N_sparse_reducido  , fluxes_N_sparse_reducido]   = Variable_recovery(TN_sparse,YN_sparse,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);

Results_fluxes_N_sparse_reducido = table(TN_sparse, fluxes_N_sparse_reducido(:,bio_met_DFVB_N_sparse), fluxes_N_sparse_reducido(:,glc_met_DFVB_N_sparse), fluxes_N_sparse_reducido(:,o2_met_DFVB_N_sparse), fluxes_N_sparse_reducido(:,ac_met_DFVB_N_sparse), fluxes_N_sparse_reducido(:,nh4_met_DFVB_N_sparse), fluxes_N_sparse_reducido(:,etoh_met_DFVB_N_sparse), fluxes_N_sparse_reducido(:,co2_met_DFVB_N_sparse), VariableNames=tags);

%% Average DFVB time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
%sparseNull
tic;
[~,~] = ode113(@DFVB_o2_limited,t,y0,options,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);
t_N_sparse(i) = toc;
end

fprintf('\nMean CPU time of DFVB reduced = %.4f s\n',mean(t_N_sparse))
fprintf('Std CPU time of DFVB reduced = %.4f s\n',std(t_N_sparse)*2)


%% Execution of DFBA


% ODE configuration
t         = [0, 3];                  % Time interval
y0        = ones(67,1)*0;            % Initial conditions
options   = odeset('AbsTol',1e-9, 'RelTol',1e-8,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(bio_met_DFBA) = 0.5;                   % Initial glucose condition
y0(glc_met_DFBA) = 10;                  % Initial biomass condition
y0(o2_met_DFBA)  = 0.3125;              % Initial oxygen condition
y0(ac_met_DFBA)  = 0;                   % Initial acetate condition
y0(co2_met_DFBA) = 0;                   % Initial co2 condition
y0(nh4_met_DFBA) = 27.7;                  % Initial nh4 condition
params     = [10 0.015 20 0.003 50];                  % Glucose consumption parameters (V_max, K_S)
idx_met    = [bio_met_DFBA, glc_met_DFBA, o2_met_DFBA, ac_met_DFBA, co2_met_DFBA, nh4_met_DFBA, etoh_met_DFBA];
idx_flux   = [bio_flux_DFBA, glc_flux_DFBA, o2_flux_DFBA, ac_flux_DFBA, co2_flux_DFBA, nh4_flux_DFBA, etoh_flux_DFBA];

tic; %Start the timer
[T2,Y2] = ode113(@DFBA_o2_limited_gurobi,t,y0,options,modelo_FBA_reducido_gurobi,idx_met,idx_flux,params);

tiempo_transcurrido = toc; %Save total execution time
disp(['DFBA execution time: ' num2str(tiempo_transcurrido) ' seconds']);

%% Average DFBA time calculation
n = 5;
t_DFBA = zeros(n,1);

for i = 1:n
tic;
[~,~] = ode113(@DFBA_o2_limited_gurobi,t,y0,options,modelo_FBA_reducido_gurobi,idx_met,idx_flux,params);
t_DFBA(i) = toc;
end
fprintf('\nMean CPU time of DFBA = %.4f s\n',mean(t_DFBA))
fprintf('Std CPU time of DFBA = %.4f s\n',std(t_DFBA)*2)


%% Comparison plot

% sparseNull
figure();
subplot(1,2,1);
plot(TN_sparse,YN_sparse(:,1),'m','LineWidth',1.5)
hold on
plot(T2,Y2(:,bio_met_DFBA),'b--','LineWidth',1.5)
plot(TN_sparse,YN_sparse(:,2),'r','LineWidth',1.5)
plot(T2,Y2(:,glc_met_DFBA),'g--','LineWidth',1.5)
grid on
legend('Biomass DFVB','Biomass DFBA','Glucose DFVB','Glucose DFVB','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Biomass and Glucose')  

subplot(1,2,2)

yyaxis left
plot(TN_sparse,YN_sparse(:,3),'m',T2,Y2(:,o2_met_DFBA),'b--','LineWidth',1.5)


yyaxis right
hold on
plot(T2, Y2(:,ac_met_DFBA), 'g', 'LineWidth', 1.5)
plot(TN_sparse, YN_sparse(:,4), 'r--', 'LineWidth', 1.5)
xlabel('Time')
legend('DFVB', 'DFBA')
grid on
legend('Oxygen DFVB','Oxygen DFBA','Acetate DFVB','Acetate DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Oxygen and Acetate') 

%% Export data

data_export  = struct();
data_export.Results_N_sparse = Results_N_sparse;
data_export.Results_N_sparse_alpha = alpha_N_sparse;
data_export.modelo_DFVB_N_sparse   = modelo_DFVB_N_sparse;
data_export.Results_N_sparse_reducido    = Results_N_sparse_reducido;
data_export.Results_N_sparse_alpha_reducido     = alpha_N_sparse_reducido;
data_export.modelo_DFVB_N_sparse_reducido = modelo_N_sparse_reducido;
data_export.Results_DFBA = Results_DFBA;
data_export.flux_DFBA = Results_flux_DFBA;
data_export.flux_DFVB_N_sparse = Results_fluxes_N_sparse;
data_export.flux_DFVB_N_sparse_reducido = Results_fluxes_N_sparse_reducido;
data_export.FVA = fva_range;

save('Results\ecolicore\limited_o2.mat', "data_export")
