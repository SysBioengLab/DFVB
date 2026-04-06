%------------------------------------------------------%
%-------- File to test DFVB yeast9 --------------------%
%------------------------------------------------------%

changeCobraSolverParams('LP','optTol',1e-8);
changeCobraSolverParams('LP','feasTol',1e-8);

%% Load the model

% initCobraToolbox(false)
load yeast9-GEM.mat

%% Import bases
tic;
N_sparse  = sparseNull(model.S, 1e-8); % Test using matlab function sparseNull
fprintf('\nCPU time (s) of sparseNull = %.4f\n', toc)
fprintf('\n Sparse degree of sparseNull = %.4f%%\n', Sparse_degree(N_sparse)*100)

%% Search exchange names
% Extract reaction names
rxns = {'growth', 'D-glucose exchange' ,'oxygen exchange','acetate exchange', 'carbon dioxide exchange', 'ammonium exchange', 'ethanol exchange'};
exchanges     = model.rxns(findExcRxns(model, true));
idx_exchanges = findRxnIDs(model, exchanges);
rxn_full_names     = model.rxnNames(idx_exchanges);
[~, idx_positions] = ismember(rxns, rxn_full_names);

% Extract coded names
rxn_code_names = exchanges(idx_positions);

% Extract metabolite names

[~, index_rxn_in_S] = ismember(rxn_code_names, model.rxns);
mets_index     = zeros(length(index_rxn_in_S), 1);

for i=1:length(index_rxn_in_S)
    mets_index(i) = find(model.S(:, index_rxn_in_S(i)));
end

mets_names = model.mets(mets_index);
%% Lower and upper bounds 
 
model = changeRxnBounds(model,'r_1992',-1000,'l'); %oxygen exchange
model = changeRxnBounds(model,'r_1761',-0,'l'); %ethanol exchange

%% Compute model with conversions for DFVB

modelo_DFVB_N_sparse   = build_DFVB(model,N_sparse);
modelo_DFBA            = model2gurobi(model);

%% Reaction indices

bio_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(1)));
bio_met_DFBA      = mets_index(1);
bio_flux_DFBA     = findRxnIDs(model, rxn_code_names(1));

glc_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(2)));
glc_met_DFBA      = mets_index(2);
glc_flux_DFBA     = findRxnIDs(model, rxn_code_names(2));
glc_ex_fba        = idx_positions(2);

o2_met_DFVB_N_sparse     = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(3)));
o2_met_DFBA       = mets_index(3);
o2_flux_DFBA      = find(startsWith(model.rxns, rxn_code_names(3)));
o2_ex_fba        = idx_positions(3);

ac_met_DFVB_N_sparse     = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(4)));
ac_met_DFBA       = mets_index(4);
ac_flux_DFBA      = findRxnIDs(model, rxn_code_names(4));
ac_ex_fba         = idx_positions(4);

co2_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(5)));
co2_met_DFBA      = mets_index(5);
co2_flux_DFBA     = findRxnIDs(model, rxn_code_names(5));
co2_ex_fba        = idx_positions(5);

nh4_met_DFVB_N_sparse      = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(6)));
nh4_met_DFBA      = mets_index(6);
nh4_flux_DFBA     = findRxnIDs(model, rxn_code_names(6));
nh4_ex_fba        = idx_positions(6);

etoh_met_DFVB_N_sparse       = find(startsWith(modelo_DFVB_N_sparse.rxns, rxn_code_names(7)));
etoh_met_DFBA     = mets_index(7);
etoh_flux_DFBA    = findRxnIDs(model, rxn_code_names(7));
etoh_ex_fba       = idx_positions(7);

tags = {'Time', 'Biomass', 'Glucose', 'Oxygen', 'Acetate', 'Ammonium', 'Ethanol', 'CO2'};
%% DFVB execution with real conversions

% ODE configuration
t         = [0 3];                  % Time interval
y0        = zeros(7,1);              % Initial conditions
options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 0.5;                    % Initial biomass condition
y0(2)      = 10;                     % Initial glucose condition
y0(3)      = 0.3125;                 % Initial oxygen condition
y0(4)      = 0;                      % Initial acetate condition
y0(5)      = 27.7;                   % Initial ammonium condition
y0(6)      = 0;                      % Initial ethanol condition
y0(7)      = 0;                      % Initial CO2 condition
params     = [10 0.015 20 0.003 100];                  % Glucose consumption parameters (V_max, K_S)

index_N_sparse   = [bio_met_DFVB_N_sparse, glc_met_DFVB_N_sparse, o2_met_DFVB_N_sparse, ac_met_DFVB_N_sparse, nh4_met_DFVB_N_sparse, etoh_met_DFVB_N_sparse, co2_met_DFVB_N_sparse];

alpha_dim_N_sparse = size(modelo_DFVB_N_sparse.nullS,2);

% nullSparse execution
tic; % Start timer
[TN_sparse,YN_sparse] = ode113(@DFVB,t,y0,options,modelo_DFVB_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);
tiempo_transcurrido_N_sparse = toc; % Save total execution time
disp(['DFVB sparse execution time: ' num2str(tiempo_transcurrido_N_sparse) ' seconds']);

Results_N_sparse = table(TN_sparse, YN_sparse(:,1), YN_sparse(:,2), YN_sparse(:,3), YN_sparse(:,4), YN_sparse(:,5), YN_sparse(:,6), YN_sparse(:,7), VariableNames=tags);

%% Average DFVB time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
%sparseNull
tic;
[~,~] = ode113(@DFVB,t,y0,options,modelo_DFVB_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);
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
fprintf('\nConversiones totales de sparseNull = %d', size(alpha_N_sparse,2))
fprintf('\nConversiones usadas con sparseNull = %d\n', size(alpha_N_sparse,2) - length(inactivos_N_sparse))

% %% Gráfico alphas
% figure;
% % sparseNull
% hold on
% for i=1:length(alpha_no_nulo_N_sparse)
%     plot(TN_sparse,alpha_N_sparse(:,alpha_no_nulo_N_sparse(i)), DisplayName= strcat('\alpha',string(alpha_no_nulo_N_sparse(i))),LineWidth=1.2)
% end
% grid on
% ylabel('Actividad de reacción')
% xlabel('Tiempo (h)')
% title('Evolución de alphas usando sparseNull')

%% DFBA execution

% ODE configuration
t         = [0 3];                  % Time interval
y0        = zeros(73,1);          % Initial conditions

% Initial conditions and parameters
y0(bio_met_DFBA) = 0.5;                 % Initial glucose condition
y0(glc_met_DFBA) = 10;                  % Initial biomass condition
y0(o2_met_DFBA)  = 0.3125;              % Initial oxygen condition
y0(ac_met_DFBA)  = 0;                   % Initial acetate condition
y0(co2_met_DFBA) = 0;                   % Initial co2 condition
y0(nh4_met_DFBA) = 27.7;                % Initial nh4 condition
params     = [10 0.015 20 0.003 100];                  % Glucose consumption parameters (V_max, K_S)
idx_met    = [bio_met_DFBA, glc_met_DFBA, o2_met_DFBA, ac_met_DFBA, co2_met_DFBA, nh4_met_DFBA, etoh_met_DFBA];
idx_flux   = [bio_flux_DFBA, glc_flux_DFBA, o2_flux_DFBA, ac_flux_DFBA, co2_flux_DFBA, nh4_flux_DFBA, etoh_flux_DFBA];

options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% DFBA execution (Gurobi warmup)
tic; % Start timer
[TFBA,YFBA] = ode113(@DFBA_gurobi,t,y0,options,modelo_DFBA,idx_met,idx_flux,params);
tiempo_transcurrido = toc; % Save total execution time
disp(['DFBA execution time with warmup: ' num2str(tiempo_transcurrido) ' seconds']);

Results_DFBA = table(TFBA, YFBA(:,bio_met_DFBA), YFBA(:,glc_met_DFBA), YFBA(:,o2_met_DFBA), YFBA(:,ac_met_DFBA), YFBA(:,nh4_met_DFBA), YFBA(:,etoh_met_DFBA), YFBA(:,co2_met_DFBA), VariableNames=tags);

%% DFBA average time calculation
n = 5;
t_DFBA = zeros(n,1);

for i = 1:n
    tic;
    [~,~] = ode113(@DFBA_gurobi,t,y0,options,modelo_DFBA,idx_met,idx_flux,params);
    t_DFBA(i) = toc;
end
fprintf('\nMean CPU time of DFBA = %.4f s\n',mean(t_DFBA))
fprintf('Std CPU time of DFBA = %.4f s\n',std(t_DFBA)*2)

%% Recuperar fluxes FBA

[fluxes_fba, fva_min, fva_max] = Variable_recovery_FBA(TFBA,YFBA,model,idx_met,idx_flux,params, true);
fva_region.fva_min = fva_min;
fva_region.fva_max = fva_max;

% fva_aerobico = fva_aerobico.min;
% fva_max = fva_aerobico.max;
Results_flux_DFBA = table(TFBA, fluxes_fba(bio_flux_DFBA,:)', fluxes_fba(glc_flux_DFBA,:)', fluxes_fba(o2_flux_DFBA,:)', fluxes_fba(ac_flux_DFBA,:)', fluxes_fba(nh4_flux_DFBA,:)', fluxes_fba(etoh_flux_DFBA,:)', fluxes_fba(co2_flux_DFBA,:)', VariableNames=tags);


%% Comparison plot

% sparseNull
figure;
subplot(1,2,1);
plot(TN_sparse,YN_sparse(:,1),'m','LineWidth',1.5)
hold on
plot(TFBA,YFBA(:,bio_met_DFBA),'b--','LineWidth',1.5)
plot(TN_sparse,YN_sparse(:,2),'r','LineWidth',1.5)
plot(TFBA,YFBA(:,glc_met_DFBA),'g--','LineWidth',1.5)
grid on
legend('Biomass DFVB','Biomass DFBA','Glucose DFVB','Glucose DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Biomass and Glucose')  

subplot(1,2,2)
plot(TN_sparse,YN_sparse(:,3),'m',TFBA,YFBA(:,o2_met_DFBA),'b--','LineWidth',1.5)
hold on
plot(TN_sparse,YN_sparse(:,6),'r',TFBA,YFBA(:,etoh_met_DFBA),'g--','LineWidth',1.5)

grid on
legend('Oxygen DFVB','Oxygen DFBA','Ethanol DFVB','Ethanol DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Oxygen and Acetate') 
sgtitle('DFBA and DFVB comparison using sparseNull');

% %% Graficos FVA
% 
% % sparseNull
% figure;
% subplot(2,2,1)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(2,:), fliplr(fva_max(2,:))],'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth',0.1);
% plot(TFBA, fluxes_fba(glc_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,glc_met_DFVB_N_sparse),'--','LineWidth',2)
% grid on
% legend('FVA range','DFBA','DFVB','Location','best')
% ylabel('Consumption rate mmol/mmol s')
% xlabel('Time (h)')
% title('Glucose consumption flux')
% 
% subplot(2,2,2)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(3,:), fliplr(fva_max(3,:))], 'k', 'FaceAlpha', 0.2, 'LineWidth',0.1);
% plot(TFBA, fluxes_fba(o2_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,o2_met_DFVB_N_sparse), '--', 'LineWidth',2)
% grid on
% legend('FVA range','DFBA','DFVB','Location','best')
% % ylabel('Velocidad de consumo mmol/mmol s')
% % ylim([-1 1])
% xlabel('Tiempo (h)')
% title('Flux de consumo de Oxígeno')
% 
% subplot(2,2,3)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(7,:), fliplr(fva_max(7,:))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
% plot(TFBA, fluxes_fba(etoh_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,etoh_met_DFVB_N_sparse),'--','LineWidth',2)
% grid on
% % ylim([-1 1])
% legend('FVA range','DFBA','DFVB','Location','best')
% ylabel('Velocidad de consumo mmol/mmol s')
% xlabel('Tiempo (h)')
% title('Flux de consumo de etanol')
% 
% subplot(2,2,4)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(6,:), fliplr(fva_max(6,:))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth',0.1);
% plot(TFBA, fluxes_fba(nh4_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,nh4_met_DFVB_N_sparse),'--','LineWidth',2)
% grid on
% legend('FVA range','DFBA','DFVB','Location','best')
% % ylabel('Velocidad de consumo mmol/mmol s')
% xlabel('Tiempo (h)')
% title('Flux de consumo de Amonio')
% sgtitle('DFBA and DFVB comparison using sparseNull');


%% Identify inactive
inactivos_N_sparse = setdiff(1:size(alpha_N_sparse,2), alpha_no_nulo_N_sparse);
inactivos_N_sparse = inactivos_N_sparse(2:end);
activos_N_sparse   = setdiff(1:size(modelo_DFVB_N_sparse.nullS,2),inactivos_N_sparse);

% %% Identificar inactivos con FVA
% 
% [min_f, max_f] = fluxVariability(model,99);
% fva_inactivadas = find(min_f == 0 & max_f == 0);
% exch_inactivas_fva = intersect(fva_inactivadas, findRxnIDs(model,exchanges));
% int_inactivas_fva  = setdiff(fva_inactivadas, exch_inactivas_fva);


%% Identify removed metabolites
modelo_N_sparse_reducido   = build_DFVB(model, N_sparse(:,activos_N_sparse));

Inactive_rxn_N_sparse = [modelo_N_sparse_reducido.inac_ex; modelo_N_sparse_reducido.inac_int];
Idx_act_rxn_N_sparse = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_N_sparse));
Idx_inact_mets_N_sparse    = met_reduced(model, Idx_act_rxn_N_sparse);
modelo_FBA_reducido        = Reduce_S_model(model,Idx_act_rxn_N_sparse);
modelo_FBA_reducido_gurobi = model2gurobi(modelo_FBA_reducido);


fprintf('\nMetabolites removed with DFVB - SparseNull = %d\n', size(Idx_inact_mets_N_sparse,1))
fprintf('\nReactions removed with DFVB - SparseNull = %d\n', size(modelo_N_sparse_reducido.inac_ex,1) + size(modelo_N_sparse_reducido.inac_int,1))


%% Calculate MSE for biomass and glucose

% bio_DFVB_N          = YN(:,1);
bio_DFVB_N_sparse   = YN_sparse(:,1);
bio_DFBA            = YFBA(:,bio_met_DFBA);

% Interpolation of DFVB at DFBA times
% bio_DFVB_N_interp = interp1(TN, bio_DFVB_N, TFBA, 'linear', 'extrap');
bio_DFVB_N_sparse_interp = interp1(TN_sparse, bio_DFVB_N_sparse, TFBA, 'linear', 'extrap');

% MSE calculation
% mse_bio_N           = mean((bio_DFBA - bio_DFVB_N_interp).^2);
mse_bio_N_sparse    = sqrt(mean((bio_DFBA - bio_DFVB_N_sparse_interp).^2));

% Display result
disp('MSE Biomass')
% disp(['El MSE entre DFVB-NullS y DFBA es: ', num2str(mse_bio_N)]);
disp(['The MSE between DFVB-sparseNull and DFBA is: ', num2str(mse_bio_N_sparse)]);

%% Run with reduced model

% Reaction indices

bio_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(1)));

glc_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(2)));

o2_met_DFVB_N_sparse     = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(3)));

ac_met_DFVB_N_sparse     = 1; %find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(4)));

co2_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(5)));

nh4_met_DFVB_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(6)));

etoh_met_DFVB_N_sparse   = 2; %find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(7)));

bio_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, rxn_code_names(1));
bio_met_DFBA      = find(modelo_FBA_reducido.S(:, bio_flux_DFBA));

glc_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, rxn_code_names(2));
glc_met_DFBA      = find(modelo_FBA_reducido.S(:, glc_flux_DFBA));
glc_ex_fba        = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(2)));

o2_flux_DFBA      = find(startsWith(modelo_FBA_reducido.rxns, rxn_code_names(3)));
o2_met_DFBA       = find(modelo_FBA_reducido.S(:, o2_flux_DFBA));
o2_ex_fba         = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(3)));

ac_flux_DFBA      = 1; %findRxnIDs(modelo_FBA_reducido, rxn_code_names(4));
ac_met_DFBA       = 1; %find(modelo_FBA_reducido.S(:, ac_flux_DFBA));
ac_ex_fba         = 1; %find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(4)));

co2_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, rxn_code_names(5));
co2_met_DFBA      = find(modelo_FBA_reducido.S(:, co2_flux_DFBA));
co2_ex_fba        = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(5)));

nh4_flux_DFBA     = findRxnIDs(modelo_FBA_reducido, rxn_code_names(6));
nh4_met_DFBA      = find(modelo_FBA_reducido.S(:, nh4_flux_DFBA));
nh4_ex_fba        = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(6)));

etoh_flux_DFBA    = 2; %findRxnIDs(modelo_FBA_reducido, rxn_code_names(7));
etoh_met_DFBA     = 2; %find(modelo_FBA_reducido.S(:, etoh_flux_DFBA));
etoh_ex_fba       = 2; %find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(7)));



%% DFVB execution (real conversions)

% ODE configuration
t         = [0 3];                  % Time interval
y0        = zeros(7,1);              % Initial conditions
options   = odeset('AbsTol',1e-6, 'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 0.5;                      % Initial biomass condition
y0(2)      = 10;                     % Initial glucose condition
y0(3)      = 0.3125;                     % Initial oxygen condition
y0(4)      = 0;                      % Initial acetate condition
y0(5)      = 27.7;                     % Initial ammonium condition
y0(6)      = 0;                      % Initial ethanol condition
y0(7)      = 0;                      % Initial CO2 condition
params     = [10 0.015 20 0.003 100];                  % Glucose consumption parameters (V_max, K_S)
index_N_sparse   = [bio_met_DFVB_N_sparse, glc_met_DFVB_N_sparse, o2_met_DFVB_N_sparse, ac_met_DFVB_N_sparse, nh4_met_DFVB_N_sparse, etoh_met_DFVB_N_sparse, co2_met_DFVB_N_sparse];
alpha_dim_N_sparse = size(modelo_N_sparse_reducido.nullS,2);

% nullSparse execution
tic; % Start timer
[TN_sparse,YN_sparse] = ode113(@DFVB,t,y0,options,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);
tiempo_transcurrido_N_sparse = toc; % Save total execution time
disp(['DFVB sparse execution time: ' num2str(tiempo_transcurrido_N_sparse) ' seconds']);

% Results_N_reducido = table(TN, YN(:,1), YN(:,2), YN(:,3), YN(:,4), YN(:,5), YN(:,6), YN(:,7), VariableNames=tags);
Results_N_sparse_reducido = table(TN_sparse, YN_sparse(:,1), YN_sparse(:,2), YN_sparse(:,3), YN_sparse(:,4), YN_sparse(:,5), YN_sparse(:,6), YN_sparse(:,7), VariableNames=tags);


%% Average DFVB time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
%sparseNull
tic;
[~,~] = ode113(@DFVB,t,y0,options,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);
t_N_sparse(i) = toc;
end

fprintf('\nMean CPU time of DFVB reduced = %.4f s\n',mean(t_N_sparse))
fprintf('Std CPU time of DFVB reduced = %.4f s\n',std(t_N_sparse)*2)



%% DFBA execution

% ODE configuration
t         = [0, 3];                  % Time interval
y0        = ones(73,1)*0;          % Initial conditions
options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(bio_met_DFBA) = 0.5;                   % Initial glucose condition
y0(glc_met_DFBA) = 10;                  % Initial biomass condition
y0(o2_met_DFBA)  = 0.3125;              % Initial oxygen condition
y0(ac_met_DFBA)  = 0;                   % Initial acetate condition
y0(co2_met_DFBA) = 0;                   % Initial co2 condition
y0(nh4_met_DFBA) = 27.7;                  % Initial nh4 condition
params     = [10 0.015 20 0.003 100];                  % Glucose consumption parameters (V_max, K_S)
idx_met    = [bio_met_DFBA, glc_met_DFBA, o2_met_DFBA, ac_met_DFBA, co2_met_DFBA, nh4_met_DFBA, etoh_met_DFBA];
idx_flux   = [bio_flux_DFBA, glc_flux_DFBA, o2_flux_DFBA, ac_flux_DFBA, co2_flux_DFBA, nh4_flux_DFBA, etoh_flux_DFBA];


tic; % Start timer

[TFBA,YFBA] = ode113(@DFBA_gurobi,t,y0,options,modelo_FBA_reducido_gurobi,idx_met,idx_flux,params);
tiempo_transcurrido = toc; % Save total execution time
disp(['DFBA execution time with warmup: ' num2str(tiempo_transcurrido) ' seconds']);

%% Cálculo tiempo promedio DFBA
n = 5;
t_DFBA = zeros(n,1);

for i = 1:n
tic;
[~,~] = ode113(@DFBA_gurobi,t,y0,options,modelo_FBA_reducido_gurobi,idx_met,idx_flux,params);
t_DFBA(i) = toc;
end
fprintf('\nMean CPU time of DFBA = %.4f s\n',mean(t_DFBA))
fprintf('Std CPU time of DFBA = %.4f s\n',std(t_DFBA))


%% Comparison plot

% sparseNull
figure;
subplot(1,2,1);
plot(TN_sparse,YN_sparse(:,1),'m','LineWidth',1.5)
hold on
plot(TFBA,YFBA(:,bio_met_DFBA),'b--','LineWidth',1.5)
plot(TN_sparse,YN_sparse(:,2),'r','LineWidth',1.5)
plot(TFBA,YFBA(:,glc_met_DFBA),'g--','LineWidth',1.5)
grid on
legend('Biomass DFVB','Biomass DFBA','Glucose DFVB','Glucose DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Biomass and Glucose')  

subplot(1,2,2)
plot(TN_sparse,YN_sparse(:,3),'m',TFBA,YFBA(:,o2_met_DFBA),'b--','LineWidth',1.5)
hold on
plot(TN_sparse,YN_sparse(:,6),'r',TFBA,YFBA(:,etoh_met_DFBA),'g--','LineWidth',1.5)

grid on
legend('Oxygen DFVB','Oxygen DFBA','Ethanol DFVB','Ethanol DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Oxygen and Acetate') 
sgtitle('DFBA and DFVB comparison using sparseNull');

%% Calculate MSE for biomass and glucose

% bio_DFVB_N          = YN(:,1);
bio_DFVB_N_sparse   = YN_sparse(:,1);
bio_DFBA            = YFBA(:,bio_met_DFBA);

% Interpolation of DFVB at DFBA times
bio_DFVB_N_sparse_interp = interp1(TN_sparse, bio_DFVB_N_sparse, TFBA, 'linear', 'extrap');

% MSE calculation
mse_bio_N_sparse    = mean((bio_DFBA - bio_DFVB_N_sparse_interp).^2);

% Display result
disp('MSE Biomass')
% disp(['El MSE entre DFVB-NullS y DFBA es: ', num2str(mse_bio_N)]);
disp(['The MSE between DFVB-sparseNull and DFBA is: ', num2str(mse_bio_N_sparse)]);

%% Exportar datos

data_export  = struct();
data_export.Results_N_sparse = Results_N_sparse;
data_export.Results_N_sparse_alpha = alpha_N_sparse;
data_export.Results_N_sparse_reducido    = Results_N_sparse_reducido;
data_export.modelo_DFVB_N_sparse_reducido = modelo_N_sparse_reducido;
data_export.Results_DFBA = Results_DFBA;

save('Results\Yeast9\aerobic.mat', "data_export")
