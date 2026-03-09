%------------------------------------------------------%
%-------- File to test DFNA yeast9 --------------------%
%------------------------------------------------------%

changeCobraSolverParams('LP','optTol',1e-8);
changeCobraSolverParams('LP','feasTol',1e-8);

%% Load the model

% initCobraToolbox(false)
load yeast9-GEM.mat
% load yeast9_anaerobic.mat

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
 
model = changeRxnBounds(model,'r_1992',-10,'l'); %oxygen exchange
model = changeRxnBounds(model,'r_1761',-0,'l'); %ethanol exchange

%% Compute model with conversions for DFNA

modelo_DFNA_N_sparse   = build_DFNA(model,N_sparse);
modelo_DFBA            = model2gurobi(model);

%% Reaction indices

bio_met_DFNA_N_sparse      = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(1)));
bio_met_DFBA      = mets_index(1);
bio_flux_DFBA     = findRxnIDs(model, rxn_code_names(1));

glc_met_DFNA_N_sparse      = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(2)));
glc_met_DFBA      = mets_index(2);
glc_flux_DFBA     = findRxnIDs(model, rxn_code_names(2));
glc_ex_fba        = idx_positions(2);

o2_met_DFNA_N_sparse     = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(3)));
o2_met_DFBA       = mets_index(3);
o2_flux_DFBA      = find(startsWith(model.rxns, rxn_code_names(3)));
o2_ex_fba        = idx_positions(3);

ac_met_DFNA_N_sparse     = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(4)));
ac_met_DFBA       = mets_index(4);
ac_flux_DFBA      = findRxnIDs(model, rxn_code_names(4));
ac_ex_fba         = idx_positions(4);

co2_met_DFNA_N_sparse      = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(5)));
co2_met_DFBA      = mets_index(5);
co2_flux_DFBA     = findRxnIDs(model, rxn_code_names(5));
co2_ex_fba        = idx_positions(5);

nh4_met_DFNA_N_sparse      = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(6)));
nh4_met_DFBA      = mets_index(6);
nh4_flux_DFBA     = findRxnIDs(model, rxn_code_names(6));
nh4_ex_fba        = idx_positions(6);

etoh_met_DFNA_N_sparse       = find(startsWith(modelo_DFNA_N_sparse.rxns, rxn_code_names(7)));
etoh_met_DFBA     = mets_index(7);
etoh_flux_DFBA    = findRxnIDs(model, rxn_code_names(7));
etoh_ex_fba       = idx_positions(7);

tags = {'Time', 'Biomass', 'Glucose', 'Oxygen', 'Acetate', 'Ammonium', 'Ethanol', 'CO2'};
%% DFNA execution with real conversions

% ODE configuration
t         = linspace(0,2,1000);                  % Time interval
y0        = zeros(7,1);              % Initial conditions
options   = odeset('AbsTol',1e-6);  % ODE options

% Initial conditions and parameters
y0(1)      = 1;                      % Initial biomass condition
y0(2)      = 10;                     % Initial glucose condition
y0(3)      = 30;                     % Initial oxygen condition
y0(4)      = 0;                      % Initial acetate condition
y0(5)      = 10;                     % Initial ammonium condition
y0(6)      = 0;                      % Initial ethanol condition
y0(7)      = 0;                      % Initial CO2 condition
params     = [10 0.5];               % Glucose consumption parameters (V_max, K_S)

index_N_sparse   = [bio_met_DFNA_N_sparse, glc_met_DFNA_N_sparse, o2_met_DFNA_N_sparse, ac_met_DFNA_N_sparse, nh4_met_DFNA_N_sparse, etoh_met_DFNA_N_sparse, co2_met_DFNA_N_sparse];

alpha_dim_N_sparse = size(modelo_DFNA_N_sparse.nullS,2);

% nullSparse execution
tic; % Start timer
[TN_sparse,YN_sparse] = ode113(@DFNA,t,y0,options,modelo_DFNA_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);
tiempo_transcurrido_N_sparse = toc; % Save total execution time
disp(['DFNA sparse execution time: ' num2str(tiempo_transcurrido_N_sparse) ' seconds']);

Results_N_sparse = table(TN_sparse, YN_sparse(:,1), YN_sparse(:,2), YN_sparse(:,3), YN_sparse(:,4), YN_sparse(:,5), YN_sparse(:,6), YN_sparse(:,7), VariableNames=tags);

%% DFNA average time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
    % sparseNull
    alpha_no_nulo_N_sparse = [];
    for i=1:size(alpha_N_sparse,2)
        if any(alpha_N_sparse(:,i))
            alpha_no_nulo_N_sparse = [alpha_no_nulo_N_sparse i];   
        end
    end
end
fprintf('Std CPU time of sparseNull = %.4f s\n',std(t_N_sparse))

% Recover non-forced fluxes and alphas

% [alpha_N         , fluxes_N]          = Variable_recovery(TN,YN,modelo_DFNA_N,index_N,params,alpha_dim_N);
[alpha_N_sparse  , fluxes_N_sparse]   = Variable_recovery(TN_sparse,YN_sparse,modelo_DFNA_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);

% Results_fluxes_N        = table(TN, fluxes_N(:,bio_met_DFNA_N), fluxes_N(:,glc_met_DFNA_N), fluxes_N(:,o2_met_DFNA_N), fluxes_N(:,ac_met_DFNA_N), fluxes_N(:,nh4_met_DFNA_N), fluxes_N(:,etoh_met_DFNA_N), fluxes_N(:,co2_met_DFNA_N), VariableNames=tags);
Results_fluxes_N_sparse = table(TN_sparse, fluxes_N_sparse(:,bio_met_DFNA_N_sparse), fluxes_N_sparse(:,glc_met_DFNA_N_sparse), fluxes_N_sparse(:,o2_met_DFNA_N_sparse), fluxes_N_sparse(:,ac_met_DFNA_N_sparse), fluxes_N_sparse(:,nh4_met_DFNA_N_sparse), fluxes_N_sparse(:,etoh_met_DFNA_N_sparse), fluxes_N_sparse(:,co2_met_DFNA_N_sparse), VariableNames=tags);

% % NullS
% alpha_no_nulo_N = [];
% for i=1:size(alpha_N,2)
%     if any(alpha_N(:,i))
%         alpha_no_nulo_N = [alpha_no_nulo_N i];   
%     end
% end
% inactivos_N = setdiff(1:size(alpha_N,2), alpha_no_nulo_N);

% sparseNull
alpha_no_nulo_N_sparse = [];
for i=1:size(alpha_N_sparse,2)
    if any(alpha_N_sparse(:,i))
        alpha_no_nulo_N_sparse = [alpha_no_nulo_N_sparse i];   
    end
end

inactivos_N_sparse = setdiff(1:size(alpha_N_sparse,2), alpha_no_nulo_N_sparse);



% alpha_no_nulo_final = [];
% for i=1:size(alpha_final,2)
%     if any(alpha_final(:,i))
%         alpha_no_nulo_final = [alpha_no_nulo_final i];   
%     end
% end
% 
% inactivos_final= setdiff(1:size(alpha_final,2), alpha_no_nulo_final);

% %%
% 
% % fprintf('\nConversiones totales de NullS = %d', size(alpha_N,2))
% % fprintf('\nConversiones usadas con NullS = %d', size(alpha_N,2) - length(inactivos_N))
fprintf('\nConversiones totales de sparseNull = %d', size(alpha_N_sparse,2))
fprintf('\nConversiones usadas con sparseNull = %d\n', size(alpha_N_sparse,2) - length(inactivos_N_sparse))

% %% Gráfico alphas
% figure;
% %NullS
% % subplot(1,2,1)
% % hold on
% % for i=1:length(alpha_no_nulo_N)
% %     plot(TN,alpha_N(:,alpha_no_nulo_N(i)), DisplayName= strcat('\alpha',string(alpha_no_nulo_N(i))),LineWidth=1.2)
% % end
% % grid on
% % ylabel('Actividad de reacción')
% % xlabel('Tiempo (h)')
% % title('Evolución de alphas usando NullS')
% 
% % sparseNull
% % subplot(1,2,2)
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
t         = linspace(0,2,1000);                  % Time interval
y0        = zeros(73,1);          % Initial conditions

% Initial conditions and parameters
y0(bio_met_DFBA)  = 1;                   % Initial glucose condition
y0(glc_met_DFBA)  = 10;                  % Initial biomass condition
y0(o2_met_DFBA)   = 30;                  % Initial oxygen condition
y0(ac_met_DFBA)   = 0;                   % Initial acetate condition
y0(co2_met_DFBA)  = 0;                   % Initial CO2 condition
y0(nh4_met_DFBA)  = 10;                  % Initial NH4 condition
y0(etoh_met_DFBA) = 0;                  % Initial ethanol condition
params     = [10 0.5];                  % Glucose consumption parameters (V_max, K_S)
idx_met    = [bio_met_DFBA, glc_met_DFBA, o2_met_DFBA, ac_met_DFBA, co2_met_DFBA, nh4_met_DFBA, etoh_met_DFBA];
idx_flux   = [bio_flux_DFBA, glc_flux_DFBA, o2_flux_DFBA, ac_flux_DFBA, co2_flux_DFBA, nh4_flux_DFBA, etoh_flux_DFBA];

options   = odeset('AbsTol',1e-6);  % ODE options

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
fprintf('Std CPU time of DFBA = %.4f s\n',std(t_DFBA))

% %% Recuperar fluxes FBA
% 
% [fluxes_fba, fva_min, fva_max] = Variable_recovery_FBA(TFBA,YFBA,model,idx_met,idx_flux,params, true);
% fva_region.fva_min = fva_min;
% fva_region.fva_max = fva_max;
% 
% % fva_aerobico = fva_aerobico.min;
% % fva_max = fva_aerobico.max;
% Results_flux_DFBA = table(TFBA, fluxes_fba(bio_flux_DFBA,:)', fluxes_fba(glc_flux_DFBA,:)', fluxes_fba(o2_flux_DFBA,:)', fluxes_fba(ac_flux_DFBA,:)', fluxes_fba(nh4_flux_DFBA,:)', fluxes_fba(etoh_flux_DFBA,:)', fluxes_fba(co2_flux_DFBA,:)', VariableNames=tags);
% 

%% Comparison plot

%NullS
% figure;
% subplot(1,2,1);
% plot(TN,YN(:,1),'m','LineWidth',1.5)
% hold on
% plot(TFBA,YFBA(:,bio_met_DFBA),'b--','LineWidth',1.5)
% plot(TN,YN(:,2),'r','LineWidth',1.5)
% plot(TFBA,YFBA(:,glc_met_DFBA),'g--','LineWidth',1.5)
% grid on
% legend('Biomass DFNA','Biomass DFBA','Glucose DFNA','Glucose DFNA','Location','best')
% ylabel('Concentration (mM)')
% xlabel('Time (h)')
% title('Biomass and Glucose')  
% 
% subplot(1,2,2)
% plot(TN,YN(:,3),'m',TFBA,YFBA(:,o2_met_DFBA),'b--','LineWidth',1.5)
% hold on
% plot(TN,YN(:,4),'r',TFBA,YFBA(:,ac_met_DFBA),'g--','LineWidth',1.5)
% grid on
% legend('Oxygen DFNA','Oxygen DFBA','Acetate DFNA','Acetate DFBA','Location','best')
% ylabel('Concentration (mM)')
% xlabel('Time (h)')
% title('Oxygen and Acetate') 
% sgtitle('DFBA and DFNA comparison using NullS');

% sparseNull
figure;
subplot(1,2,1);
plot(TN_sparse,YN_sparse(:,1),'m','LineWidth',1.5)
hold on
plot(TFBA,YFBA(:,bio_met_DFBA),'b--','LineWidth',1.5)
plot(TN_sparse,YN_sparse(:,2),'r','LineWidth',1.5)
plot(TFBA,YFBA(:,glc_met_DFBA),'g--','LineWidth',1.5)
grid on
legend('Biomass DFNA','Biomass DFBA','Glucose DFNA','Glucose DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Biomass and Glucose')  

subplot(1,2,2)
plot(TN_sparse,YN_sparse(:,3),'m',TFBA,YFBA(:,o2_met_DFBA),'b--','LineWidth',1.5)
hold on
plot(TN_sparse,YN_sparse(:,6),'r',TFBA,YFBA(:,etoh_met_DFBA),'g--','LineWidth',1.5)

grid on
legend('Oxygen DFNA','Oxygen DFBA','Ethanol DFNA','Ethanol DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Oxygen and Acetate') 
sgtitle('DFBA and DFNA comparison using sparseNull');

% %% Graficos FVA
% % % NullS
% % figure;
% % subplot(2,2,1)
% % hold on
% % fill([TFBA', fliplr(TFBA')], [fva_min(2,:), fliplr(fva_max(2,:))],'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth',0.1);
% % plot(TFBA, fluxes_fba(glc_flux_DFBA,:),'LineWidth',2)
% % plot(TN, fluxes_N(:,glc_met_DFNA_N),'--','LineWidth',2)
% % grid on
% % legend('FVA range','DFBA','DFNA','Location','best')
% % ylabel('Velocidad de consumo mmol/mmol s')
% % xlabel('Tiempo (h)')
% % title('Flux de consumo de glucosa')
% % 
% % subplot(2,2,2)
% % hold on
% % fill([TFBA', fliplr(TFBA')], [fva_min(3,:), fliplr(fva_max(3,:))], 'k', 'FaceAlpha', 0.2, 'LineWidth',0.1);
% % plot(TFBA, fluxes_fba(o2_flux_DFBA,:),'LineWidth',2)
% % plot(TN, fluxes_N(:,o2_met_DFNA_N), '--', 'LineWidth',2)
% % grid on
% % legend('FVA range','DFBA','DFNA','Location','best')
% % % ylabel('Velocidad de consumo mmol/mmol s')
% % % ylim([-1 1])
% % xlabel('Tiempo (h)')
% % title('Flux de consumo de Oxígeno')
% % 
% % subplot(2,2,3)
% % hold on
% % fill([TFBA', fliplr(TFBA')], [fva_min(4,:), fliplr(fva_max(4,:))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
% % plot(TFBA, fluxes_fba(ac_flux_DFBA,:),'LineWidth',2)
% % plot(TN, fluxes_N(:,ac_met_DFNA_N),'--','LineWidth',2)
% % grid on
% % % ylim([-1 1])
% % legend('FVA range','DFBA','DFNA','Location','best')
% % ylabel('Velocidad de consumo mmol/mmol s')
% % xlabel('Tiempo (h)')
% % title('Flux de consumo de Acetato')
% % 
% % subplot(2,2,4)
% % hold on
% % fill([TFBA', fliplr(TFBA')], [fva_min(6,:), fliplr(fva_max(6,:))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth',0.1);
% % plot(TFBA, fluxes_fba(nh4_flux_DFBA,:),'LineWidth',2)
% % plot(TN, fluxes_N(:,nh4_met_DFNA_N),'--','LineWidth',2)
% % grid on
% % legend('FVA range','DFBA','DFNA','Location','best')
% % % ylabel('Velocidad de consumo mmol/mmol s')
% % xlabel('Tiempo (h)')
% % title('Flux de consumo de Amonio')
% % sgtitle('DFBA and DFNA comparison using NullS');
% 
% % sparseNull
% figure;
% subplot(2,2,1)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(2,:), fliplr(fva_max(2,:))],'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth',0.1);
% plot(TFBA, fluxes_fba(glc_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,glc_met_DFNA_N_sparse),'--','LineWidth',2)
% grid on
% legend('FVA range','DFBA','DFNA','Location','best')
% ylabel('Consumption rate mmol/mmol s')
% xlabel('Time (h)')
% title('Glucose consumption flux')
% 
% subplot(2,2,2)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(3,:), fliplr(fva_max(3,:))], 'k', 'FaceAlpha', 0.2, 'LineWidth',0.1);
% plot(TFBA, fluxes_fba(o2_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,o2_met_DFNA_N_sparse), '--', 'LineWidth',2)
% grid on
% legend('FVA range','DFBA','DFNA','Location','best')
% % ylabel('Velocidad de consumo mmol/mmol s')
% % ylim([-1 1])
% xlabel('Tiempo (h)')
% title('Flux de consumo de Oxígeno')
% 
% subplot(2,2,3)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(7,:), fliplr(fva_max(7,:))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
% plot(TFBA, fluxes_fba(etoh_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,etoh_met_DFNA_N_sparse),'--','LineWidth',2)
% grid on
% % ylim([-1 1])
% legend('FVA range','DFBA','DFNA','Location','best')
% ylabel('Velocidad de consumo mmol/mmol s')
% xlabel('Tiempo (h)')
% title('Flux de consumo de etanol')
% 
% subplot(2,2,4)
% hold on
% fill([TFBA', fliplr(TFBA')], [fva_min(6,:), fliplr(fva_max(6,:))], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'LineWidth',0.1);
% plot(TFBA, fluxes_fba(nh4_flux_DFBA,:),'LineWidth',2)
% plot(TN_sparse, fluxes_N_sparse(:,nh4_met_DFNA_N_sparse),'--','LineWidth',2)
% grid on
% legend('FVA range','DFBA','DFNA','Location','best')
% % ylabel('Velocidad de consumo mmol/mmol s')
% xlabel('Tiempo (h)')
% title('Flux de consumo de Amonio')
% sgtitle('DFBA and DFNA comparison using sparseNull');


%% Identify inactive
% inactivos_N = setdiff(1:size(alpha_N,2), alpha_no_nulo_N);
inactivos_N_sparse = setdiff(1:size(alpha_N_sparse,2), alpha_no_nulo_N_sparse);
inactivos_N_sparse = inactivos_N_sparse(2:end);
% activos_N          = setdiff(1:size(modelo_DFNA_N.nullS,2),inactivos_N);
activos_N_sparse   = setdiff(1:size(modelo_DFNA_N_sparse.nullS,2),inactivos_N_sparse);

% %% Identificar inactivos con FVA
% 
% [min_f, max_f] = fluxVariability(model,99);
% fva_inactivadas = find(min_f == 0 & max_f == 0);
% exch_inactivas_fva = intersect(fva_inactivadas, findRxnIDs(model,exchanges));
% int_inactivas_fva  = setdiff(fva_inactivadas, exch_inactivas_fva);


%% Identify removed metabolites
% modelo_N_reducido          = build_DFNA(model,N(:,activos_N));
modelo_N_sparse_reducido   = build_DFNA(model, N_sparse(:,activos_N_sparse));

% Inactive_rxn_N        = [modelo_N_reducido.inac_ex; modelo_N_reducido.inac_int];
Inactive_rxn_N_sparse = [modelo_N_sparse_reducido.inac_ex; modelo_N_sparse_reducido.inac_int];
% Inactive_rxn_FVA      = model.rxns(fva_inactivadas);

% Idx_act_rxn_N        = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_N));
Idx_act_rxn_N_sparse = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_N_sparse));
% Idx_act_rxn_FVA      = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_FVA));

% Idx_inact_mets_N           = met_reduced(model, Idx_act_rxn_N);
Idx_inact_mets_N_sparse    = met_reduced(model, Idx_act_rxn_N_sparse);
% Idx_inact_mets_FVA         = met_reduced(model, Idx_act_rxn_FVA);
modelo_FBA_reducido        = Reduce_S_model(model,Idx_act_rxn_N_sparse);
modelo_FBA_reducido_gurobi = model2gurobi(modelo_FBA_reducido);


% fprintf('\nMetabolitos eliminados con FVA = %d\n', size(Idx_inact_mets_FVA,1))
% fprintf('\nMetabolitos eliminados con DFNA - NullS = %d\n', size(Idx_inact_mets_N,1))
fprintf('\nMetabolites removed with DFNA - SparseNull = %d\n', size(Idx_inact_mets_N_sparse,1))


% fprintf('\nReacciones eliminadas con FVA = %d\n', size(Inactive_rxn_FVA,1))
% fprintf('\nReacciones eliminadas con DFNA - NullS = %d\n', size(modelo_N_reducido.inac_ex,1) + size(modelo_N_reducido.inac_int,1))
fprintf('\nReactions removed with DFNA - SparseNull = %d\n', size(modelo_N_sparse_reducido.inac_ex,1) + size(modelo_N_sparse_reducido.inac_int,1))


%% Calculate MSE for biomass and glucose

% bio_DFNA_N          = YN(:,1);
bio_DFNA_N_sparse   = YN_sparse(:,1);
bio_DFBA            = YFBA(:,bio_met_DFBA);

% Interpolation of DFNA at DFBA times
% bio_DFNA_N_interp = interp1(TN, bio_DFNA_N, TFBA, 'linear', 'extrap');
bio_DFNA_N_sparse_interp = interp1(TN_sparse, bio_DFNA_N_sparse, TFBA, 'linear', 'extrap');

% MSE calculation
% mse_bio_N           = mean((bio_DFBA - bio_DFNA_N_interp).^2);
mse_bio_N_sparse    = sqrt(mean((bio_DFBA - bio_DFNA_N_sparse_interp).^2));

% Display result
disp('MSE Biomass')
% disp(['El MSE entre DFNA-NullS y DFBA es: ', num2str(mse_bio_N)]);
disp(['The MSE between DFNA-sparseNull and DFBA is: ', num2str(mse_bio_N_sparse)]);

%% Run with reduced model

% Reaction indices

bio_met_DFNA_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(1)));

glc_met_DFNA_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(2)));

o2_met_DFNA_N_sparse     = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(3)));

ac_met_DFNA_N_sparse     = 1; %find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(4)));

co2_met_DFNA_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(5)));

nh4_met_DFNA_N_sparse    = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(6)));

etoh_met_DFNA_N_sparse   = find(startsWith(modelo_N_sparse_reducido.rxns, rxn_code_names(7)));

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

etoh_flux_DFBA    = findRxnIDs(modelo_FBA_reducido, rxn_code_names(7));
etoh_met_DFBA     = find(modelo_FBA_reducido.S(:, etoh_flux_DFBA));
etoh_ex_fba       = find(strcmp(modelo_FBA_reducido.rxns(findExcRxns(modelo_FBA_reducido)), rxn_code_names(7)));



%% DFNA execution (real conversions)

% ODE configuration
t         = linspace(0,2,1000);                  % Time interval
y0        = zeros(7,1);              % Initial conditions
options   = odeset('AbsTol',1e-6);  % ODE options

% Initial conditions and parameters
y0(1)      = 1;                      % Initial biomass condition
y0(2)      = 10;                     % Initial glucose condition
y0(3)      = 30;                     % Initial oxygen condition
y0(4)      = 0;                      % Initial acetate condition
y0(5)      = 10;                     % Initial ammonium condition
y0(6)      = 0;                      % Initial ethanol condition
y0(7)      = 0;                      % Initial CO2 condition
params     = [10 0.5];               % Glucose consumption parameters (V_max, K_S)
% index_N          = [bio_met_DFNA_N, glc_met_DFNA_N, o2_met_DFNA_N, ac_met_DFNA_N, nh4_met_DFNA_N, etoh_met_DFNA_N, co2_met_DFNA_N];
index_N_sparse   = [bio_met_DFNA_N_sparse, glc_met_DFNA_N_sparse, o2_met_DFNA_N_sparse, ac_met_DFNA_N_sparse, nh4_met_DFNA_N_sparse, etoh_met_DFNA_N_sparse, co2_met_DFNA_N_sparse];

% alpha_dim_N = size(modelo_N_reducido.nullS,2);
alpha_dim_N_sparse = size(modelo_N_sparse_reducido.nullS,2);

% NullS execution
% tic; % Start timer
% [TN,YN] = ode113(@DFNA,t,y0,options,modelo_N_reducido,index_N,params,alpha_dim_N);
% tiempo_transcurrido_N = toc; % Save total execution time
% disp(['DFNA execution time: ' num2str(tiempo_transcurrido_N) ' seconds']);

% nullSparse execution
tic; % Start timer
[TN_sparse,YN_sparse] = ode113(@DFNA,t,y0,options,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);
tiempo_transcurrido_N_sparse = toc; % Save total execution time
disp(['DFNA sparse execution time: ' num2str(tiempo_transcurrido_N_sparse) ' seconds']);

% Results_N_reducido = table(TN, YN(:,1), YN(:,2), YN(:,3), YN(:,4), YN(:,5), YN(:,6), YN(:,7), VariableNames=tags);
Results_N_sparse_reducido = table(TN_sparse, YN_sparse(:,1), YN_sparse(:,2), YN_sparse(:,3), YN_sparse(:,4), YN_sparse(:,5), YN_sparse(:,6), YN_sparse(:,7), VariableNames=tags);

% [alpha_N_reducido  , fluxes_N_reducido]          = Variable_recovery(TN,YN,modelo_DFNA_N,index_N,params,alpha_dim_N);
% [alpha_N_sparse_reducido  , fluxes_N_sparse_reducido]   = Variable_recovery(TN_sparse,YN_sparse,modelo_DFNA_N_sparse,index_N_sparse,params,alpha_dim_N_sparse);

% Results_fluxes_N_reducido        = table(TN, fluxes_N_reducido(:,bio_met_DFNA_N), fluxes_N_reducido(:,glc_met_DFNA_N), fluxes_N_reducido(:,o2_met_DFNA_N), fluxes_N_reducido(:,ac_met_DFNA_N), fluxes_N_reducido(:,nh4_met_DFNA_N), fluxes_N_reducido(:,etoh_met_DFNA_N), fluxes_N_reducido(:,co2_met_DFNA_N), VariableNames=tags);
% Results_fluxes_N_sparse_reducido = table(TN_sparse, fluxes_N_sparse_reducido(:,bio_met_DFNA_N_sparse), fluxes_N_sparse_reducido(:,glc_met_DFNA_N_sparse), fluxes_N_sparse_reducido(:,o2_met_DFNA_N_sparse), fluxes_N_sparse_reducido(:,ac_met_DFNA_N_sparse), fluxes_N_sparse_reducido(:,nh4_met_DFNA_N_sparse), fluxes_N_sparse_reducido(:,etoh_met_DFNA_N_sparse), fluxes_N_sparse_reducido(:,co2_met_DFNA_N_sparse), VariableNames=tags);


%% Save
data_export  = struct();
data_export.Results_N_sparse = Results_N_sparse;
data_export.Results_DFBA = Results_DFBA;
data_export.Results_N_sparse_reducido    = Results_N_sparse_reducido;

%% DFNA average time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
   
    % sparseNull timing
    tic;
    [~,~] = ode113(@DFNA,t,y0,options,modelo_N_sparse_reducido,index_N_sparse,params,alpha_dim_N_sparse);
    t_N_sparse(i) = toc;
end

fprintf('\nMean CPU time of sparseNull = %.4f s\n',mean(t_N_sparse))
fprintf('Std CPU time of sparseNull = %.4f s\n',std(t_N_sparse))


%% DFBA execution

% ODE configuration
t         = [0, 2];                  % Time interval
y0        = ones(73,1)*0;          % Initial conditions
options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(bio_met_DFBA) = 1;                   % Initial biomass condition
y0(glc_met_DFBA) = 10;                  % Initial glucose condition
y0(o2_met_DFBA)  = 30;                  % Initial oxygen condition
y0(ac_met_DFBA)  = 0;                   % Initial acetate condition
y0(co2_met_DFBA) = 0;                   % Initial CO2 condition
y0(nh4_met_DFBA) = 10;                  % Initial ammonium condition
y0(etoh_met_DFBA)= 0;
params     = [10 0.5];                  % Glucose consumption parameters (V_max, K_S)
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

% %NullS
% figure;
% subplot(1,2,1);
% plot(TN,YN(:,1),'m','LineWidth',1.5)
% hold on
% plot(TFBA,YFBA(:,bio_met_DFBA),'b--','LineWidth',1.5)
% plot(TN,YN(:,2),'r','LineWidth',1.5)
% plot(TFBA,YFBA(:,glc_met_DFBA),'g--','LineWidth',1.5)
% grid on
% legend('Biomass DFNA','Biomass DFBA','Glucose DFNA','Glucose DFNA','Location','best')
% ylabel('Concentration (mM)')
% xlabel('Time (h)')
% title('Biomass and Glucose')  
% 
% subplot(1,2,2)
% plot(TN,YN(:,3),'m',TFBA,YFBA(:,o2_met_DFBA),'b--','LineWidth',1.5)
% hold on
% plot(TN,YN(:,4),'r',TFBA,YFBA(:,ac_met_DFBA),'g--','LineWidth',1.5)
% grid on
% legend('Oxygen DFNA','Oxygen DFBA','Acetate DFNA','Acetate DFBA','Location','best')
% ylabel('Concentration (mM)')
% xlabel('Time (h)')
% title('Oxygen and Acetate') 
% sgtitle('DFBA and DFNA comparison using NullS');

% sparseNull
figure;
subplot(1,2,1);
plot(TN_sparse,YN_sparse(:,1),'m','LineWidth',1.5)
hold on
plot(TFBA,YFBA(:,bio_met_DFBA),'b--','LineWidth',1.5)
plot(TN_sparse,YN_sparse(:,2),'r','LineWidth',1.5)
plot(TFBA,YFBA(:,glc_met_DFBA),'g--','LineWidth',1.5)
grid on
legend('Biomass DFNA','Biomass DFBA','Glucose DFNA','Glucose DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Biomass and Glucose')  

subplot(1,2,2)
plot(TN_sparse,YN_sparse(:,3),'m',TFBA,YFBA(:,o2_met_DFBA),'b--','LineWidth',1.5)
hold on
plot(TN_sparse,YN_sparse(:,6),'r',TFBA,YFBA(:,etoh_met_DFBA),'g--','LineWidth',1.5)

grid on
legend('Oxygen DFNA','Oxygen DFBA','Ethanol DFNA','Ethanol DFBA','Location','best')
ylabel('Concentration (mM)')
xlabel('Time (h)')
title('Oxygen and Acetate') 
sgtitle('DFBA and DFNA comparison using sparseNull');

%% Calculate MSE for biomass and glucose

% bio_DFNA_N          = YN(:,1);
bio_DFNA_N_sparse   = YN_sparse(:,1);
bio_DFBA            = YFBA(:,bio_met_DFBA);

% Interpolation of DFNA at DFBA times
% bio_DFNA_N_interp = interp1(TN, bio_DFNA_N, TFBA, 'linear', 'extrap');
bio_DFNA_N_sparse_interp = interp1(TN_sparse, bio_DFNA_N_sparse, TFBA, 'linear', 'extrap');

% MSE calculation
% mse_bio_N           = mean((bio_DFBA - bio_DFNA_N_interp).^2);
mse_bio_N_sparse    = mean((bio_DFBA - bio_DFNA_N_sparse_interp).^2);

% Display result
disp('MSE Biomass')
% disp(['El MSE entre DFNA-NullS y DFBA es: ', num2str(mse_bio_N)]);
disp(['The MSE between DFNA-sparseNull and DFBA is: ', num2str(mse_bio_N_sparse)]);

% %% Exportar datos
% 
% Results_yeast9_DFNA_aerobico  = struct();
% Results_yeast9_DFNA_aerobico.Results_N        = Results_N;
% Results_yeast9_DFNA_aerobico.Results_N_alpha  = alpha_N;
% Results_yeast9_DFNA_aerobico.modelo_DFNA_N    = modelo_DFNA_N;
% Results_yeast9_DFNA_aerobico.Results_N_sparse = Results_N_sparse;
% Results_yeast9_DFNA_aerobico.Results_N_sparse_alpha = alpha_N_sparse;
% Results_yeast9_DFNA_aerobico.modelo_DFNA_N_sparse   = modelo_DFNA_N_sparse;
% Results_yeast9_DFNA_aerobico.Results_N_reducido     = Results_N_reducido;
% Results_yeast9_DFNA_aerobico.Results_N_alpha_reducido     = alpha_N_reducido;
% Results_yeast9_DFNA_aerobico.modelo_DFNA_N_reducido = modelo_N_reducido;
% Results_yeast9_DFNA_aerobico.Results_N_sparse_reducido    = Results_N_sparse_reducido;
% Results_yeast9_DFNA_aerobico.Results_N_sparse_alpha_reducido     = alpha_N_sparse_reducido;
% Results_yeast9_DFNA_aerobico.modelo_DFNA_N_sparse_reducido = modelo_N_sparse_reducido;
% Results_yeast9_DFNA_aerobico.Results_DFBA = Results_DFBA;
% Results_yeast9_DFNA_aerobico.flux_DFBA = Results_flux_DFBA;
% Results_yeast9_DFNA_aerobico.flux_DFNA_N = Results_fluxes_N;
% Results_yeast9_DFNA_aerobico.flux_DFNA_N_sparse = Results_fluxes_N_sparse;
% Results_yeast9_DFNA_aerobico.flux_DFNA_N_reducido = Results_fluxes_N_reducido;
% Results_yeast9_DFNA_aerobico.flux_DFNA_N_sparse_reducido = Results_fluxes_N_sparse_reducido;
% Results_yeast9_DFNA_aerobico.FVA = fva_region;
