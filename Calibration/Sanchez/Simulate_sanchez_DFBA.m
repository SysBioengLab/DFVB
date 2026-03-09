%------------------------------------------------------%
%--------  Archivo para simular modelo CyT  -----------%
%------------------------------------------------------%

% initCobraToolbox(false)
changeCobraSolverParams('LP','optTol',1e-8);
changeCobraSolverParams('LP','feasTol',1e-8);

%% Plot configuration

% Plot configuration
set(groot, 'DefaultAxesFontSize', 14); % Tamaño de fuente
set(groot, 'DefaultAxesLineWidth', 1.2); % Grosor de los ejes
set(groot, 'DefaultLineLineWidth', 1.4); % Grosor de las líneas
set(groot, 'DefaultFigureColor', 'w'); % Fondo blanco para la figura
set(groot, 'DefaultAxesFontName', 'Arial'); % Fuente estándar
set(groot, 'DefaultAxesBox', 'on');
set(groot, 'DefaultAxesXColor', 'k');
set(groot, 'DefaultAxesYColor', 'k');
set(groot, 'DefaultFigurePosition', [100, 100, 900, 600]);
output_folder = 'Figures';

color1 = '#f19eec';
color2 = '#8ee1c2';
color3 = '#9190e1';
color4 = '#3c2737';

%% Cargar el modelo

load yeast9-gem.mat
parameters = table2array(readtable("R6_params.csv"));
y0 = parameters(1:6);
theta = parameters(7:end);


%% Find Reactions

rxns       = {'growth', 'D-glucose exchange' ,'ethanol exchange','glycerol exchange', 'citrate exchange', '(S)-lactate exchange', 'oxygen exchange'};
rxns_n2    = {'ammonium exchange'};
FuentesC   = {'acetaldehyde exchange', 'acetate exchange', 'D-galactose exchange' , 'D-galacturonate exchange' ,'D-fructose exchange',...
              'D-mannose exchange', 'ethyl acetate exchange' , 'fumarate exchange', 'succinate exchange', 'pyruvate exchange'};

index_rxns = zeros(size(rxns,2),1);
index_fuentesC = zeros(size(FuentesC,2),1);

for i=1:size(rxns,2)
    index_rxns(i) = find(strcmp(model.rxnNames, rxns(i)));
    mets_index(i) = find(model.S(:, index_rxns(i)));
    
end

exchanges = findExcRxns(model,1);

for i=1:size(FuentesC,2)
    index_fuentesC(i) = find(strcmp(model.rxnNames, FuentesC(i)));
    model = changeRxnBounds(model, model.rxns(index_fuentesC(i)), 0, 'b');
end

%% Update boundaries

model = changeRxnBounds(model, model.rxns(index_rxns(end)), -100, 'l'); % Oxygen exchange
model = changeRxnBounds(model, model.rxns(index_rxns(2)), 1000, 'u'); % Glucose exchange
model = changeRxnBounds(model, model.rxns(index_rxns(3)), 1000, 'u'); % Ethanol exchange
model = changeRxnBounds(model, model.rxns(index_rxns(4)), 1000, 'u'); % Glycerol exchange
model = changeRxnBounds(model, model.rxns(index_rxns(5)), 1000, 'u'); % Citrato exchange
model = changeRxnBounds(model, model.rxns(index_rxns(6)), 1000, 'u'); % Lactato exchange
model = changeRxnBounds(model, model.rxns(index_rxns(5)), 0, 'l'); % Citrato exchange
model = changeRxnBounds(model, model.rxns(index_rxns(6)), 0, 'l'); % Lactato exchange
%% Transform model to gurobi

model_DFBA = model2gurobi(model);

%% Index reactions and metabolites

bio_dfba  = find(startsWith(model_DFBA.rxns, model.rxns(index_rxns(1))));
glc_dfba  = find(startsWith(model_DFBA.rxns, model.rxns(index_rxns(2))));
o2_dfba   = find(startsWith(model_DFBA.rxns, model.rxns(index_rxns(end))));
etoh_dfba = find(startsWith(model_DFBA.rxns, model.rxns(index_rxns(3))));


%% Load Experimental data

exp_data = readtable('Aerobic_data_sanchez_1.csv');

%% DFBA execution: V1

%ODE configuration
options   = odeset('AbsTol',1e-6, 'NonNegative',[1, 2, 3, 4 ,5 ,6]);  % Opciones del ODE
% 
% % Initial conditions
% y0(1)      = exp_data.X(1);   % Initial Biomass
% y0(2)      = exp_data.GLC(2);         % Initial Ammonia
% y0(3)      = exp_data.ETOH(1);              % Initial Glucose
% y0(4)      = exp_data.GLY(2);             % Initial Fructose
% y0(5)      = exp_data.CIT(1);                               % Initial Ethanol
% y0(6)      = exp_data.LAC(1);             % Initial PAN

% Parameters
index_N     = [bio_dfba, glc_dfba, o2_dfba, etoh_dfba];

            %   v_glc  K_G  K_iE  f_etoh  f_gly  f_cit  f_lac  
% theta       =  [15     1    1e6    1.7    4    2.5   3];
t_span = [0, 11];
% Simulate
tic;
[T_sol, Y_sol] = ode113(@model_dFBA2fit_sanchez, t_span,y0,options,theta, model_DFBA, index_N);
fprintf('\nCPU Time(s) = %.2f\n', toc)

% Plot
figure();
subplot(3,2,1)
hold on
plot(T_sol, Y_sol(:,1))
plot(exp_data.TIME, exp_data.X, 'ro')
xlim(t_span)
legend('DFNA', 'Exp.')
ylabel('Concentration (g/L)')
title('Biomass')

subplot(3,2,2)
hold on
plot(T_sol, Y_sol(:,2))
plot(exp_data.TIME, exp_data.GLC, 'ro')
xlim(t_span)
legend('DFNA', 'Exp.')
title('Glucose')
ylabel('Concentration (g/L)')

subplot(3,2,3)
hold on
plot(T_sol, Y_sol(:,3))
plot(exp_data.TIME, exp_data.ETOH, 'ro')
xlim(t_span)
legend('DFNA', 'Exp.')
title('Ethanol')
ylabel('Concentration (g/L)')
% 
subplot(3,2,4)
hold on
plot(T_sol, Y_sol(:,4))
plot(exp_data.TIME, exp_data.GLY, 'ro')
xlim(t_span)
legend('DFNA', 'Exp.')
title('Glycerol')
ylabel('Concentration (g/L)')

subplot(3,2,5)
hold on
plot(T_sol, Y_sol(:,5))
plot(exp_data.TIME, exp_data.CIT, 'ro')
xlim(t_span)
legend('DFNA', 'Exp.')
title('Citrate')
ylabel('Concentration (g/L)')

subplot(3,2,6)
hold on
plot(T_sol, Y_sol(:,6))
plot(exp_data.TIME, exp_data.LAC, 'ro')
xlim(t_span)
legend('DFNA', 'Exp.')
title('Lactate')
ylabel('Concentration (g/L)')

 %% Exportar modelo para calibrar
 
modelo_dfba_struct = struct();

modelo_dfba_struct.model_DFBA  = model_DFBA;
modelo_dfba_struct.index       = index_N;
modelo_dfba_struct.y0          = y0;
modelo_dfba_struct.t_rango1    = t_span;

%% Export results
DFBA_results = table(T_sol, Y_sol(:,1), Y_sol(:,2), ...
    Y_sol(:,3), Y_sol(:,4), Y_sol(:,5), Y_sol(:,6), 'VariableNames',{'Time','Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'} );

% Experiment  = table(t_exp, biomass_exp,  ammonia_exp, glucose_exp, fructose_exp, ethanol_exp, PAN_exp,...
%     'VariableNames',{'Time','Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'} );

folder = 'Results/Calibration CyT/';

% writetable(DFBA_results, [folder,'Model_calibration_DFBA_3.csv'])
% writetable(Experiment, [folder,'Experiment_3.csv'])




