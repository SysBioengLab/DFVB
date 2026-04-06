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

%% Cargar el modelo, data experimental y parametros
load yeast9_anaerobic.mat
load('data_exp_R1.mat');
theta = table2array(readtable('Replicado1_params.csv'));

%% Composición de fuente de nitrógeno
PAN  = [find(strcmp(model.rxnNames, 'L-arginine exchange'));
        find(strcmp(model.rxnNames, 'L-glutamine exchange'));
        find(strcmp(model.rxnNames, 'L-leucine exchange'));
        find(strcmp(model.rxnNames, 'L-valine exchange'));
        find(strcmp(model.rxnNames, 'L-serine exchange'));
        find(strcmp(model.rxnNames, 'L-alanine exchange'))];

PAN_prop = [0.30, ... % Arginine
            0.25, ... % Glutamine
            0.15, ... % Leucine
            0.15, ... % Valine
            0.10, ... % Serine
            0.05];    % Alanine

MW_PAN =   [174.2, ... % Arginine
            146.1, ... % Glutamine
            131.2, ... % Leucine
            117.1, ... % Valine
            105.1, ... % Serine
            89.1];     % Alanine

model_aa_irrev = make_aa_irrev(model, PAN);

%% Find Reactions

rxns       = {'growth', 'D-glucose exchange' ,'oxygen exchange','D-fructose exchange', 'ammonium exchange', 'ethanol exchange'};
index_rxns = zeros(size(rxns,2),1);
mets_index = zeros(size(rxns,2),1);

for i=1:size(rxns,2)
    index_rxns(i) = find(strcmp(model_aa_irrev.rxnNames, rxns(i)));
    mets_index(i) = find(model_aa_irrev.S(:, index_rxns(i)));
end

exchanges = findExcRxns(model_aa_irrev,1);

%% Update boundaries

model_aa_irrev = changeRxnBounds(model_aa_irrev, model_aa_irrev.rxns(index_rxns(3)), 0, 'b'); % Oxygen exchange
model_aa_irrev = changeRxnBounds(model_aa_irrev, model_aa_irrev.rxns(index_rxns(2)), -1000, 'l'); % Glucose exchange
model_aa_irrev = changeRxnBounds(model_aa_irrev, model_aa_irrev.rxns(index_rxns(4)), -1000, 'l'); % Fructose exchange


%% Construct DFVB and DFBA models
N_sparse   = sparseNull(model_aa_irrev.S);
model_DFVB = build_DFVB(model_aa_irrev, N_sparse);

%% Index reactions and metabolites

bio_DFVB  = find(startsWith(model_DFVB.rxns, model_aa_irrev.rxns(index_rxns(1))));
glc_DFVB  = find(startsWith(model_DFVB.rxns, model_aa_irrev.rxns(index_rxns(2))));
fru_DFVB  = find(startsWith(model_DFVB.rxns, model_aa_irrev.rxns(index_rxns(4))));
n_DFVB    = find(startsWith(model_DFVB.rxns, model_aa_irrev.rxns(index_rxns(5))));
etoh_DFVB = find(startsWith(model_DFVB.rxns, model_aa_irrev.rxns(index_rxns(6))));
pan_DFVB  = find(startsWith(model_DFVB.rxns, model_aa_irrev.rxns(PAN)));

PAN_struct = struct('index_DFVB', pan_DFVB(end-5:end), 'PAN_prop', PAN_prop, 'MW_PAN', MW_PAN);

%% Add PAN restriction

A_row = zeros(1, size(model_DFVB.A,2));
A_row(1, pan_DFVB(end-5:end)) = ones(1,6);
b_row                         = 0;
sense_new                     = '<'; 

model_DFVB.A     = [model_DFVB.A; A_row];
model_DFVB.rhs   = [model_DFVB.rhs; b_row];
model_DFVB.sense = [model_DFVB.sense; sense_new];

%% Load Experimental data

% Temp_1         = readtable('Temperature_R2.csv');
% Conc_R1        = readtable('Replicado_R2.csv', 'VariableNamingRule','preserve');
biomass_factor = 0.035307074; % 10^6 cells/mL -> gDW/L

t_exp        = data_exp.t_exp;
t_rango1     = data_exp.t_rango1;
t_rango2     = data_exp.t_rango2;
Temperature  = data_exp.Temperature;
biomass_exp  = data_exp.conc(:,1);
ammonia_exp  = data_exp.conc(:,2);
glucose_exp  = data_exp.conc(:,3);
fructose_exp = data_exp.conc(:,4);
ethanol_exp  = data_exp.conc(:,5);
PAN_exp      = data_exp.conc(:,6);


%% DFVB execution

%ODE configuration
options   = odeset('AbsTol',1e-6, 'NonNegative',[1, 2, 3, 4, 5, 6]);  % Opciones del ODE

% Initial conditions
y0(1)      = biomass_exp(1);         % Initial Biomass
y0(2)      = ammonia_exp(1);         % Initial Ammonia
y0(3)      = glucose_exp(2);         % Initial Glucose
y0(4)      = fructose_exp(2);        % Initial Fructose
y0(5)      = 0;                      % Initial Ethanol
y0(6)      = PAN_exp(1);             % Initial PAN

% Parameters
index_N     = [bio_DFVB, glc_DFVB, fru_DFVB, n_DFVB, etoh_DFVB];
alpha_dim   = size(model_DFVB.nullS,2);

% Simulate first period
tic;
[T1, Y1] = ode113(@model_DFVB2fit,t_rango1,theta(1:6),options,theta(7:end),Temperature, model_DFVB, index_N, PAN_struct, alpha_dim);

% Simulate second period
y0_2     = Y1(end,:);
y0_2(2)  = ammonia_exp(6);
[T2, Y2] = ode113(@model_DFVB2fit,t_rango2,y0_2,options,theta(7:end),Temperature, model_DFVB, index_N, PAN_struct, alpha_dim);
fprintf('\nCPU Time(s) = %.2f\n', toc)

% Join the solution
T_sol = [T1; T2];
Y_sol = [Y1; Y2];

%% Plot

figure;
subplot(3,2,1)
hold on
plot(T_sol, Y_sol(:,1))
plot(t_exp, biomass_exp, 'ro')
% xlim([0 70])
legend('DFVB', 'Exp.')
ylabel('Concentration (mM)')
title('Biomass')

subplot(3,2,2)
hold on
plot(T_sol, Y_sol(:,3))
plot(t_exp, glucose_exp, 'ro')
% xlim([0 70])
legend('DFVB', 'Exp.')
title('Glucose')
ylabel('Concentration (g/L)')

subplot(3,2,3)
hold on
plot(T_sol, Y_sol(:,4))
plot(t_exp, fructose_exp, 'ro')
% xlim([0 70])
legend('DFVB', 'Exp.')
title('Fructose')
ylabel('Concentration (g/L)')

subplot(3,2,4)
hold on
plot(T_sol, Y_sol(:,5))
plot(t_exp, ethanol_exp, 'ro')
% xlim([0 70])
legend('DFVB', 'Exp.')
title('Ethanol')
ylabel('Concentration (mg/L)')

subplot(3,2,5)
hold on
plot(T_sol, Y_sol(:,2))
plot(t_exp, ammonia_exp, 'ro')
% xlim([0 70])
legend('DFVB', 'Exp.')
title('Ammonia')
ylabel('Concentration (mg/L)')

subplot(3,2,6)
hold on
plot(T_sol, Y_sol(:,6))
plot(t_exp, PAN_exp, 'ro')
% xlim([0 70])
legend('DFVB', 'Exp.')
title('PAN')
ylabel('Concentration (mg/L)')

%% Recuperar fluxes y alphas no forzado

alpha_N = alpha_recovery(T_sol,Y_sol,theta(7:end),Temperature, model_DFVB, index_N, PAN_struct, alpha_dim);

% NullS
alpha_no_nulo_N = [];
for i=1:size(alpha_N,2)
    if any(alpha_N(:,i))
        alpha_no_nulo_N = [alpha_no_nulo_N i];   
    end
end

%% Reducir modelo
inactivos_N        = setdiff(1:size(alpha_N,2), alpha_no_nulo_N);
activos_N_sin_crec = setdiff(1:size(model_DFVB.nullS,2),inactivos_N); % POS 1 NO SE PUEDE ELIMINAR -> COL 1
activos_N          = setdiff(1:size(model_DFVB.nullS,2),[inactivos_N(2:end)]); % POS 1 NO SE PUEDE ELIMINAR -> COL 1

modelo_N_reducido_sin_crec = build_DFVB(model_aa_irrev,N_sparse(:,activos_N_sin_crec));
modelo_N_reducido          = build_DFVB(model_aa_irrev,N_sparse(:,activos_N));

% % Excepciones: no borrar conversiones que tengan aminoacidos
% no_borrar       = PAN_struct.index_DFVB;
% cols_inactivas  = inactivos_N(2:end);
% tienen_valores  = any(N_sparse(no_borrar, cols_inactivas) ~= 0, 1);
% cols_a_eliminar = cols_inactivas(~tienen_valores);
% activos_N       = setdiff(1:size(model_DFVB.nullS,2), cols_a_eliminar);
% modelo_N_reducido = build_DFVB(model_aa_irrev,N_sparse(:, activos_N));

%% Identificar reacciones y metabolitos eliminados
Inactive_rxn_N        = [modelo_N_reducido.inac_ex; modelo_N_reducido.inac_int];
Idx_inact_rxn_N       = findRxnIDs(model_aa_irrev, Inactive_rxn_N);
Idx_act_rxn_N         = setdiff(1:size(model_aa_irrev.S,2) ,findRxnIDs(model_aa_irrev, Inactive_rxn_N));
Idx_inact_mets_N      = met_reduced(model_aa_irrev, Idx_act_rxn_N);
fprintf('\nMetabolitos eliminados con DFVB - SparseNull = %d\n', size(Idx_inact_mets_N,1))
fprintf('\nReacciones eliminadas con DFVB - SparseNull = %d\n', size(modelo_N_reducido.inac_ex,1) + size(modelo_N_reducido.inac_int,1))

% New indexes of reactions and metabolites

bio_DFVB_2  = find(startsWith(modelo_N_reducido.rxns, model_aa_irrev.rxns(index_rxns(1))));
glc_DFVB_2  = find(startsWith(modelo_N_reducido.rxns, model_aa_irrev.rxns(index_rxns(2))));
fru_DFVB_2  = find(startsWith(modelo_N_reducido.rxns, model_aa_irrev.rxns(index_rxns(4))));
n_DFVB_2    = find(startsWith(modelo_N_reducido.rxns, model_aa_irrev.rxns(index_rxns(5))));
etoh_DFVB_2 = find(startsWith(modelo_N_reducido.rxns, model_aa_irrev.rxns(index_rxns(6))));
pan_DFVB_2  = find(startsWith(modelo_N_reducido.rxns, model_aa_irrev.rxns(PAN)));

PAN_struct_2 = struct('index_DFVB', pan_DFVB_2, 'PAN_prop', PAN_prop, 'MW_PAN', MW_PAN);

%% Corregir restriccion PAN
A_row = zeros(1, size(modelo_N_reducido.A,2));
A_row(1, pan_DFVB_2) = ones(1,1);
b_row                         = 0;
sense_new                     = '<'; 

modelo_N_reducido.A     = [modelo_N_reducido.A; A_row];
modelo_N_reducido.rhs   = [modelo_N_reducido.rhs; b_row];
modelo_N_reducido.sense = [modelo_N_reducido.sense; sense_new];

%% DFVB execution

%ODE configuration
options   = odeset('AbsTol',1e-6, 'NonNegative',[1, 2, 3, 4, 5, 6]);  % Opciones del ODE

% Initial conditions
y0(1)      = biomass_exp(1);         % Initial Biomass
y0(2)      = ammonia_exp(1);         % Initial Ammonia
y0(3)      = glucose_exp(2);         % Initial Glucose
y0(4)      = fructose_exp(2);        % Initial Fructose
y0(5)      = 0;                      % Initial Ethanol
y0(6)      = PAN_exp(1);             % Initial PAN

% Parameters
index_N     = [bio_DFVB_2, glc_DFVB_2, fru_DFVB_2, n_DFVB_2, etoh_DFVB_2];
alpha_dim   = size(modelo_N_reducido.nullS,2);

% Simulate first period
tic;
[T1, Y1] = ode113(@model_DFVB2fit,t_rango1,theta(1:6),options,theta(7:end), Temperature, modelo_N_reducido, index_N, PAN_struct_2, alpha_dim);
    
% Simulate second period
y0_2     = Y1(end,:);
y0_2(2)  = ammonia_exp(6);
[T2, Y2] = ode113(@model_DFVB2fit,t_rango2,y0_2,options,theta(7:end),Temperature, modelo_N_reducido, index_N, PAN_struct_2, alpha_dim);
fprintf('\nCPU Time(s) = %.2f\n', toc)
% Join the solution
T_sol = [T1(1:end-1); T2];
Y_sol = [Y1(1:end-1, :); Y2];

%% Plot
DFVB_results = table(T_sol, Y_sol(:,1), Y_sol(:,2), ...
    Y_sol(:,3), Y_sol(:,4), Y_sol(:,5), Y_sol(:,6), 'VariableNames',{'Time','Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'} );

Experiment  = table(t_exp, biomass_exp,  ammonia_exp, glucose_exp, fructose_exp, ethanol_exp, PAN_exp,...
    'VariableNames',{'Time','Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'} );

folder = 'Results/Calibration CyT/';

% writetable(DFVB_results, [folder,'Model_calibration_3.csv'])
% writetable(Experiment, [folder,'Experiment_3.csv'])

%%
figure;

% Primer subplot (Biomass)
subplot(3,2,1)
hold on
h1 = plot(T_sol, Y_sol(:,1), Color = color1);
h2 = plot(t_exp, biomass_exp, 'o', Color = color3);
xlim([0 200])
ylim([0, 1.1*max(biomass_exp)])
ylabel('Concentration (g/L)')
title('Biomass')

% Segundo subplot (Glucose)
subplot(3,2,2)
hold on
plot(T_sol, Y_sol(:,3), Color = color1);
plot(t_exp, glucose_exp, 'o', Color = color3);
xlim([0 200])
ylim([0, 1.1*max(Y_sol(:,3))])
title('Glucose')


% Resto de los subplots
subplot(3,2,3)
hold on
plot(T_sol, Y_sol(:,4), Color = color1)
plot(t_exp, fructose_exp, 'o', Color = color3)
xlim([0 200])
ylim([0, 1.1*max(Y_sol(:,4))])
title('Fructose')
ylabel('Concentration (g/L)')

subplot(3,2,4)
hold on
plot(T_sol, Y_sol(:,5), Color = color1)
plot(t_exp, ethanol_exp, 'o', Color = color3)
xlim([0 200])
ylim([0, 1.1*max(ethanol_exp)])
title('Ethanol')

subplot(3,2,5)
hold on
plot(T_sol, Y_sol(:,2), Color = color1)
plot(t_exp, ammonia_exp, 'o', Color = color3)
xlim([0 200])
ylim([0, 1.1*max(ammonia_exp)])
xlabel('Time (h)')
title('Ammonia')
ylabel('Concentration (g/L)')

subplot(3,2,6)
hold on
plot(T_sol, Y_sol(:,6), Color = color1)
plot(t_exp, PAN_exp, 'o', Color = color3)
xlim([0 200])
ylim([0, 1.1*max(PAN_exp)])
xlabel('Time (h)')
title('PAN')

han = legend([h1 h2], {'Model', 'Experimental'});
set(han, 'Position', [0.4 0.01 0.2 0.03], 'Orientation', 'horizontal');


