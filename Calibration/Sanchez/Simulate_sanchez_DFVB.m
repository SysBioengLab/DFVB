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
parameters_all = table2array(readtable("Bootstraping_DFVB_Sanchez.csv"));
parameters = mean(parameters_all, 1);
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

%% Construct DFVB and DFBA models
N_sparse = sparseNull(model.S);
model_DFVB = build_DFVB(model, N_sparse);

%% Index reactions and metabolites

bio_DFVB  = find(startsWith(model_DFVB.rxns, model.rxns(index_rxns(1))));
glc_DFVB  = find(startsWith(model_DFVB.rxns, model.rxns(index_rxns(2))));
o2_DFVB   = find(startsWith(model_DFVB.rxns, model.rxns(index_rxns(end))));
etoh_DFVB = find(startsWith(model_DFVB.rxns, model.rxns(index_rxns(3))));

%% Load Experimental data

exp_data = readtable('Aerobic_data_sanchez_1.csv');

%% DFVB execution: V1

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
index_N     = [bio_DFVB, glc_DFVB, o2_DFVB, etoh_DFVB];

            %   v_glc  K_G  K_iE  f_etoh  f_gly  f_cit  f_lac  
% theta       =  [15     1    1e6    1.7    4    2.5   3];
alpha_dim   = size(model_DFVB.nullS,2);
t_span = [0, 11];
% Simulate
tic;
[T_sol, Y_sol] = ode113(@DFVB2fit_sanchez, t_span,y0,options,theta, model_DFVB, index_N, alpha_dim);
fprintf('\nCPU Time(s) = %.2f\n', toc)

% Plot
figure();
subplot(3,2,1)
hold on
plot(T_sol, Y_sol(:,1))
plot(exp_data.TIME, exp_data.X, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.')
ylabel('Concentration (g/L)')
title('Biomass')

subplot(3,2,2)
hold on
plot(T_sol, Y_sol(:,2))
plot(exp_data.TIME, exp_data.GLC, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.')
title('Glucose')
ylabel('Concentration (g/L)')

subplot(3,2,3)
hold on
plot(T_sol, Y_sol(:,3))
plot(exp_data.TIME, exp_data.ETOH, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.')
title('Ethanol')
ylabel('Concentration (g/L)')
% 
subplot(3,2,4)
hold on
plot(T_sol, Y_sol(:,4))
plot(exp_data.TIME, exp_data.GLY, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.')
title('Glycerol')
ylabel('Concentration (g/L)')

subplot(3,2,5)
hold on
plot(T_sol, Y_sol(:,5))
plot(exp_data.TIME, exp_data.CIT, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.')
title('Citrate')
ylabel('Concentration (g/L)')

subplot(3,2,6)
hold on
plot(T_sol, Y_sol(:,6))
plot(exp_data.TIME, exp_data.LAC, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.')
title('Lactate')
ylabel('Concentration (g/L)')

    %% Recuperar fluxes y alphas no forzado

alpha_N = alpha_recovery_sanchez(T_sol,Y_sol,theta, model_DFVB, index_N, alpha_dim);

% NullS
alpha_no_nulo_N = [];
for i=1:size(alpha_N,2)
    if any(alpha_N(:,i))
        alpha_no_nulo_N = [alpha_no_nulo_N i];   
    end
end

%% Reducir modelo
inactivos_N     = setdiff(1:size(alpha_N,2), alpha_no_nulo_N);
activos_N_con_1 = setdiff(1:size(model_DFVB.nullS,2),inactivos_N); % POS 1 NO SE PUEDE ELIMINAR -> COL 1
activos_N       = setdiff(1:size(model_DFVB.nullS,2),[inactivos_N(2:end)]); % POS 1 NO SE PUEDE ELIMINAR -> COL 1

modelo_N_reducido_con_1   = build_DFVB(model,N_sparse(:,activos_N_con_1));
modelo_N_reducido         = build_DFVB(model,N_sparse(:,activos_N));

%% Análisis primera conversión
exchanges_eliminados = model_DFVB.rxns(find(model_DFVB.nullS(:, inactivos_N(1))~=0));
internas_eliminadas  = model_DFVB.intrxns(find(model_DFVB.intN(:, inactivos_N(1))~=0));
total_rxns = [exchanges_eliminados; internas_eliminadas];

rxn_names_eliminadas = cell(length(total_rxns),1);
for i=1:length(total_rxns)
    idx = find(startsWith(model.rxns, total_rxns(i)));
    rxn_names_eliminadas(i) = model.rxnNames(idx);
end
table(rxn_names_eliminadas)

%% Identificar reacciones y metabolitos eliminados
Inactive_rxn_N        = [modelo_N_reducido.inac_ex; modelo_N_reducido.inac_int];
Idx_inact_rxn_N       = findRxnIDs(model, Inactive_rxn_N);
Idx_act_rxn_N         = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_N));
Idx_inact_mets_N      = met_reduced(model, Idx_act_rxn_N);
fprintf('\nMetabolitos eliminados con DFVB - SparseNull = %d\n', size(Idx_inact_mets_N,1))
fprintf('\nReacciones eliminadas con DFVB - SparseNull = %d\n', size(modelo_N_reducido.inac_ex,1) + size(modelo_N_reducido.inac_int,1))


%% Index reactions and metabolites

bio_DFVB  = find(startsWith(modelo_N_reducido.rxns, model.rxns(index_rxns(1))));
glc_DFVB  = find(startsWith(modelo_N_reducido.rxns, model.rxns(index_rxns(2))));
o2_DFVB   = find(startsWith(modelo_N_reducido.rxns, model.rxns(index_rxns(end))));
etoh_DFVB = find(startsWith(modelo_N_reducido.rxns, model.rxns(index_rxns(3))));

%% DFVB execution: V1

%ODE configuration
options   = odeset('AbsTol',1e-6, 'NonNegative',[1, 2, 3, 4 ,5 ,6]);  % Opciones del ODE

% Initial conditions
% y0(1)      = exp_data.X(1);   % Initial Biomass
% y0(2)      = exp_data.GLC(2);         % Initial Ammonia
% y0(3)      = exp_data.ETOH(1);              % Initial Glucose
% y0(4)      = exp_data.GLY(2);             % Initial Fructose
% y0(5)      = exp_data.CIT(1);                               % Initial Ethanol
% y0(6)      = exp_data.LAC(1);             % Initial PAN

% Parameters
index_N     = [bio_DFVB, glc_DFVB, o2_DFVB, etoh_DFVB];

            %   v_glc  K_G  K_iE  f_etoh  f_gly  f_cit  f_lac  
% theta       =  [15     1    1e6    1.7    4    2.5   3];
alpha_dim   = size(modelo_N_reducido.nullS,2);
t_span = [0, 11];
% Simulate
tic;
[T_sol, Y_sol] = ode113(@DFVB2fit_sanchez, t_span,y0,options,theta, modelo_N_reducido, index_N, alpha_dim);
fprintf('\nCPU Time(s) = %.2f\n', toc)

% Plot
figure();
subplot(3,2,1)
hold on
plot(T_sol, Y_sol(:,1))
plot(exp_data.TIME, exp_data.X, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.', Location='southeast')
ylabel('Concentration (g/L)')
title('Biomass')

subplot(3,2,2)
hold on
plot(T_sol, Y_sol(:,2))
plot(exp_data.TIME, exp_data.GLC, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.', Location='southeast')
title('Glucose')
ylabel('Concentration (g/L)')

subplot(3,2,3)
hold on
plot(T_sol, Y_sol(:,3))
plot(exp_data.TIME, exp_data.ETOH, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.', Location='southeast')
title('Ethanol')
ylabel('Concentration (g/L)')
% 
subplot(3,2,4)
hold on
plot(T_sol, Y_sol(:,4))
plot(exp_data.TIME, exp_data.GLY, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.', Location='southeast')
title('Glycerol')
ylabel('Concentration (g/L)')

subplot(3,2,5)
hold on
plot(T_sol, Y_sol(:,5))
plot(exp_data.TIME, exp_data.CIT, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.', Location='southeast')
title('Citrate')
ylabel('Concentration (g/L)')

subplot(3,2,6)
hold on
plot(T_sol, Y_sol(:,6))
plot(exp_data.TIME, exp_data.LAC, 'ro')
xlim(t_span)
legend('DFVB', 'Exp.', Location='southeast')
title('Lactate')
ylabel('Concentration (g/L)')


%% Exportar modelo para calibrar

modelo = struct();

modelo.model_DFVB  = modelo_N_reducido;
modelo.index       = index_N;
modelo.alpha_dim   = alpha_dim;
modelo.y0          = y0;


%% Export experimental data

biomass_exp  = exp_data.X;
glucose_exp  = exp_data.GLC;
ethanol_exp  = exp_data.ETOH;
glycerol_exp = exp_data.GLY;
citrate_exp  = exp_data.CIT;
lactate_exp  = exp_data.LAC;
t_exp        = exp_data.TIME;

data_exp = struct();
data_exp.t_exp = t_exp;
data_exp.conc  = [biomass_exp, glucose_exp, ethanol_exp, ...
                     glycerol_exp, citrate_exp, lactate_exp];


%% Sensitivity analysis
epsilon   = 1e-3;
n_theta   = length(theta);
n_vars    = 6;    % Biomass, Glucose, Ethanol, Glycerol, Citrate, Lactate
n_time    = length(t_exp);

S_mat     = zeros(n_time * n_vars, n_theta); 
S_rel_avg = zeros(n_theta, n_vars);          

% 1. Nominal simulation evaluated at experimental time points
Y_nom = simulate_reduced_sanchez_model(theta, y0, t_span, t_exp, options, modelo_N_reducido, index_N, alpha_dim);

% Avoid division by zero in relative calculations
Y_nom_safe = Y_nom;
Y_nom_safe(Y_nom_safe == 0) = 1e-8;

for j = 1:n_theta
    theta_plus  = theta;
    theta_minus = theta;
    
    delta = epsilon * abs(theta(j));
    if delta == 0
        delta = 1e-6; 
    end
    
    % Theoretical perturbations
    theta_plus(j)  = theta(j) + delta;
    theta_minus(j) = theta(j) - delta;
    
    % Prevent non-negativity constraint violation
    if theta_minus(j) < 0
        theta_minus(j) = 0; 
    end
    
    % Recalculate effective deltas (handles truncation at 0)
    delta_plus  = theta_plus(j) - theta(j);
    delta_minus = theta(j) - theta_minus(j);
    
    % Simulate perturbed models
    Y_plus  = simulate_reduced_sanchez_model(theta_plus, y0, t_span, t_exp, options, modelo_N_reducido, index_N, alpha_dim);
    Y_minus = simulate_reduced_sanchez_model(theta_minus, y0, t_span, t_exp, options, modelo_N_reducido, index_N, alpha_dim);
    
    % Absolute sensitivity (adjusted for asymmetric deltas)
    dY_dtheta = (Y_plus - Y_minus) / (delta_plus + delta_minus);
    
    % Relative sensitivity
    S_rel = dY_dtheta .* (theta(j) ./ Y_nom_safe);
    
    % Store metrics
    S_rel_avg(j, :) = mean(abs(S_rel), 1, 'omitnan');
    S_mat(:, j) = dY_dtheta(:);
end

% 3. Correlation Matrix via FIM
FIM = S_mat' * S_mat;
Cov = pinv(FIM); % Pseudoinverse due to potential singularity in metabolic models

Corr_matrix = zeros(n_theta, n_theta);
for i = 1:n_theta
    for j = 1:n_theta
        Corr_matrix(i,j) = Cov(i,j) / sqrt(Cov(i,i) * Cov(j,j));
    end
end

%% Visualization
% Full Parameter Correlation Heatmap
figure;
imagesc(abs(Corr_matrix));
colormap(whiteredcmap);
colorbar;
clim([0 1]);
title('Parameter Correlation Matrix DFVB (|\rho|)');
xlabel('Parameter Index');
ylabel('Parameter Index');

% Parameter Correlation Heatmap (Threshold >= 0.95)
figure;
imagesc(abs(Corr_matrix));
m = 1000; % High resolution for a precise cutoff
cmap = ones(m, 3); % Initialize all as white [1 1 1]
threshold_idx = round(0.95 * m); % Calculate index for 0.95
cmap(threshold_idx:end, :) = repmat([1 0 0], m - threshold_idx + 1, 1); % Red [1 0 0] for >= 0.95
colormap(cmap);
colorbar;
clim([0 1]);
title('Parameter Correlation Matrix DFVB (|\rho| \geq 0.95)');
xlabel('Parameter Index');
ylabel('Parameter Index');

% Average Relative Sensitivity Bar Chart
figure;
bar_labels = {'Biomass', 'Glucose', 'Ethanol', 'Glycerol', 'Citrate', 'Lactate'};
b = bar(S_rel_avg', 'grouped');
set(gca, 'XTickLabel', bar_labels);
title('Average Relative Sensitivity by Variable');
ylabel('Sensitivity |(\partial Y / \partial \theta) \cdot (\theta / Y)|');
legend_labels = arrayfun(@(x) sprintf('\\theta_{%d}', x), 1:n_theta, 'UniformOutput', false);
legend(b, legend_labels, 'Location', 'northeastoutside');

disp('Sensitivity and correlation analysis completed.');
writetable(table(Corr_matrix), "Calibration/Results_Bootstraping/Correlation_matrix_DFVB_Sanchez.csv")

% ----------------------------------------------------------------------- %
% HELPER FUNCTIONS
% ----------------------------------------------------------------------- %

function Y_interp = simulate_reduced_sanchez_model(p_kin, y_init, t_range, t_points, opt, model, idx_N, alpha_d)
    % Simulates a single period and interpolates to experimental time points
    [T_sim, Y_sim] = ode113(@DFVB2fit_sanchez, t_range, y_init, opt, p_kin, model, idx_N, alpha_d);
    Y_interp = interp1(T_sim, Y_sim, t_points, 'linear', 'extrap');
end

function cmap = whiteredcmap(m)
    % Sequential White-Red colormap generator
    if nargin < 1
        m = 64; 
    end
    c1 = [1 1 1]; % White (0)
    c2 = [1 0 0]; % Red (1)
    
    % Linear interpolation for each RGB channel
    cmap = [linspace(c1(1), c2(1), m)', ...
            linspace(c1(2), c2(2), m)', ...
            linspace(c1(3), c2(3), m)'];
end