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
load('data_exp_R2.mat');
theta_all = table2array(readtable('Bootstraping_DFBA_CYT.csv'));
theta = mean(theta_all, 1);

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

%% Transform model to gurobi

model_DFBA = model2gurobi(model_aa_irrev);
pan_dfba  = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(PAN)));

%% Add restricion

A_row = zeros(1, size(model_DFBA.A,2));
A_row(1, pan_dfba(end-5:end)) = ones(1,6);
b_row                         = 0;
sense_new                     = '<'; 

model_DFBA.A     = [model_DFBA.A; A_row];
model_DFBA.rhs   = [model_DFBA.rhs; b_row];
model_DFBA.sense = [model_DFBA.sense; sense_new];


%% Index reactions and metabolites

bio_dfba  = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(index_rxns(1))));
glc_dfba  = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(index_rxns(2))));
fru_dfba  = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(index_rxns(4))));
n_dfba    = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(index_rxns(5))));
etoh_dfba = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(index_rxns(6))));
pan_dfba  = find(startsWith(model_DFBA.rxns, model_aa_irrev.rxns(PAN)));

PAN_struct = struct('index_dfna', pan_dfba(end-5:end), 'PAN_prop', PAN_prop, 'MW_PAN', MW_PAN);


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


%% DFNA execution

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
index_N     = [bio_dfba, glc_dfba, fru_dfba, n_dfba, etoh_dfba];

% Simulate first period
tic;
[T1, Y1] = ode113(@model_dFBA2fit,t_rango1,theta(1:6),options,theta(7:end), Temperature, model_DFBA, index_N, PAN_struct);

% Simulate second period
y0_2     = Y1(end,:);
y0_2(2)  = ammonia_exp(6);
[T2, Y2] = ode113(@model_dFBA2fit,t_rango2,y0_2,options,theta(7:end),Temperature, model_DFBA, index_N, PAN_struct);
fprintf('\nCPU Time(s) = %.2f\n', toc)

% Join the solution
T_sol = [T1; T2];
Y_sol = [Y1; Y2];

% Plot

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


 %% Exportar modelo para calibrar
 
modelo_struct = struct();

modelo_struct.model_DFBA  = model_DFBA;
modelo_struct.index       = index_N;
modelo_struct.PAN_struct  = PAN_struct;
modelo_struct.Temperature = Temperature;
modelo_struct.y0          = y0;
modelo_struct.y0_2        = y0_2;
modelo_struct.t_rango1    = t_rango1;
modelo_struct.t_rango2    = t_rango2;

%% Export results
DFBA_results = table(T_sol, Y_sol(:,1), Y_sol(:,2), ...
    Y_sol(:,3), Y_sol(:,4), Y_sol(:,5), Y_sol(:,6), 'VariableNames',{'Time','Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'} );

% Experiment  = table(t_exp, biomass_exp,  ammonia_exp, glucose_exp, fructose_exp, ethanol_exp, PAN_exp,...
%     'VariableNames',{'Time','Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'} );

folder = 'Results/Calibration CyT/';

% writetable(DFBA_results, [folder,'Model_calibration_DFBA_3.csv'])
% writetable(Experiment, [folder,'Experiment_3.csv'])




%% Sensitivity and Correlation Analysis
epsilon   = 1e-3;
n_theta   = length(theta);
n_vars    = 6;    % Biomass, Ammonia, Glucose, Fructose, Ethanol, PAN
n_time    = length(t_exp);

S_mat     = zeros(n_time * n_vars, n_theta); 
S_rel_avg = zeros(n_theta, n_vars);          

Y_nom = simulate_dfba_model(theta, t_rango1, t_rango2, t_exp, ammonia_exp(6), options, Temperature, model_DFBA, index_N, PAN_struct);

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
    Y_plus  = simulate_dfba_model(theta_plus, t_rango1, t_rango2, t_exp, ammonia_exp(6), options, Temperature, model_DFBA, index_N, PAN_struct);
    Y_minus = simulate_dfba_model(theta_minus, t_rango1, t_rango2, t_exp, ammonia_exp(6), options, Temperature, model_DFBA, index_N, PAN_struct);
    
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
% Parameter Correlation Heatmap
figure;
imagesc(abs(Corr_matrix));
colormap(whiteredcmap);
colorbar;
clim([0 1]);
title('Parameter Correlation Matrix DFBA (\theta)');
xlabel('Parameter Index');
ylabel('Parameter Index');

% Parameter Correlation Heatmap
figure;
imagesc(abs(Corr_matrix));

% Create a custom colormap with a hard threshold at 0.95
m = 1000; % High resolution for a precise cutoff
cmap = ones(m, 3); % Initialize all as white [1 1 1]
threshold_idx = round(0.95 * m); % Calculate index for 0.95
cmap(threshold_idx:end, :) = repmat([1 0 0], m - threshold_idx + 1, 1); % Red [1 0 0] for >= 0.95

colormap(cmap);
colorbar;
clim([0 1]);
title('Parameter Correlation Matrix DFBA (|\rho| \geq 0.95)');
xlabel('Parameter Index');
ylabel('Parameter Index');

% Average Relative Sensitivity Bar Chart
figure;
bar_labels = {'Biomass', 'Ammonia', 'Glucose', 'Fructose', 'Ethanol', 'PAN'};
b = bar(S_rel_avg', 'grouped');
set(gca, 'XTickLabel', bar_labels);
title('Average Relative Sensitivity by Variable');
ylabel('Sensitivity |(\partial Y / \partial \theta) \cdot (\theta / Y)|');
legend_labels = arrayfun(@(x) sprintf('\\theta_{%d}', x), 1:n_theta, 'UniformOutput', false);
legend(b, legend_labels, 'Location', 'northeastoutside');

disp('Sensitivity and correlation analysis completed.');
writetable(table(Corr_matrix(7:end, 7:end)), "Calibration/Results_Bootstraping/Correlation_matrix_DFBA_CYT.csv")

% ----------------------------------------------------------------------- %
% HELPER FUNCTIONS
% ----------------------------------------------------------------------- %

function Y_interp = simulate_dfba_model(p, t1, t2, t_exp, NH3_mid, opt, Temp, model, idx_N, PAN_s)
    % Simulates both periods and interpolates to experimental time points.
    
    % p(1:6): initial conditions. p(7:end): kinetic parameters.
    [T1_sim, Y1_sim] = ode113(@model_dFBA2fit, t1, p(1:6), opt, p(7:end), Temp, model, idx_N, PAN_s);
    
    y0_2     = Y1_sim(end,:);
    y0_2(2)  = NH3_mid; % Experimental correction between periods
    [T2_sim, Y2_sim] = ode113(@model_dFBA2fit, t2, y0_2, opt, p(7:end), Temp, model, idx_N, PAN_s);
    
    % Merge solutions
    T_sol_full = [T1_sim(1:end-1); T2_sim];
    Y_sol_full = [Y1_sim(1:end-1, :); Y2_sim];
    
    % Interpolate to experimental time points
    Y_interp = interp1(T_sol_full, Y_sol_full, t_exp, 'linear', 'extrap');
end

function cmap = whiteredcmap(m)
    % Generador de colormap secuencial Blanco-Rojo
    if nargin < 1
        m = 64; % Resolución por defecto
    end
    
    c1 = [1 1 1]; % Blanco (0)
    c2 = [1 0 0]; % Rojo (1)
    
    % Interpolación lineal para cada canal RGB
    cmap = [linspace(c1(1), c2(1), m)', ...
            linspace(c1(2), c2(2), m)', ...
            linspace(c1(3), c2(3), m)'];
end