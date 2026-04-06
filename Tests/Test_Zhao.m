% Plot configuration
set(groot, 'DefaultAxesFontSize', 14); % Font size
set(groot, 'DefaultAxesLineWidth', 1.2); % Axes line width
set(groot, 'DefaultLineLineWidth', 1.4); % Line width
set(groot, 'DefaultFigureColor', 'w'); % White figure background
set(groot, 'DefaultAxesFontName', 'Arial'); % Default font
set(groot, 'DefaultAxesBox', 'on');
set(groot, 'DefaultAxesXColor', 'k');
set(groot, 'DefaultAxesYColor', 'k');
set(groot, 'DefaultFigurePosition', [100, 100, 800, 600]);
output_folder = 'Figures';

color1 = '#f19eec';
color2 = '#8ee1c2';
color3 = '#9190e1';

changeCobraSolverParams('LP','optTol',1e-8);
changeCobraSolverParams('LP','feasTol',1e-8);
%%

load("zhao_toymodel.mat")
N = sparseNull(model.S);

model_DFVB = build_DFVB_example(model, N);

%% DFVB execution

% ODE configuration
t         = linspace(0,1,1000);                  % Time interval
y0        = zeros(3,1);              % Initial conditions
options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 1;                      % Initial biomass condition
y0(2)      = 20;                     % Initial substrate condition
y0(3)      = 0 ;                     % Initial product condition
params     = [3.8 1];               % Glucose consumption parameters (V_max, K_S)
alpha_dim_N = size(model_DFVB.nullS,2);
index       = [2, 1, 3];

% NullS execution
tic; % Start timer
[T1,Y1] = ode113(@DFVB_example,t,y0,options,model_DFVB,index, params,alpha_dim_N);
tiempo_transcurrido_DFVB = toc; % Save total execution time
disp(['DFVB execution time: ' num2str(tiempo_transcurrido_DFVB) ' seconds']);

%% DFVB average time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
    % sparseNull
    tic;
    [~,~] = ode113(@DFVB_example,t,y0,options,model_DFVB,index, params,alpha_dim_N);
    t_N_sparse(i) = toc;
end

fprintf('\nMean CPU time of sparseNull = %.4f s\n',mean(t_N_sparse))


%% DFBA execution

% ODE configuration
t         = linspace(0,1,1000);                  % Time interval
y0        = zeros(3,1);              % Initial conditions
options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 1;                      % Initial biomass condition
y0(2)      = 20;                     % Initial substrate condition
y0(3)      = 0 ;                     % Initial product condition
params     = [3.8 1];               % Glucose consumption parameters (V_max, K_S)

index_N     = [6, 1, 7];

% DFBA execution
tic; % Start timer
[T2,Y2] = ode113(@DFBA_example,t,y0,options,model,index_N,params);
tiempo_transcurrido_DFBA = toc; % Save total execution time
disp(['DFBA execution time: ' num2str(tiempo_transcurrido_DFBA) ' seconds']);

%% DFBA average time calculation
n = 5;
t_DFBA   = zeros(n,1);

for i=1:n
    % sparseNull
    tic;
    [~,~] = ode113(@DFBA_example,t,y0,options,model,index_N,params);
    t_DFBA(i) = toc;
end

fprintf('\nMean CPU time of DFBA = %.4f s\n',mean(t_DFBA))


%%
figure;
hold on
plot(T1, Y1(:,1),  'Color', color1, 'LineWidth', 3)
plot(T1, Y1(:,2),  'Color', color2, 'LineWidth', 3)
plot(T1, Y1(:,3),  'Color', color3, 'LineWidth', 3)
legend({'X', 'S', 'P'})
xlabel('Time (h)')
ylabel('Concentration (mM)')
% ylim([-0.5, 30])
title('DFVB Simulation')

% exportgraphics(gcf, strcat(output_folder,'/DFVB_results.png'), 'Resolution',300 )


%%
figure;
subplot(1,2,1)
plot(T1, Y1, LineWidth=1.2)
legend({'c_x(t)', 'c_S(t)', 'c_P(t)'})
grid on
ylim([-5, 30])
title('DFBA simulation')

subplot(1,2,2)
plot(T2, Y2, LineWidth=1.2)
legend({'c_x(t)', 'c_S(t)', 'c_P(t)'})
grid on
ylim([-5, 30])
title('DFVB simulation')

figure;
hold on
plot(T1, Y1(:,1), LineWidth=1.2)
plot(T1, Y1(:,1), LineWidth=1.2)

plot(T1, Y1(:,2), LineWidth=1.2)
plot(T1, Y2(:,2), LineWidth=1.2)

plot(T1, Y1(:,3), LineWidth=1.2)
plot(T1, Y2(:,3), LineWidth=1.2)
ylim([-5, 30])
ylabel('Concentration (mM)')
xlabel('Time (h)')
legend({'X, DFVB', 'S, DFVB', 'P, DFVB', 'X, DFBA', 'S, DFBA', 'P, DFBA'}, 'NumColumns',2)

%% Calculate RMSE for biomass, substrate and product

bio_DFVB   = Y1(:,1);
bio_DFBA   = Y2(:,1);

S_DFVB    = Y1(:,2);
S_DFBA    = Y2(:,2);

P_DFVB    = Y1(:,3);
P_DFBA    = Y2(:,3);
% Interpolation of DFVB at DFBA times
bio_DFVB_interp = interp1(T1, bio_DFVB, T2, 'linear', 'extrap');
S_DFVB_interp   = interp1(T1, S_DFVB, T2, 'linear', 'extrap');
P_DFVB_interp   = interp1(T1, P_DFVB, T2, 'linear', 'extrap');

% RMSE calculation
mse_bio         = sqrt(mean((bio_DFBA - bio_DFVB_interp).^2));
mse_S           = sqrt(mean((S_DFBA - S_DFVB_interp).^2));
mse_P           = sqrt(mean((P_DFBA - P_DFVB_interp).^2));

% Display results
disp('RMSE Biomass')
disp(['The RMSE between DFVB and DFBA is: ', num2str(mse_bio)]);
disp('RMSE Substrate')
disp(['The RMSE between DFVB and DFBA is: ', num2str(mse_S)]);
disp('RMSE Product')
disp(['The RMSE between DFVB and DFBA is: ', num2str(mse_P)]);
%% Alphas recovery (conversions)

alphas   = zeros(length(T1),alpha_dim_N);
fluxes   = zeros(length(T1),size(model_DFVB.lb,1)-alpha_dim_N);

for i=1:length(T1)
    [~, a, v]    = DFVB_example(T1(i),Y1(i,:),model_DFVB,index, params,alpha_dim_N);
    alphas(i,:)  = a;
    fluxes(i,:)  = v;
end


figure;
hold on
plot(T1, alphas(:,2), 'Color', color2, 'LineWidth', 3)
plot(T1, alphas(:,1), 'Color', color1, 'LineWidth', 3)
legend({'\alpha_1', '\alpha_2'}, 'Location', 'northeast')
% grid on
title('Weight dynamics')
ylabel('Activity')
xlabel('Time (h)')
% exportgraphics(gcf, strcat(output_folder,'/Alpha_Zhao.png'), 'Resolution',300 )
%% Inactive reactions
modelo_reducido       = build_DFVB_example(model,N(:, 1));

Inactive_rxn_N        = [modelo_reducido.inac_ex; modelo_reducido.inac_int];

Idx_act_rxn_N        = setdiff(1:size(model.S,2) ,findRxnIDs(model, Inactive_rxn_N));

Idx_inact_mets_N     = met_reduced(model, Idx_act_rxn_N);

fprintf('\nMetabolites removed with DFVB = %d\n', size(Idx_inact_mets_N,1))
fprintf('\nReactions removed with DFVB = %d\n', size(modelo_reducido.inac_ex,1) + size(modelo_reducido.inac_int,1))

%% DFVB execution (reduced model)

% ODE configuration
t         = [0,1];                  % Time interval
y0        = zeros(3,1);              % Initial conditions
options   = odeset('AbsTol',1e-6,'NonNegative',1:length(y0));  % ODE options

% Initial conditions and parameters
y0(1)      = 1;                      % Initial biomass condition
y0(2)      = 20;                     % Initial substrate condition
y0(3)      = 0 ;                     % Initial product condition
params     = [3.8 1];               % Glucose consumption parameters (V_max, K_S)
alpha_dim_N = size(model_DFVB.nullS,2);
index_reducido = [2, 1, 1];


% NullS execution
tic; % Start timer
[T_red,Y_red] = ode113(@DFVB_example,t,y0,options,modelo_reducido,index_reducido,params, alpha_dim_N);
tiempo_transcurrido_DFVB = toc; % Save total execution time
disp(['DFVB execution time: ' num2str(tiempo_transcurrido_DFVB) ' seconds']);

%% DFVB average time calculation
n = 5;
t_N_sparse   = zeros(n,1);

for i=1:n
    % sparseNull
    tic;
    [~,~] = ode113(@DFVB_example,t,y0,options,modelo_reducido,index_reducido,params, alpha_dim_N);
    t_N_sparse(i) = toc;
end

fprintf('\nMean CPU time of sparseNull = %.4f s\n',mean(t_N_sparse))


%% Calculate MSE for biomass, substrate and product (reduced model)

bio_DFVB   = Y_red(:,1);
S_DFVB    = Y_red(:,2);
P_DFVB    = Y_red(:,3);

% Interpolation of DFVB (reduced) at DFBA times
bio_DFVB_interp = interp1(T_red, bio_DFVB, T2, 'linear', 'extrap');
S_DFVB_interp   = interp1(T_red, S_DFVB, T2, 'linear', 'extrap');
P_DFVB_interp   = interp1(T_red, P_DFVB, T2, 'linear', 'extrap');

% MSE calculation
mse_bio         = mean((bio_DFBA - bio_DFVB_interp).^2);
mse_S           = mean((S_DFBA - S_DFVB_interp).^2);
mse_P           = mean((P_DFBA - P_DFVB_interp).^2);

% Display result
disp('MSE Biomass')
disp(['The MSE between DFVB and DFBA is: ', num2str(mse_bio)]);
disp('MSE Substrate')
disp(['The MSE between DFVB and DFBA is: ', num2str(mse_S)]);
disp('MSE Product')
disp(['The MSE between DFVB and DFBA is: ', num2str(mse_P)]);
%%
figure;
subplot(1,2,1)
plot(T_red, Y_red, LineWidth=1.2)
legend({'c_x(t)', 'c_S(t)', 'c_P(t)'})
grid on
ylim([-5, 30])
title('DFVB reduced simulation')

subplot(1,2,2)
plot(T2, Y2, LineWidth=1.2)
legend({'c_x(t)', 'c_S(t)', 'c_P(t)'})
grid on
ylim([-5, 30])
title('DFBA simulation')

%% Export datos
T_DFVB = table(T1, Y1(:,1), Y1(:,2), Y1(:,3), alphas(:,1), alphas(:,2), ...
    VariableNames={'Time', 'X', 'S', 'P', 'Alpha1', 'Alpha2'});

T_DFBA = table(T2, Y2(:,1), Y2(:,2), Y2(:,3), ...
    VariableNames={'Time', 'X', 'S', 'P',});

T_DFVB_red = table(T_red, Y_red(:,1), Y_red(:,2), Y_red(:,3), ...
    VariableNames={'Time', 'X', 'S', 'P',});

% 
writetable(T_DFVB, 'Results/toymodel/DFVB_results.csv')
writetable(T_DFBA, 'Results/toymodel/DFBA_results.csv')
writetable(T_DFVB_red, 'Results/toymodel/DFVB_red_results.csv')