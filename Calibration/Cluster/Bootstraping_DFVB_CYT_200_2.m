addpath('/home/GUROBI.LIC/gurobi1300/linux64/matlab');
gurobi_setup;

addpath(genpath(pwd));

%% CyT datos
load data_exp_R1.mat
load model_DFVB_R1.mat

modelo     = modelo_struct; 
theta_init = [0.3, 0.5, 100, 80, 0, 0.2, 0.7, 0.15, 150, 200, 1000, 50, 10000, 15, 0.5, 0.9, 0.8];

% ----------------------- Parámetros del usuario ------------------------- %
Nboots = 200;            % número de réplicas bootstrap
includeFirst = true;      % forzar inclusión del primer tiempo en cada réplica
rngBase = 1232;           % semilla base para reproducibilidad

% --------------------------- Cotas originales -------------------------- %
C_LB = [0.212	0.047	92.086	58.465	0.000	0.049];
C_UB = [0.419	0.052	121.601	92.555	0.000	0.358];

p_LB = [1e-5   1e-8    1e-2  1e-2  1e-3   1e-2   1e-2   1e-5 1e-5 1e-5 1e-5];
p_UB = [2.0      1     1e3   1e3   1e4    5e3    1e4    1e4  5    2    2];

% ------------------- Opciones del optimizador (local) ------------------ %
opts = optimset('Display', 'off', ...                 
    'MaxIter', 150, ...           
    'TolFun', 1e-4); 
% ---------------------- Preparación y prealocación --------------------- %
rng(rngBase, 'twister');
seeds = randi(1e9, Nboots, 1);

npar   = numel([C_LB, p_LB]);
Xopt   = NaN(Nboots, npar);
Fval   = NaN(Nboots, 1);
Exit   = NaN(Nboots, 1);
Timing = NaN(Nboots, 1);
Theta0 = NaN(Nboots, npar);


% ======================== Bucle bootstrap (par) ======================== %
for b = 1:Nboots
    % Semilla independiente por réplica
    s = RandStream('Threefry','Seed',seeds(b));
    RandStream.setGlobalStream(s);

    % --- Muestreo alineado con reemplazo (mismas filas para todas las columnas) ---
    [data_b, idx_b] = resampleAlignedWithReplacement(data_exp, includeFirst);

    % Función objetivo
    fun_b = @(x) LossFunction(x, modelo, data_b);

    % Optimización (secuencial dentro de cada réplica)
    t0 = tic;
    [xopt_b, fval_b, exit_b, output] = fminsearchbnd(fun_b, theta_init, [C_LB, p_LB], [C_UB, p_UB], opts);
    Timing(b) = toc(t0);

    % Guardado de resultados
    Xopt(b,:)   = xopt_b;
    Fval(b)     = fval_b;
    Exit(b)     = exit_b;
    Theta0(b,:) = theta_init;
end

%%
% ------------------------- Consolidación/guardado ----------------------- %
varsParam = "p" + (1:npar);
Tresult = table((1:Nboots).', Exit, Fval, Timing, 'VariableNames', {'boot','exitflag','fval','time_s'});
Ttheta0 = array2table(Theta0, 'VariableNames', cellstr("theta0_"+varsParam));
Tparams = array2table(Xopt,  'VariableNames', cellstr(varsParam));
ResultadoBootstrap = [Tresult Ttheta0 Tparams];

save('Results_Bootstraping/Bootstraping_200_DFVB_CYT_2.mat', 'Xopt','Fval','Exit','Timing','Theta0','p_LB','p_UB');

fprintf("Guardado CYT!\n")
