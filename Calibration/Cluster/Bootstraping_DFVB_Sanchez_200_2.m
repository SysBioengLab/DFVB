addpath('/home/GUROBI.LIC/gurobi1300/linux64/matlab');
gurobi_setup;

addpath(genpath(pwd));

%% Sanchez datos

load data_exp_R6.mat
load model_DFVB_R6.mat
%%
theta = [14  0.01  100   0.5   5   2  4.5];
theta_init = [modelo.y0, theta];

% ----------------------- Parámetros del usuario ------------------------- %
Nboots = 200;              % número de réplicas bootstrap
includeFirst = true;     % forzar inclusión del primer tiempo en cada réplica
rngBase = 1232;          % semilla base para reproducibilidad

% --------------------------- Cotas originales -------------------------- %
% Parámetros originales
C_LB = modelo.y0*0.7;
C_UB = modelo.y0*1.3;

p_LB = [1   1e-4  1e-3   0.1    0.01   0.01   0.01  ];
p_UB = [50  1e4   2e4    100    10     10     10    ];

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
tic;
for b = 1:Nboots
    % Semilla independiente por réplica
    s = RandStream('Threefry','Seed',seeds(b));
    RandStream.setGlobalStream(s);

    % --- Muestreo alineado con reemplazo (mismas filas para todas las columnas) ---
    [data_b, idx_b] = resampleAlignedWithReplacement(data_exp, includeFirst);

    % Función objetivo
    fun_b = @(x) LossFunction_sanchez(x, modelo, data_b);

    % Optimización (secuencial dentro de cada réplica)
    [xopt_b, fval_b, exit_b, output] = fminsearchbnd(fun_b, theta_init, [C_LB, p_LB], [C_UB, p_UB], opts);

    % Guardado de resultados
    Xopt(b,:)   = xopt_b;
    Fval(b)     = fval_b;
    Exit(b)     = exit_b;
    Theta0(b,:) = theta_init;
    fprintf("Iteración %d/%d completada\n", b, Nboots);
end
final_time = toc;

%%
% ------------------------- Consolidación/guardado ----------------------- %
varsParam = "p" + (1:npar);
Tresult = table((1:Nboots).', Exit, Fval, ...
    'VariableNames', {'boot','exitflag','fval'});
Ttheta0 = array2table(Theta0, 'VariableNames', cellstr("theta0_"+varsParam));
Tparams = array2table(Xopt,  'VariableNames', cellstr(varsParam));
ResultadoBootstrap = [Tresult Ttheta0 Tparams];

save('Results_Bootstraping/Bootstraping_200_2_DFVB_Sanchez.mat', 'Xopt','Fval','Exit','Timing','Theta0','p_LB','p_UB');

fprintf("Guardado Sanchez!\n")
