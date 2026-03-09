// addpath('/home/GUROBI.LIC/gurobi1300/linux64/matlab');
// gurobi_setup;

// currentScriptPath = pwd;
// rootPathDceba = fileparts(currentScriptPath); 
// addpath(genpath(rootPathDceba));

numWorkers = 8;
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

%% Sanchez datos

load data_exp_R6.mat
load modelo_DFBA.mat
modelo = modelo_dfba_struct;

theta = [14  0.01  1000   0.5   5   2  4.5];
theta_init = [modelo.y0, theta];

% ----------------------- Parámetros del usuario ------------------------- %
Nboots = 1000;            % número de réplicas bootstrap
includeFirst = true;     % forzar inclusión del primer tiempo en cada réplica
rngBase = 1234;          % semilla base para reproducibilidad

% --------------------------- Cotas originales -------------------------- %
% Parámetros originales
C_LB = modelo.y0*0.7;
C_UB = modelo.y0*1.3;

p_LB = [1   1e-4  1e-3   0.1    0.01   0.01   0.01  ];
p_UB = [50  1e4   1e5    100    10     10     10    ];

% ------------------- Opciones del optimizador (local) ------------------ %
opts = optimoptions('patternsearch', ...
    'SearchFcn', [], ...                   % Desactivamos la búsqueda de GA
    'PollMethod', 'GSSPositiveBasis2N', ...% Usamos GSS para una exploración más robusta en el sondeo.
    'MeshRotate', 'on', ...                % Rotación para evitar mínimos locales alineados.
    'InitialMeshSize', 10, ...             % Pasos iniciales grandes.
    'ScaleMesh', true, ...                 % Escalado adaptativo de variables.
    'FunctionTolerance', 1e-6, ...         % Criterio de parada por tolerancia.
    'StepTolerance', 1e-6, ...
    'MeshTolerance', 1e-6, ...
    'UseParallel', false, ...
    'MaxTime', 1800, ...                   % Límite estricto de 20 minutos (1200 segundos).
    'MaxFunctionEvaluations', 10000);   % Ajusta este número: es el máximo de evaluaciones antes de que salte el MaxTime.;

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
parfor b = 1:Nboots
    % Semilla independiente por réplica
    s = RandStream('Threefry','Seed',seeds(b));
    RandStream.setGlobalStream(s);

    % --- Muestreo alineado con reemplazo (mismas filas para todas las columnas) ---
    [data_b, idx_b] = resampleAlignedWithReplacement(data_exp, includeFirst);

    % Función objetivo
    fun_b = @(x) LossFunction_DFBA_sanchez(x, modelo, data_b);

    % Optimización (secuencial dentro de cada réplica)
    [xopt_b, fval_b, exit_b] = patternsearch( ...
        fun_b, theta_init, [], [], [], [], [C_LB, p_LB], [C_UB, p_UB], [], opts);

    % Guardado de resultados
    Xopt(b,:)   = xopt_b;
    Fval(b)     = fval_b;
    Exit(b)     = exit_b;
    Theta0(b,:) = theta_init;
    % fprintf('Iteracion %d\n', b)
end
final_time = toc;

% fprintf('CPU time = %.4f \n', final_time)

%%
% ------------------------- Consolidación/guardado ----------------------- %
varsParam = "p" + (1:npar);
Tresult = table((1:Nboots).', Exit, Fval, ...
    'VariableNames', {'boot','exitflag','fval'});
Ttheta0 = array2table(Theta0, 'VariableNames', cellstr("theta0_"+varsParam));
Tparams = array2table(Xopt,  'VariableNames', cellstr(varsParam));
ResultadoBootstrap = [Tresult Ttheta0 Tparams];

save('Results_Bootstraping/Bootstraping_1000_DFBA_Sanchez.mat', 'Xopt','Fval','Exit','Timing','Theta0','p_LB','p_UB');

fprintf("Guardado Sanchez!\n")






