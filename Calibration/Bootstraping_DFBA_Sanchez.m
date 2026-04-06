// addpath('/home/GUROBI.LIC/gurobi1300/linux64/matlab');
// gurobi_setup;

// currentScriptPath = pwd;
// rootPathDceba = fileparts(currentScriptPath); 
// addpath(genpath(rootPathDceba));

numWorkers = 8;
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

%% Sanchez data

load data_exp_R6.mat
load modelo_DFBA.mat
modelo = modelo_dfba_struct;

theta = [14  0.01  1000   0.5   5   2  4.5];
theta_init = [modelo.y0, theta];

% ----------------------- User parameters ------------------------- %
Nboots = 1000;            % number of bootstrap replicates
includeFirst = true;     % force inclusion of the first time point in each replicate
rngBase = 1234;          % base seed for reproducibility

% --------------------------- Original bounds -------------------------- %
% Original parameters
C_LB = modelo.y0*0.7;
C_UB = modelo.y0*1.3;

p_LB = [1   1e-4  1e-3   0.1    0.01   0.01   0.01  ];
p_UB = [50  1e4   1e5    100    10     10     10    ];

% ------------------- Optimizer options (local) ------------------ %
opts = optimoptions('patternsearch', ...
    'SearchFcn', [], ...                   % Disable GA-based search
    'PollMethod', 'GSSPositiveBasis2N', ...% Use GSS for a more robust exploration in polling
    'MeshRotate', 'on', ...                % Rotation to avoid aligned local minima
    'InitialMeshSize', 10, ...             % Large initial steps
    'ScaleMesh', true, ...                 % Adaptive scaling of variables
    'FunctionTolerance', 1e-6, ...         % Stopping criterion by tolerance
    'StepTolerance', 1e-6, ...
    'MeshTolerance', 1e-6, ...
    'UseParallel', false, ...
    'MaxTime', 1800, ...                   % Strict limit of 20 minutes (1200 seconds)
    'MaxFunctionEvaluations', 10000);   % Adjust this number: maximum evaluations before MaxTime is reached

% ---------------------- Preparation and preallocation --------------------- %
rng(rngBase, 'twister');
seeds = randi(1e9, Nboots, 1);

npar   = numel([C_LB, p_LB]);
Xopt   = NaN(Nboots, npar);
Fval   = NaN(Nboots, 1);
Exit   = NaN(Nboots, 1);
Timing = NaN(Nboots, 1);
Theta0 = NaN(Nboots, npar);

% ======================== Bootstrap loop (parallel) ======================== %
tic;
parfor b = 1:Nboots
    % Independent seed per replicate
    s = RandStream('Threefry','Seed',seeds(b));
    RandStream.setGlobalStream(s);

    % --- Aligned resampling with replacement (same rows for all columns) ---
    [data_b, idx_b] = resampleAlignedWithReplacement(data_exp, includeFirst);

    % Objective function
    fun_b = @(x) LossFunction_DFBA_sanchez(x, modelo, data_b);

    % Optimization (sequential within each replicate)
    [xopt_b, fval_b, exit_b] = patternsearch( ...
        fun_b, theta_init, [], [], [], [], [C_LB, p_LB], [C_UB, p_UB], [], opts);

    % Save results
    Xopt(b,:)   = xopt_b;
    Fval(b)     = fval_b;
    Exit(b)     = exit_b;
    Theta0(b,:) = theta_init;
    % fprintf('Iteration %d\n', b)
end
final_time = toc;

% fprintf('CPU time = %.4f \n', final_time)

%%
% ------------------------- Consolidation/saving ----------------------- %
varsParam = "p" + (1:npar);
Tresult = table((1:Nboots).', Exit, Fval, ...
    'VariableNames', {'boot','exitflag','fval'});
Ttheta0 = array2table(Theta0, 'VariableNames', cellstr("theta0_"+varsParam));
Tparams = array2table(Xopt,  'VariableNames', cellstr(varsParam));
ResultadoBootstrap = [Tresult Ttheta0 Tparams];

save('Results_Bootstraping/Bootstraping_1000_DFBA_Sanchez.mat', 'Xopt','Fval','Exit','Timing','Theta0','p_LB','p_UB');

fprintf("Saved Sanchez!\n")
