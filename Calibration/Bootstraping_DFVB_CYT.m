addpath('/home/GUROBI.LIC/gurobi1300/linux64/matlab');
gurobi_setup;

addpath(genpath(pwd));

%% CyT data
load data_exp_R1.mat
load model_DFVB_R1.mat

modelo     = modelo_struct; 
theta_init = [0.3, 0.5, 100, 80, 0, 0.2, 0.7, 0.15, 150, 200, 1000, 50, 10000, 15, 0.5, 0.9, 0.8];

% ----------------------- User parameters ------------------------- %
Nboots = 1000;            % number of bootstrap replicates
includeFirst = true;      % force inclusion of the first time point in each replicate
rngBase = 1234;           % base seed for reproducibility

% --------------------------- Original bounds -------------------------- %
C_LB = [0.212	0.047	92.086	58.465	0.000	0.049];
C_UB = [0.419	0.052	121.601	92.555	0.000	0.358];

p_LB = [1e-5   1e-8    1e-2  1e-2  1e-3   1e-2   1e-2   1e-5 1e-5 1e-5 1e-5];
p_UB = [2.0      1     1e3   1e3   1e4    5e3    1e4    1e4  5    2    2];

% ------------------- Optimizer options (local) ------------------ %
opts = optimset('Display', 'off', ...                 
    'MaxIter', 150, ...           
    'TolFun', 1e-4); 
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
for b = 1:Nboots
    % Independent seed for each replicate
    s = RandStream('Threefry','Seed',seeds(b));
    RandStream.setGlobalStream(s);

    % --- Aligned resampling with replacement (same rows for all columns) ---
    [data_b, idx_b] = resampleAlignedWithReplacement(data_exp, includeFirst);

    % Objective function
    fun_b = @(x) LossFunction(x, modelo, data_b);

    % Optimization (sequential within each replicate)
    t0 = tic;
    [xopt_b, fval_b, exit_b, output] = fminsearchbnd(fun_b, theta_init, [C_LB, p_LB], [C_UB, p_UB], opts);
    Timing(b) = toc(t0);

    % Save results
    Xopt(b,:)   = xopt_b;
    Fval(b)     = fval_b;
    Exit(b)     = exit_b;
    Theta0(b,:) = theta_init;
end

%%
% ------------------------- Consolidation/saving ----------------------- %
varsParam = "p" + (1:npar);
Tresult = table((1:Nboots).', Exit, Fval, Timing, 'VariableNames', {'boot','exitflag','fval','time_s'});
Ttheta0 = array2table(Theta0, 'VariableNames', cellstr("theta0_"+varsParam));
Tparams = array2table(Xopt,  'VariableNames', cellstr(varsParam));
ResultadoBootstrap = [Tresult Ttheta0 Tparams];

save('Results_Bootstraping/Bootstraping_DFVB_CYT.mat', 'Xopt','Fval','Exit','Timing','Theta0','p_LB','p_UB');

fprintf("Saved CYT!\n")
