function problem = build_DFVB_example(model, N)
%--------------------------------------------------------------------------
% Builds the optimization problem for the DFVB framework based on the input metabolic model
% and the null basis matrix computed with sparseNull function.
% This function is specifically for the toy model of Zhao et al. (2017), Example 1: Spirallus.
%
% Inputs:
%   model - COBRA metabolic model structure.
%   N     - Null basis matrix computed from the stoichiometric matrix of the model.
%
% Outputs:
%   problem - Structure containing the optimization problem for DFVB, including matrices, 
%             bounds, objective function, and other relevant information.
%--------------------------------------------------------------------------
alpha_model = construct_alphas_example(model, N);

% Extract data from model

lb = model.lb;  % Lower boundaries      
ub = model.ub;  % Upper boundaries
c  = model.c ;  % Objective function

exchanges  = [1, 6, 7];
% Extract relevant rows from intN
[intN, lb_int, ub_int, nullS, lb_ex, ub_ex, c_correct, inac_int, inac_ex]  = forced_idx(alpha_model);

% Matrix of optimization problem
Iex     = eye(size(nullS,1) , size(nullS,1));
Zero_ex = zeros(size(intN,1), size(Iex,2)  );

A  = sparse([   -Iex, nullS; ...    % Nint * alpha - v_ex = 0
              Zero_ex, intN; ...    % Nint * alpha <= v_max
              Zero_ex, intN]);      % Nint * alpha >= v_min

% Alpha boundaries
% tic;
% [alpha_lb, alpha_ub] = alpha_range(alpha_model);
% disp(toc)

% ub2       = [ ub_ex; alpha_ub];
% lb2       = [ lb_ex; alpha_lb];
% 
ub2       = [ ub_ex; ones(size(alpha_model.alphalb))*  1e3];
lb2       = [ lb_ex; ones(size(alpha_model.alphalb))* -1e3];

% Vector b (right size of A*x = b)

b      = [zeros(size(ub_ex)); ub_int; lb_int];

% csense = [repmat('=', size(ub_ex,1), 1)];
csense = [repmat('=', size(ub_ex,1), 1); repmat('<', size(ub_int,1), 1); repmat('>', size(lb_int,1), 1)];


%Objective function
c2  = [ c_correct ; 0 * ones(size(alpha_model.alphalb))];

%Optimization model construction
problem.A               = A;
problem.S               = A;
problem.ub              = ub2;
problem.lb              = lb2;
problem.obj             = c2;
problem.rhs             = b;
problem.modelsense      = 'Max'; 
problem.sense           = csense;
problem.nullS           = nullS;
problem.intN            = intN;
problem.rxns            = alpha_model.rxns(setdiff(1:size(alpha_model.rxns,1), inac_ex));
problem.intrxns         = alpha_model.intrxns(setdiff(1:size(alpha_model.intrxns,1), inac_int));
problem.inac_ex         = alpha_model.rxns(inac_ex);
problem.inac_int        = alpha_model.intrxns(inac_int);
problem.rxns_sin_proc   = alpha_model.rxns;
end