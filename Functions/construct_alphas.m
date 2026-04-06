function problem = construct_alphas(model,N)
%--------------------------------------------------------------------------
% Dynamic conversion balance analysis model constructor
% It constructs the optimization model using sparse null basis
%
% Inputs:
%        model  : metabolic model in standard form
%        N      : null basis of the stoichiometric matrix
%--------------------------------------------------------------------------
%Extract data from model

lb = model.lb;  % Lower boundaries      
ub = model.ub;  % Upper boundaries
c  = model.c ;  % Objective function

%Build optimization problem and extract parameters

if isempty(N)
    N = sparseNull(model.S, 1e-8);
end

exchanges  = find(findExcRxns(model, true));
nullS      = N(exchanges,:);
intS       = N(setdiff(1:length(ub), exchanges),:);

% Matrix of optimization problem
Iex     = eye(size(nullS,1),size(nullS,1));
Zero_ex = zeros(size(intS,1), size(Iex,2));

A  = sparse([   -Iex, nullS; ...    % Nint * alpha - v_ex = 0
              Zero_ex, intS; ...    % Nint * alpha <= v_max
              Zero_ex, intS]);      % Nint * alpha >= v_min

% Boundaries
ub_ex = ub(exchanges);
lb_ex = lb(exchanges);

int_ub = ub(setdiff(1:length(ub), exchanges));
int_lb = lb(setdiff(1:length(lb), exchanges));

alfa_ub      = ones(size(nullS,2),1)* 1e6;
alfa_lb      = ones(size(nullS,2),1)*-1e6;

ub2       = [ ub_ex; alfa_ub];
lb2       = [ lb_ex; alfa_lb];

% Vector b (right size of A*x = b)

b = [zeros(size(ub_ex)); int_ub; int_lb];
csense = [repmat('=', size(ub_ex,1), 1); repmat('<', size(int_ub,1), 1); repmat('>', size(int_lb,1), 1)];

%Objective function
c2  = [ c(exchanges) ; zeros(size(nullS,2),1)];

%Optimization model construction
problem.A          = A;
problem.S          = A;
problem.ub         = ub2;
problem.lb         = lb2;
problem.obj        = c2;
problem.rhs        = b;
problem.modelsense = 'Max'; 
problem.metNames   = model.rxns(exchanges); 
problem.sense      = csense;
problem.nullS      = nullS;
problem.alphalb    = alfa_lb;
problem.alphaub    = alfa_ub;
problem.rxns       = model.rxns(exchanges); 
problem.intrxns    = model.rxns(setdiff(1:length(ub), exchanges));
problem.intN       = intS;

end

