function [dy,da, v] = DFNA_example(t,y,model,index,params,alpha_dim)
%--------------------------------------------------------------------------
% Dynamic conversion balance analysis 
% Resolve an optimization problem for pseudoestationary state S*v, with 
% ODEs for biomass, substrate (glucose) and product (protein). Based on
% Zhao et al 2017, Example 1: Spirallus
%
% Inputs:
%        t      : vector of time
%        y      : vector with actual state of metabolites
%        model  : standard form of metabolic network
%        index  : vector with metabolic index in solution
%        params : (qs_max, k_S) 
%
%--------------------------------------------------------------------------

%% Rutine for WarmUp

persistent sol_factible
if t == 0
    sol_factible = struct();
elseif strcmp(sol_factible.status,'OPTIMAL')
    model.vbasis = sol_factible.vbasis;
    model.cbasis = sol_factible.cbasis;
end

gurobiparams.OutputFlag  = 0;
gurobiparams.LPWarmStart = 1;

%% Parameters
v_max  = params(1);
K_S    = params(2);

%% States
X    = y(1) ;   % Biomass
S    = y(2) ;   % Substrate
P    = y(3) ;   % Product


%% MODEL

% Calculate specific rates based on the substrates present
if (S>1e-2)
    rS   = v_max*S/(K_S+S); % mmol/gdcw/h
else
    rS   = 0;
end

% Update lower bound of uptake substrate
model.lb(index(2))   = -rS;  


% Solution 
sol = gurobi(model, gurobiparams);
sol_factible = sol;

% persistent alphas
% if t == 0
%     alphas = [];
% end

if strcmp(sol.status, 'OPTIMAL')
    mu    = sol.x(index(1));
    rS    = sol.x(index(2));
    rP    = sol.x(index(3));

    
    % Define ODE system
    dy(1) = mu*X; 
    dy(2) = rS*X;
    dy(3) = rP*X;
    dy    = dy(:);
    da    = sol.x(end-alpha_dim+1:end)';
    v     = sol.x(1:end-alpha_dim);  
else
    dy = zeros(size(y,1),1);
    da = zeros(alpha_dim,1);
    v =  zeros(size(dy,1),1);

end

end