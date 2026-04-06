function [dy,a, v] = DFVB_example(t,y,model,index,params,alpha_dim)
%--------------------------------------------------------------------------
% Dynamic Flux Vector Balancing (DFVB)
% Solves a pseudo-steady-state flux balance problem coupled to ODEs for
% biomass, substrate, and product dynamics.
% The formulation is adapted from Zhao et al. (2017), Example 1: Spirallus.
%
% Inputs:
%   t      : Current time
%   y      : Current state vector [biomass; substrate; product]
%   model  : Metabolic network in standard form
%   index  : Indices of the relevant fluxes in the solution vector
%   params : Model parameters [qS_max, K_S]
%   alpha_dim : Number of alpha variables in the optimization solution
%
% Outputs:
%   dy     : Time derivatives of the state variables
%   a      : State of alpha variables
%   v      : Flux vector of the system
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

if strcmp(sol.status, 'OPTIMAL')
    mu    = sol.x(index(1));
    rS    = sol.x(index(2));
    rP    = sol.x(index(3));

    
    % Define ODE system
    dy(1) = mu*X; 
    dy(2) = rS*X;
    dy(3) = rP*X;
    dy    = dy(:);
    a     = sol.x(end-alpha_dim+1:end)';
    v     = sol.x(1:end-alpha_dim);  
else
    dy = zeros(size(y,1),1);
    a  = zeros(alpha_dim,1);
    v  =  zeros(size(dy,1),1);

end

end