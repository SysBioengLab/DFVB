function [dy,da, v] = DFVB_o2_limited(t,y,model,index,params,alpha_dim)
%--------------------------------------------------------------------------
% Dynamic Flux Vector Balancing (DFVB) - Oxygen-Limited Condition
%
% Computes the time derivatives of the macroscopic states by solving a
% pseudo-steady-state metabolic optimization problem at each time point.
% This routine is specifically formulated for oxygen-limited conditions,
% where oxygen uptake is constrained by the extracellular O2 availability
% and mass transfer through kLa.
%
% Inputs:
%   t         - Current simulation time
%   y         - Current state vector of extracellular metabolites
%   model     - Gurobi model describing the metabolic network
%   index     - Indices of the flux variables in the solution vector
%   params    - Model parameters [v_max, K_S, vO2_max, K_O2, kla]
%   alpha_dim - Number of alpha variables in the optimization solution
%
% Outputs:
%   dy        - Time derivative of the state vector
%   a         - Alpha variables extracted from the optimization solution
%   v         - Flux vector extracted from the optimization solution
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
v_max  = params(1);   % vmax glucosa
K_S    = params(2);   % Ks glucosa
vO2_max = params(3);  % vmax O2 (mmol/gDCW/h) 
K_O2    = params(4);  % Ks O2 (mM)
kla     = params(5);

%% States
X    = y(1) ;   % Biomass
G    = y(2) ;   % Substrate
O2   = y(3) ;   % Oxigeno
Ac   = y(4) ;   % Acetato
Nh4  = y(5) ;   % Amonio
EtoH = y(6) ;   % Etanol
CO2  = y(7) ;   % CO2


%% MODEL

if (Ac<1e-4)
    model.lb(index(4)) = 0;
end 

% if (EtoH<1e-4)
%     model.lb(index(6)) = 0;
% end 


% Update lower bound of uptake substrate
model.lb(index(2))   = -v_max*G/(K_S+G);
model.lb(index(3))   = - vO2_max * O2 / (K_O2 + O2);


% Solution 
sol = gurobi(model, gurobiparams);
sol_factible = sol;


if strcmp(sol.status, 'OPTIMAL')
    mu    = sol.x(index(1));
    rS    = sol.x(index(2));
    rO2   = sol.x(index(3));
    rAc   = sol.x(index(4));
    rN    = sol.x(index(5));
    rEt   = sol.x(index(6));
    rCO2  = sol.x(index(7));
    da    = sol.x(end-alpha_dim+1:end)';
    v     = sol.x(1:end-alpha_dim);
    
    % Define ODE system
    dy(1) = mu*X; 
    dy(2) = rS*X;
    dy(3) = rO2*X + kla * (0.3125 - O2);
    dy(4) = rAc*X;
    dy(5) = rN*X;
    dy(6) = rEt*X;
    dy(7) = rCO2*X;
    dy    = dy(:);
    
else
    dy    = zeros(size(y,1),1);
    dy(3) = kla * (0.3125 - O2);
    da    = zeros(alpha_dim,1);
    v     =  zeros(size(dy,1),1);

end


end