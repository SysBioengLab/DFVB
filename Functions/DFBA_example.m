function [dy, v] = DFBA_example(t,y,model,index,params)
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
% Outputs:
%        dy     : vector with derivatives of metabolites
%        v      : vector with fluxes of the system
%
%--------------------------------------------------------------------------

%% Parameters
v_max  = params(1);
K_S    = params(2);

%% States
X    = y(1) ;   % Biomass
S    = y(2) ;   % Substrate
P    = y(2) ;   % Product


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
sol = optimizeCbModel(model);



if strcmp(sol.origStat, 'OPTIMAL')
    mu    = sol.x(index(1));
    rS    = sol.x(index(2));
    rP    = sol.x(index(3));

    
    % Define ODE system
    dy(1) = mu*X; 
    dy(2) = rS*X;
    dy(3) = rP*X;
    dy    = dy(:);
    
else
    dy = zeros(size(y,1),1);
    v =  zeros(size(dy,1),1);

end

end