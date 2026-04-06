function [dy, v] = DFBA_gurobi(t,y,model,idx_met,idx_flux,params)
%--------------------------------------------------------------------------
% Dynamic flux balance analysis
% Solves a pseudo-steady-state optimization problem S*v = 0 to obtain the
% intracellular flux distribution and couple it to extracellular dynamics
% for biomass, substrate (glucose), oxygen, acetate, CO2, NH4 and ethanol.
% The model is solved using Gurobi (doens't require COBRA Toolbox).
%
% Inputs:
%   t        : current time
%   y        : vector of extracellular state concentrations
%   model    : Gurobi model describing the metabolic network
%   idx_met  : indices of metabolites in y
%   idx_flux : indices of fluxes in the model solution vector
%   params   : model parameters [q_max, K_S, vO2_max, K_O2, kla]
%
% Outputs:
%   dy       : time derivative of the state vector
%   v        : optimal flux vector returned by the LP solver 
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
q_max   = params(1);
K_S     = params(2);
vO2_max = params(3);  
K_O2    = params(4);  
kla     = params(5);
%% Index
bio_met = idx_met(1);
glc_met = idx_met(2);
o2_met  = idx_met(3);
ac_met  = idx_met(4);
co2_met = idx_met(5);
nh4_met = idx_met(6);
etoh_met = idx_met(7);

bio_flux = idx_flux(1);
glc_flux = idx_flux(2);
o2_flux  = idx_flux(3);
ac_flux  = idx_flux(4);
co2_flux = idx_flux(5);
nh4_flux = idx_flux(6);
etoh_flux = idx_flux(7);

%% States
X  = y(bio_met);   % Biomass
G  = y(glc_met);   % Substrate
O2 = y(o2_met);
Ac = y(ac_met);

%% MODEL

% Calculate specific rates based on the substrates present

if (Ac<1e-4)
    model.lb(ac_flux) = 0;
end

% Update lower bound of uptake substrate
model.lb(glc_flux)  = -q_max*G/(K_S+G);


% Solution 
sol = gurobi(model, gurobiparams);
sol_factible = sol;


if strcmp(sol.status, 'OPTIMAL')
    mu    = sol.x(bio_flux);
    rS    = sol.x(glc_flux);
    rO2   = sol.x(o2_flux);
    rAc   = sol.x(ac_flux);
    rCO2  = sol.x(co2_flux);
    rNH4  = sol.x(nh4_flux);
    rEtoh = sol.x(etoh_flux);
    v    = sol.x';

    dy(bio_met) = mu*X;   %Biomass
    dy(glc_met) = rS*X;   %Glucose
    dy(o2_met)  = rO2*X + kla*(0.3125 - O2);  %Oxygen
    dy(ac_met)  = rAc*X;  %Acetate
    dy(co2_met) = rCO2*X;
    dy(nh4_met) = rNH4*X;
    dy(etoh_met)= rEtoh*X;
    dy          = dy(:);

    
else
    dy = zeros(size(y,1),1);
    v  = zeros(size(model.lb))';
    if (0.3125 - O2) > 1e-4
        dy(o2_met)  = kla*(0.3125 - O2);
    end
end


end