function [dy, v] = DFBA_gurobi(t,y,model,idx_met,idx_flux,params)
%--------------------------------------------------------------------------
% Dynamic flux balance analysis 
% Resolve an optimization problem for pseudoestationary state S*v, with 
% ODEs for biomass, substrate (glucose) and product (protein). Based on
% Zhao et al 2017, Example 1: Spirallus
%
% Inputs:
%        t      : vector of time
%        y      : vector with actual state of metabolites
%        model  : standard form of metabolic network
%        ix_glc : index of glucose flux in matrix S
%        ix_bio : index of biomass flux in matrix S
%        ix_prot: index of protein flux in matrix S
%        params : (qs_max, k_S) for Monod
%        FVA    : true for run flux variability analysis at each iteration
%--------------------------------------------------------------------------

% %% Rutine for extract iteration time
% % Variable persistente para almacenar el tiempo anterior
% persistent prevTime
% persistent iteration_time
% 
% if isempty(prevTime)
%     prevTime = tic;
% end
% 
% if t == 0
%     iteration_time = [];
% end
% % Si no hay flag, estamos en una iteración
% iteration_time = [iteration_time toc(prevTime)];  % Tiempo desde la última iteración
% 
% if t == 2
%     figure
%     plot(iteration_time)
%     grid on
%     xlabel('Iteration')
%     ylabel('Time (s)')
%     title('Iteration Time dCBA')
%     fprintf('Cantidad de iteraciones = %d \n',size(iteration_time,2))
%     % disp(sum(iteration_time))
%     % disp(size(iteration_time,2))
% end   
% 
% % Inicializar variable persistente
% prevTime = tic;  % Iniciar cronómetro
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
q_max  = params(1);
K_S    = params(2);

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

if (G>1e-2)
    rS   = q_max*G/(K_S+G); % mmol/gdcw/h
else
    rS   = 0;
end

if (Ac<1e-2)
    model.lb(ac_flux) = 0;
end


% Update lower bound of uptake substrate
model.lb(glc_flux)  = -rS;  


% Solution 
sol = gurobi(model, gurobiparams);
sol_factible = sol;

% persistent alphas
% if t == 0
%     alphas = [];
% end

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
    dy(o2_met)  = rO2*X;  %Oxygen
    dy(ac_met)  = rAc*X;  %Acetate
    dy(co2_met) = rCO2*X;
    dy(nh4_met) = rNH4*X;
    dy(etoh_met)= rEtoh*X;
    dy          = dy(:);

    
else
    dy = zeros(size(y,1),1);
    v  = zeros(size(model.lb))';
end


end