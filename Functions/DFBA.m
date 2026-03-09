function [dy, v, fva_min, fva_max] = DFBA(t,y,model,idx_met,idx_flux,params, FVA)
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
% persistent prevTime1
% persistent iteration_time1
% 
% if isempty(prevTime1)
%     prevTime1 = tic;
% end
% 
% if t == 0
%     iteration_time1 = [];
% end
% % Si no hay flag, estamos en una iteración
% iteration_time1 = [iteration_time1 toc(prevTime1)];  % Tiempo desde la última iteración
% 
% if t == 2
%     figure
%     plot(iteration_time1)
%     grid on
%     xlabel('Iteration')
%     ylabel('Time (s)')
%     title('Iteration Time dFBA')
%     % disp(sum(iteration_time1))
%     fprintf('Cantidad de iteraciones = %d \n',size(iteration_time1,2))
% end   
% 
% % Inicializar variable persistente
% prevTime1 = tic;  % Iniciar cronómetro

if nargin < 7
    FVA = false;
end

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
etoh_met= idx_met(7);

bio_flux  = idx_flux(1);
glc_flux  = idx_flux(2);
o2_flux   = idx_flux(3);
ac_flux   = idx_flux(4);
co2_flux  = idx_flux(5);
nh4_flux  = idx_flux(6);
etoh_flux = idx_flux(7);


%% States
X  = y(bio_met);   % Biomass
G  = y(glc_met);   % Substrate
O2 = y(o2_met);
Ac = y(ac_met);
EtOH = y(etoh_met);


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

if (EtOH<1e-2)
    model.lb(etoh_flux) = 0;
end

if (O2<1e-2)
    model.lb(o2_flux) = 0;
end

% Update lower bound of uptake substrate
model.lb(glc_flux)  = -rS;  


% Solution 

sol = optimizeCbModel(model,'max');

% FVA
if FVA
    if sol.stat == 1
        [fva_min, fva_max] = fluxVariability(model, 100, 'max', model.rxns(idx_flux),0);
    else
        fva_min = zeros(length(model.rxns(idx_flux)),1);
        fva_max = zeros(length(model.rxns(idx_flux)),1);
    end
else
    fva_min = zeros(length(model.rxns(idx_flux)),1);
    fva_max = zeros(length(model.rxns(idx_flux)),1);
end

if sol.stat==1
    mu    = sol.v(bio_flux);
    rS    = sol.v(glc_flux);
    rO2   = sol.v(o2_flux);
    rAc   = sol.v(ac_flux);
    rCO2  = sol.v(co2_flux);
    rNH4  = sol.v(nh4_flux);
    rEtOH = sol.v(etoh_flux);
    v     = sol.v';

    dy(bio_met) = mu*X;   %Biomass
    dy(glc_met) = rS*X;   %Glucose
    dy(o2_met)  = rO2*X;  %Oxygen
    dy(ac_met)  = rAc*X;  %Acetate
    dy(co2_met) = rCO2*X;
    dy(nh4_met) = rNH4*X;
    dy(etoh_met)= rEtOH*X;
    dy          = dy(:);

else
    dy = zeros(size(y,1),1);
    v  = zeros(size(model.lb))';
end
end