function [dy,da, v] = DFNA(t,y,model,index,params,alpha_dim)
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
% persistent alphas
% persistent flag
% flag = true;

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
G    = y(2) ;   % Substrate
O2   = y(3) ;   % Oxigeno
Ac   = y(4) ;   % Acetato
Nh4  = y(5) ;   % Amonio
EtoH = y(6) ;   % Etanol
CO2  = y(7) ;   % CO2


%% MODEL

% Calculate specific rates based on the substrates present
if (G>1e-2)
    rS   = v_max*G/(K_S+G); % mmol/gdcw/h
else
    rS   = 0;
end

if (Ac<1e-2)
    model.lb(index(4)) = 0;
end

% Update lower bound of uptake substrate
model.lb(index(2))   = -rS;  


% Solution 
sol = gurobi(model, gurobiparams);
sol_factible = sol;


% if t == 0
%     alphas = [];
% end

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
    dy(3) = rO2*X;
    dy(4) = rAc*X;
    dy(5) = rN*X;
    dy(6) = rEt*X;
    dy(7) = rCO2*X;
    dy    = dy(:);
    
else
    dy = zeros(size(y,1),1);
    da = zeros(alpha_dim,1);
    v =  zeros(size(dy,1),1);

end

% if size(da,2) == 1
%     alphas = [alphas da];
% else
%     alphas = [alphas da'];
% end

% if t == 2 & flag
%     assignin('base', 'alpha_final', alphas);
%     flag = false;
% end

end