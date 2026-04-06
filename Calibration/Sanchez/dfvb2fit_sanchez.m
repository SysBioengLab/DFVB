function [dxdt, da, v] = dfvb2fit_sanchez(t, x, theta, model, index, alpha_dim)

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

X = x(1); %Biomass
G = x(2); %Glucose
E = x(3); %Ethanol

%% Update Metabolic Model
index_mu   = index(1);
index_glc  = index(2);
index_o2   = index(3);
index_etoh = index(4);

MW_glc  = 180.156/1e3;
MW_etoh = 46.07/1e3; 
MW_lac  = 90.08/1e3;
MW_gly  = 92.09382/1e3;
MW_cit  = 192.124/1e3;


% Estimated parameters
v_glc   = theta(1);
K_G     = theta(2);
K_iE    = theta(3);
lb_o2   = theta(4);
f_gly   = theta(5);
f_cit   = theta(6);
f_lac   = theta(7);


% Limiting substrate consumption if depleted
if G > 1e-3
    rS = v_glc*G/((K_G+G)*(1+E/K_iE));
else
    rS = 0;
end

% if E>=0
%     rE = v_et*E/(K_E + E);
% else
%     rE = 0;
% end

model.lb(index_glc)  = -rS;
model.lb(index_o2)   = -lb_o2;


%% Solve DFNA

sol = gurobi(model, gurobiparams);
sol_factible = sol;

% Extract mu
if strcmp(sol.status, 'OPTIMAL')
    mu   = sol.x(index_mu);
    rS   = sol.x(index_glc);
    rE   = sol.x(index_etoh);


    da = sol.x(end-alpha_dim+1:end)';
    v  = sol.x(1:end-alpha_dim);

    % Differential equations
    dXdt   =  mu * X;              % Biomass 
    dGdt   =  rS * X * MW_glc ;      % Glucose
    dEdt   =  rE * X * MW_etoh ; % Ethanol
    dGly   =  f_gly  * mu * X * MW_gly  ; % Glycerol
    dCit   =  f_cit  * mu * X * MW_cit  ; % Citrate
    dLac   =  f_lac  * mu * X * MW_lac ; % Lactate
    
    dxdt   = [dXdt;dGdt;dEdt;dGly;dCit;dLac];


else
    dxdt = zeros(6,1);
    da  = zeros(alpha_dim,1);
    v   = zeros(262,1);
end




end




