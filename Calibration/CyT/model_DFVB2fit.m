function [dxdt, da, v] = model_DFVB2fit(t, x, theta, T_data, model, index, PAN_struct, alpha_dim)

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

%% Zenteno Model

% Time varying temperature
T   = interp1(T_data.Time, T_data.Temperature, t, "nearest", "extrap")+273.15;     %Linear data interpolation

% State Variables
X    = x(1);           % Biomass  (g/L)
N    = x(2);           % Ammonia  (g/L)
G    = x(3);           % Glucose  (g/L)
F    = x(4);           % Fructose (g/L)
E    = x(5);           % Ethanol  (g/L)
PAN  = x(6);           % PAN      (g/L)

% Estimated Parameters
mu0    = theta(1);       % 1/h
Kn0    = theta(2);       % kg N/m3 (Cramer 2002)
Kg0    = theta(3);       % kg G/m3 (Cramer 2002)
Kf0    = theta(4);       % kg F/m3
Kie0   = theta(5);       % kg E/m3 (Boulton 1979)
Yxn    = theta(6);       % kg Bio/kg N (Zenteno 2010), biomass/nitrogen yield coeff
Yxg    = theta(7);       % kg bio/kg G (Zenteno 2010), biomass/glucose yield coeff
Yxf    = theta(8);       % kg bio/kg F (zenteno 2010), biomass/fructose yield coeff
Yef    = theta(9);       % kg E/ kg F (Zenteno 2010), ethanol/fructose yield coeff
betaG0 = theta(10);
betaF0 = theta(11);

%Fixed Parameters 
Kd0     = 0.00044;    % 1/h     (Zenteno 2010)
Kig0    = 42.5283;    % kg G/m3 (estimated salmon 2003, zenteno 2010)
Yeg     = 0.4517;     % kg E/ kg G (Zenteno 2010), ethanol/glucose yield coeff
Cde     = 0.0415;     % m3/kg E (Salmon 2003)
Etd     = 130000.0;   % kj/kmol
R       = 8.314;      % kJ/kmol/K (Boulton 1979)
Eac     = 59453.0;    % kJ/kmol (Boulton 1979)
Eafe    = 11000.0;    % kJ/kmol (Zenteno 2010, trial and error)
EaKn    = 46055.0;    % kJ/kmol (Boulton 1979)
EaKg    = 46055.0;    % kJ/kmol (Boulton 1979)
EaKf    = 46055.0;    % kJ/kmol (Boulton 1979)
EaKig   = 46055.0;    % kJ/kmol (Boulton 1979)
EaKie   = 46055.0;    % kJ/kmol (Boulton 1979)
Eam     = 37681.0;    % kJ/kmol (Boulton 1979)
m0      = 0.01;       % kgS/kg bio/h

%Molecular weights
MW_etoh     = 46.07; % g/mol
MW_nitrogen = 17.03; % g/mol

%% Constitutive Equations

% maximum specific growth rate, 1/h
mu_max    = mu0*exp(Eac*(T-300)/(300*R*T));

% maximum ethanol specific production rate from glucose, kg E/kg bio/h
betaG_max = betaG0*exp(Eafe*(T-296.15)/(296.15*R*T));

% maximum ethanol specific production rate from fructose, kg E/kg bio/h
betaF_max = betaF0*exp(Eafe*(T-296.15)/(296.15*R*T));

% Nitrogen compounds saturation constant, kg K/m3
Kn        = Kn0*exp(EaKn*(T-293.15)/(293.15*R*T));

% Glucose saturation constant, kg G/m3
Kg        = Kg0*exp(EaKg*(T-293.15)/(293.15*R*T));

% Fructose saturation constant, kg F/m3
Kf        = Kf0*exp(EaKf*(T-293.15)/(293.15*R*T));

% Glucose inhibition constant, kg G/m3
Kig       = Kig0*exp(EaKig*(T-293.15)/(293.15*R*T));

% Ethanol inhibition constant, kg E/m3
Kie       = Kie0*exp(EaKie*(T-293.15)/(293.15*R*T));

% maintenace coefficient, kg S/kg bio/h
m         = m0*exp(Eam*(T-293.3)/(293.3*R*T));

% ethanol production rate from glucose
beta_G    = betaG_max*(G/(G+Kg))*(Kie/(E+Kie));

% ethanol production rate from fructose
beta_F    = betaF_max*(F/(F+Kf))*(Kig/(G+Kig))*(Kie/(E+Kie));

% threshold temperature of thermal death
Td        = -0.0001*E^3+0.0049*E^2-0.1279*E+315.89;

% specific death rate
if T>=Td
    Kd    = Kd0*exp((Cde*E)+(Etd*(T-305.65))/(305.65*R*T));
else
    Kd    = 0;
end

%% Update Metabolic Model
index_mu   = index(1);
index_glc  = index(2);
index_fru  = index(3);
index_N    = index(4);
index_etoh = index(5);

% Nitrogen source
index_PAN = PAN_struct.index_DFVB;
MW_PAN    = PAN_struct.MW_PAN;

% Converting units
v_etoh_G = beta_G / MW_etoh * 1000;  % [mmol/gDW/h]
v_etoh_F = beta_F / MW_etoh * 1000;  % [mmol/gDW/h]
model.ub(index_etoh) = v_etoh_G + v_etoh_F;           % Ethanol production

% Limiting substrate consumption if depleted
if G <= 0
    model.lb(index_glc) = 0;
end

if F <= 0
    model.lb(index_fru) = 0;
end

% Distributing the consumption rate in NH4 (N) and PAN
vN_ggDW = (mu_max / Yxn) * ((N+PAN) / ((N+PAN) + Kn));   % [g/gDW/h]

if N <= 0
    model.lb(index_N) = 0;
else
    model.lb(index_N) = - vN_ggDW * (N/(N + PAN)) * 1000/ MW_nitrogen; 
end

if PAN <= 0
    model.rhs(end) = 0;
else 
    model.rhs(end) =  vN_ggDW * (PAN/(N + PAN)) * 1000 / mean(MW_PAN);
end

%% Solve DFVB
sol = gurobi(model, gurobiparams);
sol_factible = sol;

% Extract mu
if strcmp(sol.status, 'OPTIMAL')
    mu           = sol.x(index_mu);
    da           = sol.x(end-alpha_dim+1:end)';
    v_pan        = sum(sol.x(index_PAN) .* MW_PAN(1)' ./ 1000);
    
    % Convertir tasas de mmol/gDW/h → g/L/h
    v_N_gh    = sol.x(index_N) * X * MW_nitrogen / 1000;
    v_etoh_gh = sol.x(index_etoh) * X * MW_etoh     / 1000;

    % Differential equations
    dXdt   =  mu*X - Kd*X;                                  % Biomass 
    dNdt   =  v_N_gh;                                       % Nitrogen
    dGdt   = -( (mu/Yxg) + (beta_G/Yeg) + m*(G/(G+F)) )*X;  % Glucose
    dFdt   = -((mu/Yxf) + (beta_F/Yef) + m*(F/(G+F)))*X;    % Fructose
    dEdt   = v_etoh_gh;                                     % Ethanol
    dPANdt = -v_pan * X;
    
    dxdt   = [dXdt;dNdt;dGdt;dFdt;dEdt; dPANdt];
    v      = sol.x(1:end-alpha_dim);
else
    dxdt   = zeros(6,1);
    da     = zeros(alpha_dim,1);
    v      = zeros(size(dxdt,1),1);
end

end



