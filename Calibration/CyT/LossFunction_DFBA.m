function loss = LossFunction_DFBA(parameters, modelo, data_exp)

% DFNA execution

%ODE configuration
options   = odeset('AbsTol',1e-6, 'NonNegative',[1, 2, 3, 4, 5, 6]);  % Opciones del ODE

y0    = parameters(1:6);
theta = parameters(7:end);

% Simulate first period
[T1, Y1] = ode113(@model_dFBA2fit,modelo.t_rango1,y0,options,theta, modelo.Temperature, modelo.model_DFBA, modelo.index, modelo.PAN_struct);

% Simulate second period
y0_2     = Y1(end,:);
y0_2(2)  = modelo.y0_2(2);
[T2, Y2] = ode113(@model_dFBA2fit,modelo.t_rango2,y0_2,options,theta,modelo.Temperature, modelo.model_DFBA, modelo.index, modelo.PAN_struct);

% Join the solution
T_sol = [T1(1:end-1); T2];
Y_sol = [Y1(1:end-1, :); Y2];

% Calculate error

sol_interp  = interp1(T_sol, Y_sol, data_exp.t_exp, "spline", "extrap");

denom = data_exp.conc;
denom(denom == 0) = 1;  % evitar división por cero
denom(denom < 0)  = 1;

loss_vector =  [1 1 1 1 1 1] .* mean(((sol_interp - data_exp.conc) ./ max(denom)).^2, 1, 'omitnan');
% loss_vector = [10 10 0.1 0.1 0.1 0.1] .* mean((sol_interp - data_exp.conc).^2./max(data_exp.conc, [], 1, 'omitnan'), 1, 'omitnan') ; 
loss        = sum(loss_vector);

fprintf('Loss value = %.4f\n', loss)


%%


