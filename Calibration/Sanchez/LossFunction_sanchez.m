function loss = LossFunction_sanchez(parameters, modelo, data_exp)

% DFNA execution

%ODE configuration
options   = odeset('AbsTol',1e-6, 'NonNegative',[1, 2, 3, 4, 5, 6]);  % Opciones del ODE

y0     = parameters(1:6);
theta  = parameters(7:end);
t_span = [0, 12];


% Simulate first period
[T_sol, Y_sol] = ode113(@dfvb2fit_sanchez, t_span,y0,options,theta, modelo.model_DFNA, modelo.index, modelo.alpha_dim);

% Calculate error
sol_interp  = interp1(T_sol, Y_sol, data_exp.t_exp, "spline", "extrap");

denom = data_exp.conc;
denom(denom == 0) = 1;  % evitar división por cero
denom(denom < 0)  = 1;

loss_vector =  [1 1 1 1 1 1] .* mean(((sol_interp - data_exp.conc) ./ max(denom)).^2, 1, 'omitnan');
loss        = sum(loss_vector);

fprintf('Loss value = %.4f\n', loss)


%%


