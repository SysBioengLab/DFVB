function [lb_alpha, ub_alpha] = alpha_range(model)

n        = length(model.alphalb);
lb_alpha = zeros(n,1);
ub_alpha = zeros(n,1);
params.OutputFlag  = 0;
params.LPWarmStart = 1;

f_obj          = zeros(size(model.obj));
f_obj(end-n+1) = 1;
model.obj      = f_obj;
presol         = gurobi(model, params);
model.vbasis       = presol.vbasis;
model.cbasis       = presol.cbasis;

for i=1:n

    f_obj    = zeros(size(model.obj));

    % Upper bound
    f_obj(end-n+i)  = 1;
    model.obj       = f_obj;
    sol             = gurobi(model, params);
    ub_alpha(i)     = sol.x(end-n+i);

    % Lower bound
    f_obj(end-n+i)  = -1;
    model.obj       = f_obj;
    sol             = gurobi(model, params);
    lb_alpha(i)     = sol.x(end-n+i);

end

end