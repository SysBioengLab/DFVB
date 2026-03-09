function [alphas, fluxes] = alpha_recovery(T, Y, theta, T_data, model, index, PAN_struct, alpha_dim)

% Alphas recovery (conversions)

alphas   = zeros(length(T),alpha_dim);
fluxes   = zeros(length(T),size(model.lb,1)-alpha_dim);

for i=1:length(T)
    [~, a, v]    = model_dCBA2fit(T(i), Y(i,:), theta, T_data, model, index, PAN_struct, alpha_dim);
    alphas(i,:)  = a;
    fluxes(i,:)  = v;

end