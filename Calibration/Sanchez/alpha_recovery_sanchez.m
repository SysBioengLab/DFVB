function [alphas, fluxes] = alpha_recovery_sanchez(T, Y, theta, model, index, alpha_dim)

% Alphas recovery (conversions)

alphas   = zeros(length(T),alpha_dim);
fluxes   = zeros(length(T),size(model.lb,1)-alpha_dim);

for i=1:length(T)
    [~, a, v]    = dfna2fit_sanchez(T(i), Y(i,:), theta, model, index, alpha_dim);
    alphas(i,:)  = a;
    fluxes(i,:)  = v;

end