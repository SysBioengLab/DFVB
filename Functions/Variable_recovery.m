function [alphas, fluxes] = Variable_recovery(T,Y,model,index,params,alpha_dim)

% Alphas recovery (conversions)

alphas   = zeros(length(T),alpha_dim);
fluxes   = zeros(length(T),size(model.lb,1)-alpha_dim);

for i=1:length(T)
    [~, a, v]    = DFNA(T(i),Y(i,:),model,index,params,alpha_dim);
    alphas(i,:)  = a;
    fluxes(i,:)  = v;
end