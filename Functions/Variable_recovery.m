function [alphas, fluxes] = Variable_recovery(T,Y,model,index,params,alpha_dim)
%--------------------------------------------------------------------------
% Performs variable recovery for the DFVB framework, extracting both the fluxes and the alpha variables
% Inputs:
%   T         - Time points for the simulation.
%   Y         - State variables corresponding to the time points.
%   model     - COBRA metabolic model structure.
%   index     - Indices of the variables in the model.
%   params    - Additional parameters for the DFVB simulation.
%   alpha_dim - Dimension of the alpha variables.
%
% Outputs:
%   alphas    - Recovered alpha variables for each time point.
%   fluxes    - Recovered fluxes for each time point.
% --------------------------------------------------------------------------

% Alphas recovery (conversions)

alphas   = zeros(length(T),alpha_dim);
fluxes   = zeros(length(T),size(model.lb,1)-alpha_dim);

for i=1:length(T)
    [~, a, v]    = DFVB(T(i),Y(i,:),model,index,params,alpha_dim);
    alphas(i,:)  = a;
    fluxes(i,:)  = v;
end