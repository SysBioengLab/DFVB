function [fluxes, fva_min, fva_max] = Variable_recovery_FBA(T,Y,model,idx_met,idx_flux,params, fva)

% Alphas recovery (conversions)

fluxes   = [];
fva_min  = [];
fva_max  = [];

for i=1:length(T)
    [~, v, min_range, max_range] = DFBA(T(i), Y(i,:), model, idx_met, idx_flux, params, fva);
    fva_min = [fva_min min_range];
    fva_max = [fva_max max_range];
    try
        fluxes  = [fluxes  v'];
    
    catch ME
        break
    end
end