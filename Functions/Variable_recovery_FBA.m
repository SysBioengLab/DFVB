function [fluxes, fva_min, fva_max] = Variable_recovery_FBA(T,Y,model,idx_met,idx_flux,params, fva)
% --------------------------------------------------------------------------
% Performs variable recovery for flux balance analysis (FBA) using the DFBA framework.
% Inputs:
%   T         - Time points for the simulation.
%   Y         - State variables corresponding to the time points.
%   model     - COBRA metabolic model structure.
%   idx_met   - Indices of the metabolites in the model.
%   idx_flux  - Indices of the fluxes in the model.
%   params    - Additional parameters for the DFBA simulation.
%   fva       - Flag indicating whether to perform flux variability analysis (FVA).
% Outputs:
%   fluxes    - Recovered fluxes for each time point.
%   fva_min   - Minimum flux values from FVA for each time point (if FVA is performed).
%   fva_max   - Maximum flux values from FVA for each time point (if FVA is performed).
% --------------------------------------------------------------------------

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