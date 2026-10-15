function [mets] = met_reduced(model, idxs)
% -------------------------------------------------------------------------------------------
% This function identifies the metabolites that are not involved in any of the reactions specified by idxs.
% Inputs:
%   model: A COBRA model structure containing the stoichiometric matrix S.
%   idxs: A vector of reaction indices for which we want to find the associated metabolites.
%
% Output:
%   mets: A vector of metabolite indices that are not involved in any of the reactions
% specified by idxs.
% -------------------------------------------------------------------------------------------

mets = find(all(model.S(:, idxs ) == 0, 2));

end