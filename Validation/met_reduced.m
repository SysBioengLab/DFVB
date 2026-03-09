function [mets] = met_reduced(model, idxs)

mets = find(all(model.S(:, idxs ) == 0, 2));

end