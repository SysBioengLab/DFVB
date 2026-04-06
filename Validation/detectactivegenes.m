function geneTable = detectactivegenes(model, act_rxns, namefile)  
% --------------------------------------------------------------------
% Detect active genes based on active reactions in a metabolic model
% Inputs:
%   model - COBRA model structure with fields 'genes' and 'grRules'
%   act_rxns - vector of indices of active reactions
%   namefile (optional) - base name for output CSV file (without extension)
% Output:
%   geneTable - table with columns 'GeneID' and 'Status' indicating active,
%               ambiguous, or inactive genes
% --------------------------------------------------------------------
    if nargin < 3
        namefile = '';
    end

    % Initialize map to classify genes
    geneStatusMap = containers.Map(model.genes, repmat({'inactive'}, size(model.genes)));
    
    % Iterate only the active reactions
    for i = act_rxns(:)'  
        grRule = strtrim(model.grRules{i});
    
        if ~isempty(grRule)
            % Detect if contains ' and ' (deterministic case)
            if contains(grRule, ' and ')
                % Extract genes from the rule
                genesInRule = regexp(grRule, '[a-zA-Z0-9_]+', 'match');
                genesInRule = intersect(genesInRule, model.genes);
                for g = 1:length(genesInRule)
                    geneID = genesInRule{g};
                    geneStatusMap(geneID) = 'active'; 
                end
            else
                % 'or' case or single gene mark as ambiguous if not already active
                genesInRule = regexp(grRule, '[a-zA-Z0-9_]+', 'match');
                genesInRule = intersect(genesInRule, model.genes);
                for g = 1:length(genesInRule)
                    geneID = genesInRule{g};
                    if ~strcmp(geneStatusMap(geneID), 'active')
                        geneStatusMap(geneID) = 'ambiguous';
                    end
                end
            end
        end
    end
    
    % Create table with results
    geneList = model.genes(:);
    statusList = values(geneStatusMap, geneList);
    statusList = statusList(:);  
    
    geneTable = table(geneList, statusList, 'VariableNames', {'GeneID', 'Status'});
    
    % Show summary
    disp('Summary of genes by status:');
    disp(groupcounts(geneTable.Status));
    
    % save as CSV
    if ~isempty(namefile)
        writetable(geneTable, strcat(namefile,'.csv'));
    end
end