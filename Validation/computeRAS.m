function ras = computeRAS(model, expr)
    % Step 1: count how many reactions each gene participates in
    numGenes = length(model.genes);
    geneUseCount = zeros(numGenes, 1);
    for i = 1:numel(model.grRules)
        rule = model.grRules{i};
        if isempty(rule), continue, end
        genesInRule = regexp(rule, '[^\s()]+', 'match');
        for j = 1:length(genesInRule)
            g = genesInRule{j};
            idx = find(strcmp(model.genes, g));
            if ~isempty(idx)
                geneUseCount(idx) = geneUseCount(idx) + 1;
            end
        end
    end
    % Step 2: normalize gene expression by gene usage
    exprNorm = expr ./ max(geneUseCount, 1);
    geneExprMap = containers.Map(model.genes, exprNorm);

    % Step 3: compute RAS per reaction by evaluating logical rules
    ras = nan(numel(model.grRules), 1);
    for i = 1:numel(model.grRules)
        rule = model.grRules{i};
        if isempty(rule), continue, end
        try
            exprTree = parseRule(rule, geneExprMap);
            ras(i) = eval(exprTree);
        catch
            ras(i) = NaN;
        end
    end
end

function exprStr = parseRule(rule, geneExprMap)
    % Convert logical operators to symbols and substitute gene values
    rule = lower(rule);
    rule = strrep(rule, 'and', '&');
    rule = strrep(rule, 'or', '|');
    genes = unique(regexp(rule, '[\w]+', 'match'));
    for i = 1:length(genes)
        g = genes{i};
        if isKey(geneExprMap, g)
            val = geneExprMap(g);
        else
            val = 0;
        end
        rule = regexprep(rule, ['\<' g '\>'], sprintf('%.12f', val));
    end
    exprStr = evalTree(rule);
end

function expr = evalTree(expr)
    expr = strtrim(expr);
    if isempty(expr)
        expr = '0';
        return;
    end

    % Remove unnecessary outer parentheses
    while startsWith(expr, '(') && endsWith(expr, ')')
        sub = expr(2:end-1);
        if balanced(sub)
            expr = sub;
        else
            break;
        end
    end

    % Split at top-level OR (|)
    idx = splitTopLevel(expr, '|');
    if numel(idx) > 1
        parts = cellfun(@evalTree, idx, 'UniformOutput', false);
        expr = ['(' strjoin(parts, ' + ') ')'];
        return;
    end

    % Split at top-level AND (&)
    idx = splitTopLevel(expr, '&');
    if numel(idx) > 1
        parts = cellfun(@evalTree, idx, 'UniformOutput', false);
        expr = ['min(' strjoin(parts, ',') ')'];
        return;
    end

    % No operators: return as-is
    expr = expr;
end

function isBalanced = balanced(s)
    % Check if parentheses are balanced
    level = 0;
    isBalanced = true;
    for i = 1:length(s)
        if s(i) == '(', level = level + 1; end
        if s(i) == ')', level = level - 1; end
        if level < 0
            isBalanced = false;
            return;
        end
    end
    isBalanced = level == 0;
end

function parts = splitTopLevel(expr, sep)
    % Split expression at separators that are at the top parenthesis level
    parts = {};
    current = '';
    level = 0;
    i = 1;
    while i <= length(expr)
        c = expr(i);
        if c == '(', level = level + 1; end
        if c == ')', level = level - 1; end
        if level == 0 && i + length(sep) - 1 <= length(expr) && strcmp(expr(i:i+length(sep)-1), sep)
            parts{end+1} = strtrim(current);
            current = '';
            i = i + length(sep);
            continue;
        else
            current = [current c];
            i = i + 1;
        end
    end
    if ~isempty(current)
        parts{end+1} = strtrim(current);
    end
end
