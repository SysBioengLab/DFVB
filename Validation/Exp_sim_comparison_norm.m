% Plot configuration
set(groot, 'DefaultAxesFontSize', 14); % Font size
set(groot, 'DefaultAxesLineWidth', 1.2); % Axes line width
set(groot, 'DefaultLineLineWidth', 1.4); % Line width
set(groot, 'DefaultFigureColor', 'w'); % White figure background
set(groot, 'DefaultAxesFontName', 'Arial'); % Default font
set(groot, 'DefaultAxesBox', 'on');
set(groot, 'DefaultAxesXColor', 'k');
set(groot, 'DefaultAxesYColor', 'k');
set(groot, 'DefaultFigurePosition', [100, 100, 900, 600]);
output_folder = 'Figures';

color1 = '#f19eec';
color2 = '#8ee1c2';
color3 = '#9190e1';
color4 = '#3c2737';


%% Experimental analysis  
% Import the model
model = load('iML1515_mod.mat').IML1515;

%%

% Import data
rxns_fva = readtable("Active_rxn_FVA_aer.csv");
rxns_ns  = readtable('Active_rxn_N_sparse_aer.csv');
exp_data = readtable('log_tpm_glucosemedia_normalized.csv', 'VariableNamingRule','preserve');

exp_only_data = table2array(exp_data(:,2:end));

%% Compare FVA and DFNA

inactivas_fva       = find(rxns_fva.Var2==0);
inativas_ns        = find(rxns_ns.Var2==0);
intersect_activas = intersect(inactivas_fva, inativas_ns);

new_dfna = setdiff(inativas_ns, inactivas_fva);


%% Calculate threshold
perc1 = zeros(size(exp_only_data, 2),1);
perc3 = zeros(size(exp_only_data, 2),1);
perc5 = zeros(size(exp_only_data, 2),1);
perc20 = zeros(size(exp_only_data, 2),1);
perc25 = zeros(size(exp_only_data, 2),1);
perc30 = zeros(size(exp_only_data, 2),1);
perc40 = zeros(size(exp_only_data, 2),1);
perc50 = zeros(size(exp_only_data, 2),1);
perc70 = zeros(size(exp_only_data, 2),1);
perc80 = zeros(size(exp_only_data, 2),1);

for k = 1:size(exp_only_data, 2)
    ras = computeRAS(model, exp_only_data(:,k));
    ras_valid = ras(~isnan(ras));

    perc1(k) = prctile(ras_valid, 1);
    perc3(k) = prctile(ras_valid, 3);
    perc5(k) = prctile(ras_valid, 5);
    perc20(k) = prctile(ras_valid, 20);
    perc25(k) = prctile(ras_valid, 25);
    perc30(k) = prctile(ras_valid, 30);
    perc40(k) = prctile(ras_valid, 40);
    perc50(k) = prctile(ras_valid, 50);
    perc70(k) = prctile(ras_valid, 70);
    perc80(k) = prctile(ras_valid, 80);


end

fprintf('\nMean of RAS percentile 1 across cultures = %.4f', mean(perc1))
fprintf('\nMean of RAS percentile 3 across cultures = %.4f', mean(perc3))
fprintf('\nMean of RAS percentile 5 across cultures = %.4f \n', mean(perc5))
fprintf('\nMean of RAS percentile 20 across cultures = %.4f \n', mean(perc20))
fprintf('\nMean of RAS percentile 25 across cultures = %.4f \n', mean(perc25))
fprintf('\nMean of RAS percentile 30 across cultures = %.4f \n', mean(perc30))
fprintf('\nMean of RAS percentile 40 across cultures = %.4f \n', mean(perc40))
fprintf('\nMean of RAS percentile 50 across cultures = %.4f \n', mean(perc50))
fprintf('\nMean of RAS percentile 40 across cultures = %.4f \n', mean(perc70))
fprintf('\nMean of RAS percentile 50 across cultures = %.4f \n', mean(perc80))


%% Normalize data
data_vis = ras(~isnan(computeRAS(model, exp_only_data(:,2))));
data_vis(data_vis>2) = 2;

figure;
hold on
histogram(data_vis, 25, 'DisplayName', 'Counts', 'EdgeColor',color3, 'FaceColor',color3)
% xline(0.2,'--', 'Color', color2, 'DisplayName','Threshold = 0.2', LineWidth=3)
xline(mean(perc40),'--', 'Color', color2, 'DisplayName','Threshold = perc 40%', LineWidth=3)
xlim([0, 1.5])
legend()
xlabel('Reaction Activity Score')
title('RAS distribution')

% exportgraphics(gcf, strcat(output_folder,'/Clasificador_RAS_norm.png'), 'Resolution',300 )

%% Analysis of all experimental conditions

% Initialize
n_cond = size(exp_data, 2) - 1; 

results = table('Size', [n_cond 15], ...
    'VariableTypes', repmat("double", 1, 15), ...
    'VariableNames', {'CondIdx','ExpTotal','ExpInactive','SimInactive_NS',...
                      'SimInactive_FVA','Accuracy_NS', 'Accuracy_FVA', 'Precision_NS', 'Precision_FVA', ...
                      'Recall_NS', 'Recall_FVA','F1_NS','F1_FVA', 'New_rxn_found', 'Diff_FVA_intersect'});

%% Export RAS
ras_output = zeros(size(rxns_ns,1)-1,size(exp_only_data, 2));
for k = 1:size(exp_only_data, 2)
    ras_output(:,k) = computeRAS(model, exp_only_data(:,k));
end

%% Classifier results

% threshold = 0.2;
threshold = mean(perc40);

for k = 1:size(exp_only_data, 2)
    % RAS
    ras = computeRAS(model, exp_only_data(:,k));
    valid_rxn_idx = ~isnan(ras);
    ras_valid = ras(valid_rxn_idx);
    
    % Simulation
    sim_ns = table2array(rxns_ns(valid_rxn_idx, 2));
    sim_fva = table2array(rxns_fva(valid_rxn_idx, 2));
    
    % Experimental > simulation
    inactive_exp_idx = ras_valid <= threshold;
    
    % Compare
    pos_inactive_exp = find(inactive_exp_idx);
    new_rxns_found   = setdiff(pos_inactive_exp, new_dfna);

    % NS predictions
    pred_ns = sim_ns(inactive_exp_idx) == 0;
    acc_ns = mean(pred_ns);
    
    % FVA predictions
    pred_fva = sim_fva(inactive_exp_idx) == 0;
    acc_fva = mean(pred_fva);

    % Simulation > experimental
    total_inact_ns = sim_ns == 0;
    total_inact_fva = sim_fva == 0;
    
    exp_from_sim_ns = ras_valid(total_inact_ns);
    exp_from_sim_fva = ras_valid(total_inact_fva);
    
    % F1-Score calculation
    TP_ns = sum(ras_valid <= threshold & sim_ns == 0);
    FP_ns = sum(ras_valid > threshold & sim_ns == 0);
    FN_ns = sum(ras_valid <= threshold & sim_ns ~= 0);
    precision_ns = TP_ns / (TP_ns + FP_ns + eps);
    recall_ns    = TP_ns / (TP_ns + FN_ns + eps);
    f1_ns        = 2 * (precision_ns * recall_ns) / (precision_ns + recall_ns + eps);
    
    TP_fva = sum(ras_valid <= threshold & sim_fva == 0);
    FP_fva = sum(ras_valid > threshold & sim_fva == 0);
    FN_fva = sum(ras_valid <= threshold & sim_fva ~= 0);
    precision_fva = TP_fva / (TP_fva + FP_fva + eps);
    recall_fva = TP_fva / (TP_fva + FN_fva + eps);
    f1_fva = 2 * (precision_fva * recall_fva) / (precision_fva + recall_fva + eps);
    
    % Save results
    results.CondIdx(k)        = k;
    results.ExpTotal(k)       = numel(ras_valid);
    results.ExpInactive(k)    = sum(inactive_exp_idx);
    results.SimInactive_NS(k) = sum(sim_ns == 0);
    results.Accuracy_NS(k)    = acc_ns;
    results.SimInactive_FVA(k)= sum(sim_fva == 0);
    results.Accuracy_FVA(k)   = acc_fva;
    results.Precision_NS(k)  = precision_ns;
    results.Precision_FVA(k) = precision_fva;
    results.Recall_NS(k)      = recall_ns;
    results.Recall_FVA(k)     = recall_fva;
    results.F1_NS(k)          = f1_ns;
    results.F1_FVA(k)         = f1_fva;
    results.New_rxn_found(k)      = length(new_rxns_found);
    results.Diff_FVA_intersect(k) = length(new_dfna);
end

%% Show results
disp(results);
% writetable(results, 'Results/IML1515/Transcriptomics_aer_iML1515_tr50.csv')

%% Compare newly identified set

hallazgo_mean = mean(100 * results.New_rxn_found./results.Diff_FVA_intersect);
hallazgo_sd   = std(100 * results.New_rxn_found./results.Diff_FVA_intersect);

%% Plots
figure("Position",[100, 100, 1800, 300])

for i = 1:4
    subplot(1,4,i)
    
    % Select data
    switch i
        case 1
            boxplot([results.Accuracy_NS(:), results.Accuracy_FVA(:)], ...
                'Labels', {'DFNA', 'FVA'}, 'Whisker', 1.5);
            title('Accuracy');
        case 4
            boxplot([results.F1_NS(:), results.F1_FVA(:)], ...
                'Labels', {'DFNA', 'FVA'}, 'Whisker', 1.5);
            title('F1-score');
        case 2
            boxplot([results.Precision_NS(:), results.Precision_FVA(:)], ...
                'Labels', {'DFNA', 'FVA'}, 'Whisker', 1.5);
            title('Precision');
        case 3
            boxplot([results.Recall_NS(:), results.Recall_FVA(:)], ...
                'Labels', {'DFNA', 'FVA'}, 'Whisker', 1.5);
            title('Recall');
    end

    % Uniform aesthetics
    ylim([0.3 0.9])
    yticks(0.3:0.1:0.9)

    % Change colors
    ax = gca;
    boxes = findobj(ax, 'Tag', 'Box');
    whiskers = findobj(ax, 'Tag', 'Whisker');
    medians = findobj(ax, 'Tag', 'Median');
    outliers = findobj(ax, 'Tag', 'Outliers');

    set(boxes, 'Color', color4, 'LineWidth', 1.5);
    set(whiskers, 'Color', color4, 'LineWidth', 1.5);
    set(medians, 'Color', color3, 'LineWidth', 1.8);
    set(outliers, 'MarkerEdgeColor', color2);
end

% exportgraphics(gcf, strcat(output_folder,'/Transcriptomics_validation_aer.png'), 'Resolution',300 )
