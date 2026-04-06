%% Aerobic
load("ecolicore/aerobic.mat");

DFVB_result    = data_export.Results_N_sparse;
DFVB_result_reducido = data_export.Results_N_sparse_reducido;
DFBA_result    = data_export.Results_DFBA;
FVA_range      = data_export.FVA;
FVA_DFVB       = data_export.flux_DFNA_N_sparse;
FVA_DFBA       = data_export.flux_DFBA;
FVA_DFVB_red   = data_export.flux_DFNA_N_sparse_reducido;
alpha_dynamics = data_export.Results_N_sparse_alpha;

folder = 'Results/ecolicore/';
writetable(DFBA_result, [folder,'DFBA_ecolicore_aer_conc.csv'])
writetable(DFVB_result, [folder,'DFVB_ecolicore_aer_conc.csv'])
writetable(DFVB_result_reducido, [folder,'DFVB_ecolicore_reducido_aer_conc.csv'])
writetable(FVA_DFBA, [folder,'DFBA_ecolicore_aer_flux.csv'])
writetable(FVA_DFVB, [folder,'DFVB_ecolicore_aer_flux.csv'])
writetable(FVA_DFVB_red, [folder,'DFVB_reducido_ecolicore_aer_flux.csv'])
writetable(array2table(alpha_dynamics), [folder, 'DFVB_alphas_ecolicore_aer.csv'])

index_fva = [2, 3, 4, 6];
fva_range_csv = table(FVA_range.fva_min(index_fva(1),:)', FVA_range.fva_max(index_fva(1),:)',...
    FVA_range.fva_min(index_fva(2),:)', FVA_range.fva_max(index_fva(2),:)', ...
    FVA_range.fva_min(index_fva(3),:)', FVA_range.fva_max(index_fva(3),:)',...
    FVA_range.fva_min(index_fva(4),:)', FVA_range.fva_max(index_fva(4),:)',...
    'VariableNames',{'Glucose_min', 'Glucose_max','Oxygen_min', 'Oxygen_max',...
    'Acetate_min','Acetate_max', 'Ammonia_min', 'Ammonia_max'});
writetable(fva_range_csv, [folder, 'FVA_range_ecolicore_aer.csv'])

%% Anaerobic
load("ecolicore/anaerobic.mat");

DFVB_result    = data_export.Results_N_sparse;
DFVB_result_reducido = data_export.Results_N_sparse_reducido;
DFBA_result    = data_export.Results_DFBA;
FVA_range      = data_export.FVA;
FVA_DFVB       = data_export.flux_DFNA_N_sparse;
FVA_DFBA       = data_export.flux_DFBA;
FVA_DFVB_red   = data_export.flux_DFNA_N_sparse_reducido;
alpha_dynamics = data_export.Results_N_sparse_alpha;

folder = 'Results/ecolicore/';
writetable(DFBA_result, [folder,'DFBA_ecolicore_anaer_conc.csv'])
writetable(DFVB_result, [folder,'DFVB_ecolicore_anaer_conc.csv'])
writetable(DFVB_result_reducido, [folder,'DFVB_ecolicore_reducido_anaer_conc.csv'])
writetable(FVA_DFBA, [folder,'DFBA_ecolicore_anaer_flux.csv'])
writetable(FVA_DFVB, [folder,'DFVB_ecolicore_anaer_flux.csv'])
writetable(FVA_DFVB_red, [folder,'DFVB_reducido_ecolicore_anaer_flux.csv'])
writetable(array2table(alpha_dynamics), [folder, 'DFVB_alphas_ecolicore_anaer.csv'])

index_fva = [2, 3, 4, 6];
fva_range_csv = table(FVA_range.fva_min(index_fva(1),:)', FVA_range.fva_max(index_fva(1),:)',...
    FVA_range.fva_min(index_fva(2),:)', FVA_range.fva_max(index_fva(2),:)', ...
    FVA_range.fva_min(index_fva(3),:)', FVA_range.fva_max(index_fva(3),:)',...
    FVA_range.fva_min(index_fva(4),:)', FVA_range.fva_max(index_fva(4),:)',...
    'VariableNames',{'Glucose_min', 'Glucose_max','Oxygen_min', 'Oxygen_max',...
    'Acetate_min','Acetate_max', 'Ammonia_min', 'Ammonia_max'});
writetable(fva_range_csv, [folder, 'FVA_range_ecolicore_anaer.csv'])


%% Aerobico
load("ecolicore/limited.mat");

DFVB_result    = data_export.Results_N_sparse;
DFVB_result_reducido = data_export.Results_N_sparse_reducido;
DFBA_result    = data_export.Results_DFBA;
FVA_range      = data_export.FVA;
FVA_DFVB       = data_export.flux_DFNA_N_sparse;
FVA_DFBA       = data_export.flux_DFBA;
FVA_DFVB_red   = data_export.flux_DFNA_N_sparse_reducido;
alpha_dynamics = data_export.Results_N_sparse_alpha;

folder = 'Results/ecolicore/';
writetable(DFBA_result, [folder,'DFBA_ecolicore_limited_conc.csv'])
writetable(DFVB_result, [folder,'DFVB_ecolicore_limited_conc.csv'])
writetable(DFVB_result_reducido, [folder,'DFVB_ecolicore_reducido_limited_conc.csv'])
writetable(FVA_DFBA, [folder,'DFBA_ecolicore_limited_flux.csv'])
writetable(FVA_DFVB, [folder,'DFVB_ecolicore_limited_flux.csv'])
writetable(FVA_DFVB_red, [folder,'DFVB_reducido_ecolicore_limited_flux.csv'])
writetable(array2table(alpha_dynamics), [folder, 'DFVB_alphas_ecolicore_limited.csv'])

index_fva = [2, 3, 4, 6];
fva_range_csv = table(FVA_range.fva_min(index_fva(1),:)', FVA_range.fva_max(index_fva(1),:)',...
    FVA_range.fva_min(index_fva(2),:)', FVA_range.fva_max(index_fva(2),:)', ...
    FVA_range.fva_min(index_fva(3),:)', FVA_range.fva_max(index_fva(3),:)',...
    FVA_range.fva_min(index_fva(4),:)', FVA_range.fva_max(index_fva(4),:)',...
    'VariableNames',{'Glucose_min', 'Glucose_max','Oxygen_min', 'Oxygen_max',...
    'Acetate_min','Acetate_max', 'Ammonia_min', 'Ammonia_max'});
writetable(fva_range_csv, [folder, 'FVA_range_ecolicore_limited.csv'])
