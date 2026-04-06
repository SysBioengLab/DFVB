%% Aerobic
load("Yeast9/aerobic.mat");

DFVB_result    = data_export.Results_N_sparse;
DFVB_result_reducido = data_export.Results_N_sparse_reducido;
DFBA_result    = data_export.Results_DFBA;

folder = 'Results/Yeast9/';
writetable(DFBA_result, [folder,'DFBA_Yeast9_aer_conc.csv'])
writetable(DFVB_result, [folder,'DFVB_Yeast9_aer_conc.csv'])
writetable(DFVB_result_reducido, [folder,'DFVB_Yeast9_reducido_aer_conc.csv'])


%% Anaerobic
load("Yeast9/anaerobic.mat");

DFVB_result    = data_export.Results_N_sparse;
DFVB_result_reducido = data_export.Results_N_sparse_reducido;
DFBA_result    = data_export.Results_DFBA;

folder = 'Results/Yeast9/';
writetable(DFBA_result, [folder,'DFBA_Yeast9_anaer_conc.csv'])
writetable(DFVB_result, [folder,'DFVB_Yeast9_anaer_conc.csv'])
writetable(DFVB_result_reducido, [folder,'DFVB_Yeast9_reducido_anaer_conc.csv'])

%% Aerobico
load("Yeast9/limited.mat");

DFVB_result    = data_export.Results_N_sparse;
DFVB_result_reducido = data_export.Results_N_sparse_reducido;
DFBA_result    = data_export.Results_DFBA;


folder = 'Results/Yeast9/';
writetable(DFBA_result, [folder,'DFBA_Yeast9_limited_conc.csv'])
writetable(DFVB_result, [folder,'DFVB_Yeast9_limited_conc.csv'])
writetable(DFVB_result_reducido, [folder,'DFVB_Yeast9_reducido_limited_conc.csv'])
