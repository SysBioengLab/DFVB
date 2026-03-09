import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

# Read csv data
dFNA_diaux = pd.read_csv('DFNA_ecolicore_diaux_conc.csv')  
dFBA_diaux = pd.read_csv('DFBA_ecolicore_diaux_conc.csv')
dFNA_diaux_reducido = pd.read_csv('DFNA_reducido_ecolicore_diaux_conc.csv')

dFNA_aer = pd.read_csv('DFNA_ecolicore_aer_conc.csv')  
dFBA_aer = pd.read_csv('DFBA_ecolicore_aer_conc.csv')
dFNA_aer_reducido = pd.read_csv('DFNA_reducido_ecolicore_aer_conc.csv')

dFNA_anaer = pd.read_csv('DFNA_ecolicore_anaer_conc.csv')
dFBA_anaer = pd.read_csv('DFBA_ecolicore_anaer_conc.csv')
dFNA_anaer_reducido = pd.read_csv('DFNA_reducido_ecolicore_anaer_conc.csv')

# Replace constant values of oxygen with zero
dFNA_anaer['Oxygen'] = 0
dFBA_anaer['Oxygen'] = 0

# # Flux data
dFNA_flux = pd.read_csv('DFNA_ecolicore_diaux_flux.csv')
dFNA_flux_reducido = pd.read_csv('DFNA_reducido_ecolicore_diaux_flux.csv')
dFBA_flux = pd.read_csv('DFBA_ecolicore_diaux_flux.csv')
FVA_range = pd.read_csv('FVA_range_ecolicore_diaux.csv')
alphas    = pd.read_csv('DFNA_alphas_ecolicore_diaux.csv')

MW_glc  = 180.16/1000  # g/mmol
MW_ac   = 60.05/1000    # g/mmol
MW_etoh = 46.07/1000  # g/mmol
MW_o2   = 32.00/1000    # g/mmol

# Asignar colores específicos a cada variable
palette = sns.color_palette("plasma", 8)
color1  = palette[0]
color2  = sns.color_palette("plasma", 20)[-1]
color3  = palette[2]
color4  = palette[-2]
linewidth_value = 4
linewidth_seg   = 4
fontsize_legend = 16
fontsize_labels = 16
fontsize_ticks  = 16
fontsize_title  = 16
xticks = [0, 0.5, 1, 1.5, 2]
plt.rcParams['font.family'] = 'Cambria'

fig, axs = plt.subplots(2, 3, figsize=(18, 8), dpi=300)

casos = [
    ('Aerobic', dFNA_aer, dFBA_aer),
    ('Anaerobic', dFNA_anaer, dFBA_anaer),
    ('Oxygen-limited', dFNA_diaux, dFBA_diaux)
]

for col, (titulo, df_DFNA, df_DFBA) in enumerate(casos):
    # Biomass and Glucose
    axs[0, col].plot(df_DFNA['Time'], df_DFNA['Biomass'], color=color1, label='Biomass, DFVB', linewidth=linewidth_value)
    axs[0, col].plot(df_DFBA['Time'], df_DFBA['Biomass'], '--', color=color2, label='Biomass, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[0, col].plot(df_DFNA['Time'], df_DFNA['Glucose']*MW_glc, color=color3, label='Glucose, DFVB', linewidth=linewidth_value)
    axs[0, col].plot(df_DFBA['Time'], df_DFBA['Glucose']*MW_glc, '--', color=color4, label='Glucose, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[0, col].set_title(titulo, fontsize=fontsize_title)
    axs[0, col].tick_params(axis='x', labelsize=fontsize_ticks)
    axs[0, col].set_xticks(xticks)
    axs[0, col].tick_params(axis='y', labelsize=fontsize_ticks)
    axs[0, col].grid(True, alpha=0.3)
    if col == 0:
        axs[0, col].set_ylabel('Concentration (g/L)', fontsize=fontsize_labels)

    # Oxygen and Acetate
    axs[1, col].plot(df_DFNA['Time'], df_DFNA['Oxygen']*MW_o2, color=color1, label='Oxygen, DFVB', linewidth=linewidth_value)
    axs[1, col].plot(df_DFBA['Time'], df_DFBA['Oxygen']*MW_o2, '--', color=color2, label='Oxygen, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[1, col].plot(df_DFNA['Time'], df_DFNA['Acetate']*MW_ac, color=color3, label='Acetate, DFVB', linewidth=linewidth_value)
    axs[1, col].plot(df_DFBA['Time'], df_DFBA['Acetate']*MW_ac, '--', color=color4, label='Acetate, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[1, col].set_xlabel('Time (h)', fontsize=fontsize_labels)
    axs[1, col].tick_params(axis='x', labelsize=fontsize_ticks)
    axs[1, col].tick_params(axis='y', labelsize=fontsize_ticks)
    axs[1, col].set_xticks(xticks)
    axs[1, col].grid(True, alpha=0.3)
    if col == 0:
        axs[1, col].set_ylabel('Concentration (g/L)', fontsize=fontsize_labels)

custom_lines = [
    Line2D([0], [0], color=color1, lw=linewidth_value, label='DFVB'),
    Line2D([0], [0], color=color2, lw=linewidth_seg, label='DFBA'),
    Line2D([0], [0], color=color3, lw=linewidth_value, label='DFVB'),
    Line2D([0], [0], color=color4, lw=linewidth_seg, label='DFBA')
]

legend_names = [['Biomass, DFVB', 'Biomass, DFBA', 'Glucose, DFVB', 'Glucose, DFBA'],
                     ['Oxygen, DFVB', 'Oxygen, DFBA', 'Acetate, DFVB', 'Acetate, DFBA']]
for row in range(2):
    axs[row, 2].legend(
        handles=custom_lines,
        loc=(1.05, 0.3),
        labels=legend_names[row],
        ncol=1, fontsize=fontsize_legend
    )

plt.tight_layout()
plt.savefig('../Figs/ecolicore_conc.svg', bbox_inches='tight', dpi=600)


# Diauxic alphas plot
palette = sns.color_palette("gray", alphas.shape[1]+5)[:-1]
plasma  = sns.color_palette("plasma", 8)
plt.figure(figsize=(7, 5), dpi=300)
for i, col in enumerate(alphas.columns):
    plt.plot(dFNA_diaux['Time'], alphas[col], color=palette[i], alpha=1.0, linewidth=2)
plt.xlabel('Time (h)', fontsize=fontsize_labels)
plt.ylabel(r'Activation level ($\alpha$)', fontsize=fontsize_labels)
plt.axvline(0.84, color=plasma[2], linestyle='--', linewidth=2, label='Carbon source switch')
plt.axvline(1.615, color=plasma[-2], linestyle='--', linewidth=2, label='Carbon source exhausted')
plt.legend(fontsize=fontsize_legend, framealpha=1)
plt.xlim([0, 2])
plt.ylim(-15, 22)
plt.tick_params(axis='x', labelsize=fontsize_ticks)
plt.tick_params(axis='y', labelsize=fontsize_ticks)
plt.tight_layout()
plt.savefig('../Figs/ecolicore_diaux_alphas.pdf', dpi=300)

# RMSE DFNA vs DFBA

def RMSE(true, pred):
    return np.sqrt(np.mean((true - pred)**2))

# Calcular RMSE para cada variable y condición

condiciones = {
    'Aerobic': (dFNA_aer, dFBA_aer),
    'Anaerobic': (dFNA_anaer, dFBA_anaer),
    'Diauxic': (dFNA_diaux, dFBA_diaux)
}

resultados = []

for cond, (df_true, df_pred) in condiciones.items():
    for var in df_true.columns[1:]:  
        pred_interp = np.interp(df_true['Time'], df_pred['Time'], df_pred[var])
        rmse = RMSE(df_true[var], pred_interp)
        resultados.append({'Condition': cond, 'Variable': var, 'RMSE': rmse})

rmse_df = pd.DataFrame(resultados)
rmse_pivot = rmse_df.pivot(index='Variable', columns='Condition', values='RMSE')
print("RMSE DFVB vs DFBA:")
print(rmse_pivot)
print('\n')

# RMSE DFNA reduced vs DFBA
condiciones = {
    'Aerobic': (dFNA_aer_reducido, dFBA_aer),
    'Anaerobic': (dFNA_anaer_reducido, dFBA_anaer),
    'Diauxic': (dFNA_diaux_reducido, dFBA_diaux)
}

resultados = []

for cond, (df_true, df_pred) in condiciones.items():
    for var in df_true.columns[1:]: 
        pred_interp = np.interp(df_true['Time'], df_pred['Time'], df_pred[var])
        rmse = RMSE(df_true[var], pred_interp)
        resultados.append({'Condition': cond, 'Variable': var, 'RMSE': rmse})

rmse_df = pd.DataFrame(resultados)
rmse_pivot = rmse_df.pivot(index='Variable', columns='Condition', values='RMSE')
print("RMSE Reduced DFNA vs DFBA:")
print(rmse_pivot)

# Error before carbon source switch in diauxic case
mask = dFNA_diaux_reducido['Time'] < 3
for var in dFNA_diaux_reducido.columns[1:]:  
    pred_interp = np.interp(dFNA_diaux_reducido.loc[mask, 'Time'], dFBA_diaux['Time'], dFBA_diaux[var])
    rmse = RMSE(dFNA_diaux_reducido.loc[mask, var], pred_interp)
    print(f'RMSE for {var} in reduced diauxic case (Time < 0.8 h): {rmse}')