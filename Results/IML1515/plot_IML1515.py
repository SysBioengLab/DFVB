import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

# Read csv data
dFNA_diaux = pd.read_csv('DFNA_iML_diaux_conc.csv')  
dFBA_diaux = pd.read_csv('DFBA_iML_diaux_conc.csv')
dFNA_diaux_reducido = pd.read_csv('DFNA_iML_reducido_diaux_conc.csv')

dFNA_aer = pd.read_csv('DFNA_iML_aer_conc.csv')  
dFBA_aer = pd.read_csv('DFBA_iML_aer_conc.csv')
dFNA_aer_reducido = pd.read_csv('DFNA_iML_reducido_aer_conc.csv')

dFNA_anaer = pd.read_csv('DFNA_iML_anaer_conc.csv')
dFBA_anaer = pd.read_csv('DFBA_iML_anaer_conc.csv')
dFNA_anaer_reducido = pd.read_csv('DFNA_iML_reducido_anaer_conc.csv')


palette = sns.color_palette("plasma", 8)
color1  = palette[0]
color2  = sns.color_palette("plasma", 20)[-1]
color3  = palette[2]
color4  = palette[-2]
linewidth_value = 4
linewidth_seg   = 4
fontsize_legend = 14
fontsize_labels = 14
fontsize_ticks  = 14
fontsize_title  = 14
xticks = [0, 0.5, 1, 1.5, 2]


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
    axs[0, col].plot(df_DFNA['Time'], df_DFNA['Glucose'], color=color3, label='Glucose, DFVB', linewidth=linewidth_value)
    axs[0, col].plot(df_DFBA['Time'], df_DFBA['Glucose'], '--', color=color4, label='Glucose, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[0, col].set_title(titulo, fontsize=fontsize_title)
    axs[0, col].tick_params(axis='x', labelsize=fontsize_ticks)
    axs[0, col].set_xticks(xticks)
    axs[0, col].tick_params(axis='y', labelsize=fontsize_ticks)
    axs[0, col].yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))
    axs[0, col].grid(True, alpha=0.3)
    if col == 0:
        axs[0, col].set_ylabel('Concentration (mM)', fontsize=fontsize_labels)

    # Oxygen and Acetate
    axs[1, col].plot(df_DFNA['Time'], df_DFNA['Oxygen'], color=color1, label='Oxygen, DFVB', linewidth=linewidth_value)
    axs[1, col].plot(df_DFBA['Time'], df_DFBA['Oxygen'], '--', color=color2, label='Oxygen, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[1, col].plot(df_DFNA['Time'], df_DFNA['Acetate'], color=color3, label='Acetate, DFVB', linewidth=linewidth_value)
    axs[1, col].plot(df_DFBA['Time'], df_DFBA['Acetate'], '--', color=color4, label='Acetate, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[1, col].set_xlabel('Time (h)', fontsize=fontsize_labels)
    axs[1, col].tick_params(axis='x', labelsize=fontsize_ticks)
    axs[1, col].tick_params(axis='y', labelsize=fontsize_ticks)
    axs[1, col].set_xticks(xticks)
    axs[1, col].yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))
    axs[1, col].grid(True, alpha=0.3)
    if col == 0:
        axs[1, col].set_ylabel('Concentration (mM)', fontsize=fontsize_labels)

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
plt.savefig('../Figs/ecolicore_concentrations_comparison.svg', bbox_inches='tight', dpi=600)


def RMSE(true, pred):
    return np.sqrt(np.mean((true - pred)**2))


# RMSE for DFNA vs DFBA
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
rmse_df    = pd.DataFrame(resultados)
rmse_pivot = rmse_df.pivot(index='Variable', columns='Condition', values='RMSE')
print('RMSE Results DFNA vs DFBA:')
print(rmse_pivot)
print('\n')

# RMSE for Reduced DFNA vs Full DFNA

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

# Crear DataFrame de resultados
rmse_df = pd.DataFrame(resultados)
rmse_pivot = rmse_df.pivot(index='Variable', columns='Condition', values='RMSE')
print('RMSE Results Reduced DFNA vs DFNA:')
print(rmse_pivot)