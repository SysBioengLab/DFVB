import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

# Read csv data
DFVB_diaux = pd.read_csv('DFVB_ecolicore_limited_conc.csv')  
dFBA_diaux = pd.read_csv('DFBA_ecolicore_limited_conc.csv')
DFVB_diaux_reducido = pd.read_csv('DFVB_ecolicore_reducido_limited_conc.csv')

DFVB_aer = pd.read_csv('DFVB_ecolicore_aer_conc.csv')  
dFBA_aer = pd.read_csv('DFBA_ecolicore_aer_conc.csv')
DFVB_aer_reducido = pd.read_csv('DFVB_ecolicore_reducido_aer_conc.csv')

DFVB_anaer = pd.read_csv('DFVB_ecolicore_anaer_conc.csv')
dFBA_anaer = pd.read_csv('DFBA_ecolicore_anaer_conc.csv')
DFVB_anaer_reducido = pd.read_csv('DFVB_ecolicore_reducido_anaer_conc.csv')

# Replace constant values of oxygen with zero
DFVB_anaer['Oxygen'] = 0
dFBA_anaer['Oxygen'] = 0

# Flux data
DFVB_flux = pd.read_csv('DFVB_ecolicore_limited_flux.csv')
DFVB_flux_reducido = pd.read_csv('DFVB_reducido_ecolicore_limited_flux.csv')
dFBA_flux = pd.read_csv('DFBA_ecolicore_limited_flux.csv')
FVA_range = pd.read_csv('FVA_range_ecolicore_limited.csv')
alphas    = pd.read_csv('DFVB_alphas_ecolicore_limited.csv')

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
xticks = [0, 0.5, 1, 1.5, 2, 2.5, 3]
plt.rcParams['font.family'] = 'Cambria'

fig, axs = plt.subplots(2, 3, figsize=(18, 8))

casos = [
    ('Aerobic', DFVB_aer, dFBA_aer),
    ('Anaerobic', DFVB_anaer, dFBA_anaer),
    ('Oxygen-limited', DFVB_diaux, dFBA_diaux)
]

for col, (titulo, df_DFVB, df_DFBA) in enumerate(casos):
    # Biomass and Glucose
    axs[0, col].plot(df_DFVB['Time'], df_DFVB['Biomass'], color=color1, label='Biomass, DFVB', linewidth=linewidth_value)
    axs[0, col].plot(df_DFBA['Time'], df_DFBA['Biomass'], '--', color=color2, label='Biomass, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[0, col].plot(df_DFVB['Time'], df_DFVB['Glucose']*MW_glc, color=color3, label='Glucose, DFVB', linewidth=linewidth_value)
    axs[0, col].plot(df_DFBA['Time'], df_DFBA['Glucose']*MW_glc, '--', color=color4, label='Glucose, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[0, col].set_title(titulo, fontsize=fontsize_title)
    axs[0, col].tick_params(axis='x', labelsize=fontsize_ticks)
    axs[0, col].set_xticks(xticks)
    axs[0, col].tick_params(axis='y', labelsize=fontsize_ticks)
    axs[0, col].grid(True, alpha=0.3)
    if col == 0:
        axs[0, col].set_ylabel('Concentration (g/L)', fontsize=fontsize_labels)

    # Oxygen and Acetate
    axs[1, col].plot(df_DFVB['Time'], df_DFVB['Oxygen']*MW_o2, color=color1, label='Oxygen, DFVB', linewidth=linewidth_value)
    axs[1, col].plot(df_DFBA['Time'], df_DFBA['Oxygen']*MW_o2, '--', color=color2, label='Oxygen, DFBA', linewidth=linewidth_seg, dashes=(5, 5))
    axs[1, col].plot(df_DFVB['Time'], df_DFVB['Acetate']*MW_ac, color=color3, label='Acetate, DFVB', linewidth=linewidth_value)
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

plt.tight_layout(pad=1.0, w_pad=2.5)
plt.savefig('../Figs/ecolicore_conc.svg', dpi=600)
plt.close()

fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=300)

# --- Subplot izquierdo: Oxygen ---
axs[0].plot(DFVB_diaux['Time'], DFVB_diaux['Oxygen']*MW_o2, color=color1, linewidth=linewidth_value, label='DFVB')
axs[0].plot(dFBA_diaux['Time'],  dFBA_diaux['Oxygen']*MW_o2,  '--', color=color2, linewidth=linewidth_seg, dashes=(5,5), label='DFBA')
axs[0].set_title('Oxygen', fontsize=fontsize_title)
axs[0].set_xlabel('Time (h)', fontsize=fontsize_labels)
axs[0].set_ylabel('Concentration (g/L)', fontsize=fontsize_labels)
axs[0].set_xticks(xticks)
axs[0].tick_params(axis='both', labelsize=fontsize_ticks)
axs[0].grid(True, alpha=0.3)
axs[0].legend(fontsize=fontsize_legend)

# --- Subplot derecho: Acetate ---
axs[1].plot(DFVB_diaux['Time'], DFVB_diaux['Acetate']*MW_ac, color=color3, linewidth=linewidth_value, label='DFVB')
axs[1].plot(dFBA_diaux['Time'],  dFBA_diaux['Acetate']*MW_ac,  '--', color=color4, linewidth=linewidth_seg, dashes=(5,5), label='DFBA')
axs[1].set_title('Acetate', fontsize=fontsize_title)
axs[1].set_xlabel('Time (h)', fontsize=fontsize_labels)
axs[1].set_xticks(xticks)
axs[1].tick_params(axis='both', labelsize=fontsize_ticks)
axs[1].grid(True, alpha=0.3)
axs[1].legend(fontsize=fontsize_legend)

plt.tight_layout()
plt.savefig('../Figs/diauxic_shift_comparison.svg', dpi=300)
# plt.show()
plt.close()


# RMSE DFVB vs DFBA

def RMSE(true, pred):
    return np.sqrt(np.mean((true - pred)**2))

# Calcular RMSE para cada variable y condición

condiciones = {
    'Aerobic': (DFVB_aer, dFBA_aer),
    'Anaerobic': (DFVB_anaer, dFBA_anaer),
    'Diauxic': (DFVB_diaux, dFBA_diaux)
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
print(rmse_pivot.to_string(float_format=lambda x: f"{x:.2e}"))
print('\n')

# # RMSE DFVB reduced vs DFBA
# condiciones = {
#     'Aerobic': (DFVB_aer_reducido, dFBA_aer),
#     'Anaerobic': (DFVB_anaer_reducido, dFBA_anaer),
#     'Diauxic': (DFVB_diaux_reducido, dFBA_diaux)
# }

# resultados = []

# for cond, (df_true, df_pred) in condiciones.items():
#     for var in df_true.columns[1:]: 
#         pred_interp = np.interp(df_true['Time'], df_pred['Time'], df_pred[var])
#         rmse = RMSE(df_true[var], pred_interp)
#         resultados.append({'Condition': cond, 'Variable': var, 'RMSE': rmse})

# rmse_df = pd.DataFrame(resultados)
# rmse_pivot = rmse_df.pivot(index='Variable', columns='Condition', values='RMSE')
# print("RMSE Reduced DFVB vs DFBA:")
# print(rmse_pivot)


