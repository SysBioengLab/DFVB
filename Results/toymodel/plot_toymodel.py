import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import matplotlib as mpl

# Read csv data
dFNA = pd.read_csv('DFNA_results.csv')  
dFBA = pd.read_csv('DFBA_results.csv')    
dFNA_red = pd.read_csv('DFNA_red_results.csv')   


palette = sns.color_palette("plasma", 8)
color1  = palette[0]
color2  = sns.color_palette("plasma", 20)[-1]
color3  = palette[2]
color4  = palette[-2]
plt.rcParams['font.family'] = 'Cambria'

linewidth_value = 4
linewidth_seg   = 4
fontsize_legend = 14
fontsize_labels = 16
fontsize_ticks  = 16


fig, axs = plt.subplots(2, 3, figsize=(18, 8), dpi=300)


custom_lines = [
    Line2D([0], [0], color=color1, lw=linewidth_value, label='DFVB'),
    Line2D([0], [0], color=color3, lw=linewidth_value, label='DFVB'),
    Line2D([0], [0], color='k', lw=linewidth_seg, alpha = 0.4, label='DFVB'),
    Line2D([0], [0], color=color2, lw=linewidth_seg, label='DFBA'),
    Line2D([0], [0], color=color4, lw=linewidth_seg, label='DFBA'),
    Line2D([0], [0], color='k', lw=linewidth_seg, alpha = 0.2,label='DFBA')
]

axs[0,0].plot(dFNA['Time'], dFNA['X'], color=color1, label='X, DFVB', linewidth=linewidth_value)
axs[0,0].plot(dFBA['Time'], dFBA['X'], '--', color=color2, label='DFBA', linewidth=linewidth_seg, dashes=(5, 5))
axs[0,0].plot(dFNA['Time'], dFNA['S'], color=color3, label='S, DFVB', linewidth=linewidth_value)
axs[0,0].plot(dFBA['Time'], dFBA['S'], '--', color=color4, label='DFBA', linewidth=linewidth_seg, dashes=(5, 5))
axs[0,0].plot(dFBA['Time'], dFBA['P'], '--', color='k', alpha=0.4, label='DFBA', linewidth=linewidth_seg, dashes=(5, 5), zorder=0)
axs[0,0].plot(dFNA['Time'], dFNA['P'], color='k',alpha=0.2, label='P, DFVB', linewidth=linewidth_value, zorder=0)

axs[0,0].legend(custom_lines, ['X, DFVB', 'S, DFVB','P, DFVB', 'DFBA', 'DFBA', 'DFBA'], ncol=2, loc='best', columnspacing=0.7, fontsize=fontsize_legend)
axs[0,0].set_xlabel('Time (h)', fontsize = fontsize_labels)
axs[0,0].set_ylabel('Concentration (g/L)', fontsize = fontsize_labels)
axs[0,0].tick_params(axis='x', labelsize=fontsize_ticks)
axs[0,0].tick_params(axis='y', labelsize=fontsize_ticks)
axs[0,0].yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x)}"))
axs[0,0].grid(True, alpha=0.3)

axs[0, 1].plot(dFNA['Time'], dFNA['Alpha1'], color=color3, label=r'$\alpha_1$', linewidth=linewidth_value)
axs[0, 1].plot(dFNA['Time'], dFNA['Alpha2'], color=color4, label=r'$\alpha_2$', linewidth=linewidth_value)
axs[0, 1].legend(fontsize=fontsize_legend)
axs[0, 1].set_ylabel('Null-space basis weight', fontsize=fontsize_labels)
axs[0, 1].set_xlabel('Time (h)', fontsize=fontsize_labels)
axs[0, 1].grid(True, alpha=0.3)
axs[0, 1].tick_params(axis='x', labelsize=fontsize_ticks)
axs[0, 1].tick_params(axis='y', labelsize=fontsize_ticks)
legend_names = [['X', 'x', 's', 's', 'P', 'p'],
                     ['X', 'x', 's', 's', 'P', 'p'],]
for row in range(2):
    axs[row, 2].legend(
        handles=custom_lines,
        loc=(1.05, 0.3),
        labels=legend_names[row],
        ncol=1, fontsize=fontsize_legend
    )

plt.tight_layout()
plt.savefig('../Figs/toymodel_conc.svg', dpi=300)

# Calcualte and print RMSE for each variable
print("RMSE for X:", np.sqrt(np.mean((dFBA.X - dFNA.X)**2)))
print("RMSE for S:", np.sqrt(np.mean((dFBA.S - dFNA.S)**2)))
print("RMSE for P:", np.sqrt(np.mean((dFBA.P - dFNA.P)**2)))

# Calculate and print RMSE for each variable
print("RMSE for X (DFNA reduced):", np.sqrt(np.mean((dFBA.X - dFNA_red.X)**2)))
print("RMSE for S (DFNA reduced):", np.sqrt(np.mean((dFBA.S - dFNA_red.S)**2)))
print("RMSE for P (DFNA reduced):", np.sqrt(np.mean((dFBA.P - dFNA_red.P)**2)))