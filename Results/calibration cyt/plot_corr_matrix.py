import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
import statistics as stats
import numpy as np



corr_matrix_DFBA = pd.read_csv('..\\..\\Calibration\\Results_Bootstraping\\Correlation_matrix_DFBA_CYT.csv')
corr_matrix_DFVB = pd.read_csv('..\\..\\Calibration\\Results_Bootstraping\\Correlation_matrix_DFVB_CYT.csv')

fontsize_labels = 16
fontsize_legend = 14
fontsize_ticks  = 12

labels = [r"$\mu_0$",r"$K_{N,0}$",r"$K_{G,0}$",r"$K_{F,0}$",r"$K_{E,0}$",r"$Y_{XN}$",r"$Y_{XG}$",r"$Y_{XF}$",r"$Y_{EF}$",r"$\beta_{G0}$",r"$\beta_{F0}$"]

plt.figure(figsize=(11,5))

plt.subplot(1,2,1)
sns.heatmap(corr_matrix_DFBA, annot=False, cmap='coolwarm', vmin=-1, vmax=1,
            xticklabels=labels, yticklabels=labels, cbar=False, square=True)
plt.title('DFBA', fontsize=fontsize_labels)

plt.subplot(1,2,2)
sns.heatmap(corr_matrix_DFVB, annot=False, cmap='coolwarm', vmin=-1, vmax=1,
            xticklabels=labels, yticklabels=labels, cbar=True, square=True)

plt.title('DFVB', fontsize=fontsize_labels)

plt.tight_layout()
plt.savefig("../Figs/Correlation_matrix_DFBA_DFVB_CYT.svg")
plt.close()

threshold = 0.90
mask_DFBA = corr_matrix_DFBA.abs() < threshold
mask_DFVB = corr_matrix_DFVB.abs() < threshold

plt.figure(figsize=(11,5))

plt.subplot(1,2,1)
sns.heatmap(corr_matrix_DFBA.where(~mask_DFBA, 0.0), annot=False, cmap='coolwarm',
            vmin=-1, vmax=1, xticklabels=labels, yticklabels=labels,
            cbar=False, square=True)
plt.title('DFBA', fontsize=fontsize_labels)

plt.subplot(1,2,2)
sns.heatmap(corr_matrix_DFVB.where(~mask_DFVB, 0.0), annot=False, cmap='coolwarm',
            vmin=-1, vmax=1, xticklabels=labels, yticklabels=labels,
            cbar=True, cbar_kws={'label': 'Correlation coefficient'}, square=True)
plt.title('DFVB', fontsize=fontsize_labels)

plt.tight_layout()
plt.savefig("../Figs/Correlation_matrix_DFBA_DFVB_CYT_abs_ge_0p95.svg")
plt.close()

