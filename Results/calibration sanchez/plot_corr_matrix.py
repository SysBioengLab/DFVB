import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
import statistics as stats
import numpy as np




corr_matrix_DFBA = pd.read_csv('..\\..\\Calibration\\Results_Bootstraping\\Correlation_matrix_DFBA_Sanchez.csv')
corr_matrix_DFVB = pd.read_csv('..\\..\\Calibration\\Results_Bootstraping\\Correlation_matrix_DFVB_Sanchez.csv')

fontsize_labels = 16
fontsize_legend = 14
fontsize_ticks  = 12

labels = [r"$v_{glc}$",r"$k_G$",r"$K_{ie}$",r"$l_{b,o2}$",r"$f_{gly}$",r"$f_{cit}$",r"$f_{lac}$"]
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
plt.savefig("../Figs/Correlation_matrix_DFBA_DFVB_Sanchez.svg")
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
            cbar=True, square=True)
plt.title('DFVB', fontsize=fontsize_labels)

plt.tight_layout()
plt.savefig("../Figs/Correlation_matrix_DFBA_DFVB_Sanchez_abs_ge_0p95.svg")
plt.close()
