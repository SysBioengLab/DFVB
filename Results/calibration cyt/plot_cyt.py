import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
import statistics as stats
import numpy as np

ub_DFBA = pd.read_csv('UB_Bootstraping_DFBA.csv')
lb_DFBA = pd.read_csv('LB_Bootstraping_DFBA.csv')
mean_DFBA = pd.read_csv('Mean_Bootstraping_DFBA.csv')

ub_DFVB = pd.read_csv('UB_Bootstraping_DFVB.csv')
lb_DFVB = pd.read_csv('LB_Bootstraping_DFVB.csv')
mean_DFVB = pd.read_csv('Mean_Bootstraping_DFVB.csv')

data_exp = pd.read_csv('Experiment_1.csv')

palette = sns.color_palette("plasma", 8)
color1  = palette[0]
color2  = sns.color_palette("plasma", 20)[-3]
color3  = palette[2]
color4  = palette[-2]

fontsize_labels = 16
fontsize_legend = 14
fontsize_ticks  = 12

variables = ['Biomass', 'Glucose', 'Fructose', 'Ethanol', 'Ammonia', 'PAN']
ylabels = ['Biomass (g/L)', 'Glucose (g/L)', 'Fructose (g/L)', 'Ethanol (g/L)', 'Ammonia (g/L)', 'PAN (g/L)']

# Calcula el RMSE para todos los compuestos y muestra una tabla
results = []
for var, ylabel in zip(variables, ylabels):
    mask = ~data_exp[var].isna()
    x_exp = data_exp.loc[mask, 'Time']
    y_exp = data_exp.loc[mask, var]
    interp_DFBA = interp1d(mean_DFBA['Time'], mean_DFBA[var], kind='linear', fill_value='extrapolate')
    interp_DFVB = interp1d(mean_DFVB['Time'], mean_DFVB[var], kind='linear', fill_value='extrapolate')
    y_model_DFBA = interp_DFBA(x_exp)
    y_model_DFVB = interp_DFVB(x_exp)
    rmse_DFBA = np.sqrt(np.mean((y_model_DFBA - y_exp) ** 2))
    rmse_DFVB = np.sqrt(np.mean((y_model_DFVB - y_exp) ** 2))
    results.append({'Compuesto': ylabel, 'RMSE DFBA': rmse_DFBA, 'RMSE DFVB': rmse_DFVB})

tabla_errores = pd.DataFrame(results)
print('RMSE DFBA vs DFVB:')
print(tabla_errores)
print("\n")


fontsize_labels = 16
fontsize_legend = 14
fontsize_ticks  = 12

fig, axes = plt.subplots(2, 3, figsize=(18, 9), dpi=500)

for idx, (ax, var, ylabel) in enumerate(zip(axes.flat, variables, ylabels)):
    ax.fill_between(lb_DFBA['Time'], lb_DFBA[var], ub_DFBA[var], color=color2, alpha=0.3, label='DFBA 95% CI')
    ax.fill_between(lb_DFVB['Time'], lb_DFVB[var], ub_DFVB[var], color=color1, alpha=0.3, label='DFVB 95% CI')
    ax.plot(mean_DFBA['Time'], mean_DFBA[var], color=color2, label='DFBA Mean', linewidth=2)
    ax.plot(mean_DFVB['Time'], mean_DFVB[var], color=color1, label='DFVB Mean', linewidth=2)
    ax.plot(data_exp['Time'], data_exp[var], 'o', color='black', label='Experimental Data')
    ax.set_xlim(0, 194)
    ax.set_ylabel(ylabel, fontsize=fontsize_labels)
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

    # Añadir RMSE en cada panel
    rmse_dfba = tabla_errores.loc[idx, 'RMSE DFBA']
    rmse_dfvb = tabla_errores.loc[idx, 'RMSE DFVB']

    if var not in ["Biomass","Ethanol"]:
        ax.text(0.95, 0.95, f'RMSE DFVB: {rmse_dfvb:.3f}\nRMSE DFBA: {rmse_dfba:.3f}',
        transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='k'))
    else:
        ax.text(0.05, 0.95, f'RMSE DFVB: {rmse_dfvb:.3f}\nRMSE DFBA: {rmse_dfba:.3f}',
        transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='k'))

    if var in ['Ethanol', 'Ammonia', 'PAN']:
        ax.set_xlabel('Time (h)', fontsize=fontsize_labels)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=5, fontsize=fontsize_legend)
plt.tight_layout()
plt.savefig('../Figs/Calibration_CyT.pdf', bbox_inches='tight')

# Parameter analysis
parameter_dfba = pd.read_csv("Bootstraping_params_DFBA_CYT.csv")
parameter_DFVB = pd.read_csv("Bootstraping_params_DFVB_CYT.csv")

name_parameters = ['X_0', 'N_0', 'G_0', 'F_0', 'E_0', 'PAN_0',
                   r'\mu_0', 'K_n0', 'K_g0', 'K_f0', 'K_ie0', 'Y_xn', 'Y_xg', 'Y_xf', 'Y_ef', r'\beta_f0', r'\beta_ie0']

def ci_width_over_mean(df):
    num = df.select_dtypes(include=[np.number])
    n = num.count()
    means = num.mean()
    stds = num.std(ddof=1)
    se = stds / np.sqrt(n)
    ci_lower = means - 1.96 * se
    ci_upper = means + 1.96 * se
    ci_width = 2 * 1.96 * se
    ratio = ci_width / means.replace(0, 1)  
    return pd.DataFrame({
        "n": n,
        "mean": means,
        "std": stds,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "ci_width_over_mean": ratio
    })

# Uso
ci_dfba = ci_width_over_mean(parameter_dfba)
ci_DFVB = ci_width_over_mean(parameter_DFVB)

ci_DFVB.index = name_parameters
ci_dfba.index = name_parameters

comparacion_df = pd.DataFrame({
    'mean_dfba': ci_dfba['mean'],
    'std_dfba': ci_dfba['std'],
    'ci_width_over_mean_dfba': ci_dfba['ci_width_over_mean'],
    'mean_dfvb': ci_DFVB['mean'],
    'std_dfvb': ci_DFVB['std'],
    'ci_width_over_mean_dfvb': ci_DFVB['ci_width_over_mean']
})
print(comparacion_df)