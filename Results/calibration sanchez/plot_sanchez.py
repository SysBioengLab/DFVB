import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
import numpy as np
import statistics as stats

ub_DFBA   = pd.read_csv('UB_Bootstraping_DFBA.csv')
lb_DFBA   = pd.read_csv('LB_Bootstraping_DFBA.csv')
mean_DFBA = pd.read_csv('Mean_Bootstraping_DFBA.csv')

ub_DFVB   = pd.read_csv('UB_Bootstraping_DFVB.csv')
lb_DFVB   = pd.read_csv('LB_Bootstraping_DFVB.csv')
mean_DFVB = pd.read_csv('Mean_Bootstraping_DFVB.csv')
data_exp  = pd.read_csv('Aerobic_data_sanchez_1.csv')

palette = sns.color_palette("plasma", 8)
color1  = palette[0]
color2  = sns.color_palette("plasma", 20)[-3]
color3  = palette[2]
color4  = palette[-2]

fontsize_labels = 16
fontsize_legend = 14
fontsize_ticks  = 12

variables = ['X', 'GLC', 'ETOH', 'GLY', 'CIT', 'LAC']
ylabels = ['Biomass (g/L)', 'Glucose (g/L)', 'Ethanol (g/L)', 'Glycerol (g/L)', 'Citrate (g/L)', 'Lactate (g/L)']

# Calcula el RMSE para todos los compuestos y muestra una tabla
results = []
for var, ylabel in zip(variables, ylabels):
    mask = ~data_exp[var].isna()
    x_exp = data_exp.loc[mask, 'TIME']
    y_exp = data_exp.loc[mask, var]
    interp_DFBA = interp1d(mean_DFBA['TIME'], mean_DFBA[var], kind='linear', fill_value='extrapolate')
    interp_DFVB = interp1d(mean_DFVB['TIME'], mean_DFVB[var], kind='linear', fill_value='extrapolate')
    y_model_DFBA = interp_DFBA(x_exp)
    y_model_DFVB = interp_DFVB(x_exp)
    rmse_DFBA = np.sqrt(np.mean((y_model_DFBA - y_exp) ** 2))
    rmse_DFVB = np.sqrt(np.mean((y_model_DFVB - y_exp) ** 2))
    results.append({'Compuesto': ylabel, 'RMSE DFBA': rmse_DFBA, 'RMSE DFVB': rmse_DFVB})

tabla_errores = pd.DataFrame(results)
print('RMSE DFVB vs DFBA:')
print(tabla_errores)
print('\n')

fig, axes = plt.subplots(2, 3, figsize=(18, 9), dpi=500)

for idx, (ax, var, ylabel) in enumerate(zip(axes.flat, variables, ylabels)):
    ax.fill_between(lb_DFBA['TIME'], lb_DFBA[var], ub_DFBA[var], color=color2, alpha=0.3, label='DFBA 95% CI')
    ax.fill_between(lb_DFVB['TIME'], lb_DFVB[var], ub_DFVB[var], color=color1, alpha=0.3, label='DFVB 95% CI')
    ax.plot(mean_DFBA['TIME'], mean_DFBA[var], color=color2, label='DFBA Mean', linewidth=2)
    ax.plot(mean_DFVB['TIME'], mean_DFVB[var], color=color1, label='DFVB Mean', linewidth=2)
    ax.plot(data_exp['TIME'], data_exp[var], 'o', color='black', label='Experimental Data')
    ax.set_xlim(0, 11)
    ax.set_ylabel(ylabel, fontsize=fontsize_labels)
    ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)

    # Añadir RMSE en cada panel
    rmse_dfba = tabla_errores.loc[idx, 'RMSE DFBA']
    rmse_dfvb = tabla_errores.loc[idx, 'RMSE DFVB']

    if var in ["GLC"]:
        ax.text(0.95, 0.95, f'RMSE DFVB: {rmse_dfvb:.3f}\nRMSE DFBA: {rmse_dfba:.3f}',
        transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='right',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='k'))
    else:
        ax.text(0.05, 0.95, f'RMSE DFVB: {rmse_dfvb:.3f}\nRMSE DFBA: {rmse_dfba:.3f}',
        transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='k'))

    if var in ['GLY', 'CIT', 'LAC']:
        ax.set_xlabel('Time (h)', fontsize=fontsize_labels)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=5, fontsize=fontsize_legend)
plt.tight_layout()
plt.savefig('../Figs/Calibration_Sanchez.pdf', bbox_inches='tight')

# Parameter analysis
parameter_dfba = pd.read_csv("Bootstraping_params_DFBA_Sanchez.csv")
parameter_DFVB = pd.read_csv("Bootstraping_params_DFVB_Sanchez.csv")

name_parameters = ['X_0', 'G_0', 'E_0', 'Gly_0', 'Cit_0', 'Lac_0',
                   'V_max', 'K_g', 'Kie', 'lb_02', 'Yx_gly', 'Yx_cit', 'Yx_lac']

def ci_width_over_mean(df):
    num = df.select_dtypes(include=[np.number])
    n = num.count()
    means = num.mean()
    ci_lower = num.quantile(0.025)
    ci_upper = num.quantile(0.975)
    return pd.DataFrame({
        "n": n,
        "mean": means,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper
    })

ci_dfba = ci_width_over_mean(parameter_dfba)
ci_DFVB = ci_width_over_mean(parameter_DFVB)

ci_DFVB.index = name_parameters
ci_dfba.index = name_parameters

comparacion_df = pd.DataFrame({
    'mean_dfba': ci_dfba['mean'],
    'ci_dfba': ci_dfba[['ci_lower', 'ci_upper']].apply(tuple, axis=1),
    'mean_dfvb': ci_DFVB['mean'],
    'ci_dfvb': ci_DFVB[['ci_lower', 'ci_upper']].apply(tuple, axis=1),
    'width_ratio_dfvb/dfba': (
        (ci_DFVB['ci_upper'] - ci_DFVB['ci_lower']) /
        (ci_dfba['ci_upper'] - ci_dfba['ci_lower'])
    ),
    'relative_mean_difference_dfba': (
        (ci_DFVB['mean'] - ci_dfba['mean']) / ci_dfba['mean']
    )
})
print(comparacion_df)
comparacion_df.to_csv('Parameter_Comparison_Sanchez.csv', index=True)

def bootstrap_percentile_summary(boot, theta_hat=None, alpha=0.05):
    """boot: DataFrame (B réplicas x parámetros). theta_hat: Series con el ajuste a los datos originales."""
    boot = boot.select_dtypes(include=[np.number]).dropna()   # quita réplicas fallidas
    point = theta_hat if theta_hat is not None else boot.median()
    lo, hi = boot.quantile(alpha / 2), boot.quantile(1 - alpha / 2)
    return pd.DataFrame({
        "n": boot.count(),
        "estimate": point,
        "boot_se": boot.std(ddof=1),          # ESTE es el error estándar (sin dividir por sqrt(B))
        "ci_lower": lo,
        "ci_upper": hi,
        "rel_width": (hi - lo) / point.abs(),
    })

def frac_at_bounds(boot, lb, ub, rtol=1e-3):
    """Fracción de réplicas pegadas a la cota inferior o superior (lb, ub: Series por parámetro)."""
    tol = rtol * (ub - lb)
    return ((boot - lb).le(tol) | (ub - boot).le(tol)).mean()

# Asigna nombres a las columnas ANTES de resumir (así no dependes del orden del índice)
parameter_dfba.columns = name_parameters
parameter_DFVB.columns = name_parameters

s_dfba = bootstrap_percentile_summary(parameter_dfba)   # o theta_hat=theta_hat_dfba
s_dfvb = bootstrap_percentile_summary(parameter_DFVB)

comparacion_df = pd.DataFrame({
    "est_dfba": s_dfba["estimate"], "ci_dfba": s_dfba[["ci_lower", "ci_upper"]].apply(tuple, axis=1),
    "rel_width_dfba": s_dfba["rel_width"],
    "est_dfvb": s_dfvb["estimate"], "ci_dfvb": s_dfvb[["ci_lower", "ci_upper"]].apply(tuple, axis=1),
    "rel_width_dfvb": s_dfvb["rel_width"],
    "width_ratio_dfvb/dfba": s_dfvb["rel_width"] / s_dfba["rel_width"],
})
print(comparacion_df)


