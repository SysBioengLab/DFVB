import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

data = pd.read_csv('Transcriptomics_aer_iML1515_tr40.csv')

metrics = [('Accuracy', ['Accuracy_NS', 'Accuracy_FVA']),
 ('Precision', ['Precision_NS', 'Precision_FVA']),
 ('Recall', ['Recall_NS', 'Recall_FVA']),
 ('F1-score', ['F1_NS', 'F1_FVA'])]

metric_cols   = [col for _, cols in metrics for col in cols]
summary_stats = data[metric_cols].agg(['mean', 'std']).reset_index()

rows = []
for metric_name, cols in metrics:
    for method_label, col in zip(['DFNA', 'FVA'], cols):
        m = data[col].mean()
        s = data[col].std()
        rows.append({'Metric': metric_name, 'Method': method_label, 'mean': m, 'std': s})

df_stats = pd.DataFrame(rows)

table = pd.DataFrame(index=[m for m, _ in metrics], columns=['DFNA', 'FVA'])
for metric_name, cols in metrics:
    for method_label, col in zip(['DFNA', 'FVA'], cols):
        m = data[col].mean()
        s = data[col].std()
        table.loc[metric_name, method_label] = f"{m:.3f} ± {s:.3f}"

print("Summary (mean ± std), n={}".format(len(data)))
print(table.to_string())

# Preparar datos en formato largo para el plot (usa std como barra de error)
df_long = []
for metric_name, cols in metrics:
    for method_label, col in zip(['DFVB','FVA'], cols):
        for v in data[col].values:
            df_long.append({'Metric': metric_name, 'Method': method_label, 'value': v})
df_long = pd.DataFrame(df_long)

# Gráfico de barras con error (std)
plt.figure(figsize=(9, 5), dpi=200)
# sns.set_style('whitegrid')
palette = sns.color_palette('plasma', 8)
sns.barplot(data=df_long, x='Metric', y='value', hue='Method', errorbar='sd', palette=[palette[2], palette[-2]], capsize=0.08)
plt.ylabel('Score', fontsize = 12)
plt.xlabel('')
plt.legend(title='', fontsize = 12)
plt.tight_layout()
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
# plt.grid(axis='y', alpha=0.3)
# plt.savefig('../Figs/Comparison_transcriptomics_aer_iML.pdf', dpi=300)
plt.savefig('../Figs/Comparison_transcriptomics_aer_iML.svg', dpi=300)
