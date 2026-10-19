"""Re-layout existing audited outcomes without changing data or running searches."""
from pathlib import Path
import csv, json, math, hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/bma_scaling/outcomes.csv'
import os
OUT = Path(os.environ['REBUILD_OUTPUT']) / 'figures'
OUT.mkdir(parents=True,exist_ok=True)
data = list(csv.DictReader(SOURCE.open(newline='')))
assert len(data) == 1200
plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.bbox': None})
METHODS = ['ars_newton', 'mc3_frequency', 'mc3_renormalized']
LABELS = ['ARS', 'MC3 frequency', 'MC3 visited-score']
COLORS = ['#1764ab', '#d46a18', '#22845a']
records = []

def value(row, key):
    try:
        return float(row[key])
    except (ValueError, TypeError):
        return math.nan

def rows_for(method, budget, group_size, p):
    rows = [r for r in data if r['family'] == 'main' and r['method'] == method
            and int(r['budget']) == budget and int(r['p']) == p and int(r['group_size']) == group_size]
    assert len(rows) == 20
    return rows

def draw(ax, x, rows, key, color, label, context):
    vals = np.array([value(r, key) for r in rows])
    vals = vals[np.isfinite(vals)]
    if len(vals):
        low, mid, high = np.quantile(vals, [.25, .5, .75])
        ax.errorbar(x, mid, yerr=[[mid-low], [high-mid]], fmt='o', color=color,
                    markersize=4, capsize=2, label=label)
        records.append(dict(**context, metric=key, n=len(vals), q25=low, median=mid, q75=high))

for budget in (100000,):
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.75), sharex=True)
    for row_index, group_size in enumerate((2, 4)):
        for method_index, (method, label, color) in enumerate(zip(METHODS, LABELS, COLORS)):
            for p in (12, 16, 20, 24):
                rows = rows_for(method, budget, group_size, p)
                feasible = [r for r in rows if r['feasible'] == 'TRUE']
                returned = [r for r in rows if math.isfinite(value(r, 'cost'))]
                x = p + (method_index-1)*.55
                for column, key, chosen in ((0, 'retained_mass', returned), (1, 'cost_ratio', feasible)):
                    draw(axes[row_index, column], x, chosen, key, color, label if p == 12 else None,
                         dict(budget=budget, group_size=group_size, p=p, method=method))
                axes[row_index, 1].text(x, .215-.085*method_index, f'{len(feasible)}/20',
                    transform=axes[row_index, 1].get_xaxis_transform(), ha='center', va='center', fontsize=8, color=color)
        for column in range(2):
            ax = axes[row_index, column]
            ax.set_xticks([12, 16, 20, 24])
            ax.set_xlim(10, 26)
            ax.grid(axis='y', alpha=.16)
            if row_index:
                ax.set_xlabel('Predictors')
        axes[row_index, 0].set_ylabel(('Pairs' if group_size == 2 else 'Quartets') + '\nRetained mass\n(returned summaries)')
        axes[row_index, 0].axhline(.9, color='black', ls=':', lw=.8)
        axes[row_index, 1].axhline(1, color='black', ls=':', lw=.8)
        axes[row_index, 1].set_ylabel('Complexity / exact optimum\n(feasible outputs)')
        axes[row_index, 1].margins(y=.65)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=3, frameon=False)
    fig.tight_layout(pad=.55, h_pad=1.0, w_pad=1.15, rect=(0, 0, 1, .93))
    fig.savefig(OUT / f'scaling_main_{budget}.pdf', metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)


(OUT.parent/'scaling_values.json').write_text(json.dumps(records,indent=2))
