"""Replot saved network outcomes with fonts specified at manuscript width."""
from pathlib import Path
import csv
import hashlib
import itertools
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
import os
OUT = Path(os.environ['REBUILD_OUTPUT']) / 'figures'
OUT.mkdir(parents=True,exist_ok=True)
SOURCE = ROOT / 'results/bayesian_network/bayesian_network_exact.csv'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 9, 'xtick.labelsize': 9,
    'ytick.labelsize': 9, 'legend.fontsize': 9, 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
    'savefig.bbox': None})
metadata = {'CreationDate': None, 'ModDate': None}
rows = list(csv.DictReader(SOURCE.open(newline='')))
rows = [r for r in rows if r['tree_order'] == 'lexicographic']
assert len(rows) == 800
eps = [.01, .05, .10, .20]
conditions = [(n, a) for n in (100, 300) for a in (.7, 1.2)]
colors = ['#24608A', '#5796BF', '#A85A21', '#DB924E']
styles = ['-', '--', '-', '--']
markers = ['o', '^', 'D', 'o']
variables = ['cost_ratio_vs_ranked_dags', 'cost_ratio_vs_ranked_skeletons',
             'maximum_skeleton_edge_error']
titles = ['(a) Grouped vs\nindividual DAGs', '(b) Optimized vs\nranked skeletons',
          '(c) Posterior edge-\nprobability error']
ylabels = ['Complexity ratio', 'Complexity ratio',
           'Maximum absolute change']
limits = [(0, .40), (.85, 1.01), (0, .21)]
fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.15))
records = []
for panel, ax in enumerate(axes):
    ax.set_title(titles[panel], fontsize=9, pad=8)
    ax.set(xlim=(0, .21), ylim=limits[panel],
           xlabel='Omitted-mass\ntolerance', ylabel=ylabels[panel])
    ax.set_xticks(eps, ['.01', '.05', '.10', '.20'])
    ax.tick_params(length=2.5, pad=3)
    ax.grid(axis='y', color='#EAEAEA', linewidth=.7)
    ax.set_axisbelow(True)
    if panel == 1:
        ax.set_yticks([.85, .90, .95, 1.00])
        ax.axhline(1, color='#555555', ls=':', lw=.8)
    if panel == 2:
        ax.set_yticks([0, .05, .10, .15, .20])
        ax.plot([0, .21], [0, .21], color='#555555', ls=':', lw=.8)
    for i, (n, a) in enumerate(conditions):
        cells = [np.array([float(r[variables[panel]]) for r in rows
                          if int(r['sample_size']) == n and float(r['signal']) == a
                          and float(r['epsilon']) == e]) for e in eps]
        assert all(len(v) == 50 for v in cells)
        middle = np.array([(np.mean(v) if panel == 2 else np.median(v)) for v in cells])
        lo, hi = np.array([np.quantile(v, [.1, .9]) for v in cells]).T
        ax.fill_between(eps, lo, hi, color=colors[i], alpha=.11, linewidth=0)
        ax.plot(eps, middle, color=colors[i], ls=styles[i], marker=markers[i],
                markersize=2.8, linewidth=1.05, label=f'n = {n}, signal = {a:.1f}')
        records.append(dict(metric=variables[panel], sample_size=n, signal=a,
                            epsilon=eps, center=middle.tolist(), q10=lo.tolist(), q90=hi.tolist()))
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False,
           handlelength=2.2, columnspacing=2, bbox_to_anchor=(.5, .005))
fig.subplots_adjust(left=.082, right=.99, top=.81, bottom=.34, wspace=.58)
fig.savefig(OUT / 'bayesian_network_compression.pdf', metadata=metadata)
plt.close(fig)


(OUT.parent/'network_values.json').write_text(json.dumps(records,indent=2))
