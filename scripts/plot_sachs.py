from pathlib import Path
import os,json,itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ['REBUILD_OUTPUT'])/'figures'
OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':None})
metadata={'CreationDate':None,'ModDate':None}
selections = json.loads((ROOT / 'results/sachs/sachs_primary/selections.json').read_text())
comparison = json.loads((ROOT / 'results/sachs/sachs_primary/comparison.json').read_text())
names = np.load(ROOT / 'data/sachs/processed/sachs.npz')['names']
pairs = list(itertools.combinations(range(11), 2))
angles = np.linspace(0, 2*np.pi, 11, endpoint=False)
pos = np.c_[np.cos(angles), np.sin(angles)]
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.85))
sachs_records = []
for index, (ax, target) in enumerate(zip(axes, [.95, .99])):
    skeletons = selections[f'{target}_optimized']
    edges = []
    for e, (i, j) in enumerate(pairs):
        present = [bool(int(s, 16) & (1 << e)) for s in skeletons]
        if not any(present):
            continue
        shared = all(present)
        ax.plot(pos[[i, j], 0], pos[[i, j], 1], color='#14667c' if shared else '#bf7534',
                linewidth=1.3, linestyle='-' if shared else '--', zorder=1)
        edges.append(dict(nodes=[str(names[i]), str(names[j])], shared=shared))
    ax.scatter(pos[:, 0], pos[:, 1], s=700, color='#e2edf0',
               edgecolors='#14667c', linewidths=.8, zorder=3)
    for name, (x, y) in zip(names, pos):
        ax.text(x, y, str(name), ha='center', va='center', fontsize=9, zorder=4)
    row = next(r for r in comparison if r['target'] == target and r['method'] == 'optimized')
    count = len(skeletons)
    ax.set_title(f'({"ab"[index]}) Target {target:.0%}: {count} skeleton{"s" if count != 1 else ""}\n'
                 f'Retained mass {row["retained_mass"]:.5f}', fontsize=10, pad=8)
    ax.set(xlim=(-1.23, 1.23), ylim=(-1.23, 1.23), aspect='equal')
    ax.axis('off')
    sachs_records.append(dict(target=target, skeletons=skeletons, edges=edges,
                              retained_mass=row['retained_mass']))
fig.subplots_adjust(left=.015, right=.985, bottom=.025, top=.78, wspace=.12)
fig.savefig(OUT / 'sachs_real_skeletons.pdf', metadata=metadata)
plt.close(fig)

(OUT.parent/'sachs_values.json').write_text(json.dumps(sachs_records,indent=2))
