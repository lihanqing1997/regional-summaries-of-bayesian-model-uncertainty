"""Rebuild empirical tables from saved outcomes; no sampling or model fitting."""
from pathlib import Path
import csv, json, math, os
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ['REBUILD_OUTPUT'])/'tables'
OUT.mkdir(parents=True,exist_ok=True)
def save(name,rows):
    frame=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(OUT/(name+'.csv'),index=False)
    (OUT/(name+'.html')).write_text(frame.to_html(index=False,border=0),encoding='utf8')

s=pd.read_csv(ROOT/'results/bma_scaling/outcomes.csv')
methods=['ars_newton','mc3_frequency','mc3_renormalized']
def scaling(family):
    rows=[]
    for method in methods:
        g=s[(s.family==family)&(s.budget==100000)&(s.method==method)]
        rows.append(dict(method=method,attempts=len(g),feasible=int(g.feasible.sum()),
            optimum_attained=int(g.attained_optimal.sum()),
            optimum_certified=int(g.certified_optimal.sum()) if method=='ars_newton' else None,
            median_complexity_ratio=g.loc[g.feasible,'cost_ratio'].median(),
            median_seconds=g.seconds.median()))
    return rows
save('table04_main_comparison',scaling('main'))
save('table07a_stress_comparison',scaling('stress'))
rows=[]
for (family,budget,method),g in s[s.method.str.startswith('mc3')].groupby(['family','budget','method']):
    misses=g[np.isfinite(g.cost)&(g.retained_mass<.9-1e-9)]
    shortfall=.9-misses.retained_mass
    rows.append(dict(family=family,budget=budget,method=method,misses=len(misses),attempts=len(g),
        median_shortfall=shortfall.median(),maximum_shortfall=shortfall.max()))
save('table07b_mass_shortfalls',rows)

d=pd.read_csv(ROOT/'results/bma_components/outcomes.csv')
p=d[d.family=='primary']
rows=[]
for (case,epsilon),g in p.groupby(['case_id','epsilon'],sort=False):
    assert len(g)==100
    rows.append(dict(case_id=case,epsilon=epsilon,attempts=len(g),feasible=int(g.verified_feasible.sum()),certified=int(g.verified_optimal.sum())))
save('table08_component_outcomes',rows)

paired=pd.read_csv(ROOT/'results/bma_components/paired_comparisons.csv')
table9=paired[paired.case_id.isin(['stress_quartets_p20_b005','stress_quartets_p20_b020'])].copy()
assert len(table9)==8
# Recompute point estimates and paired counts, but reuse the recorded bootstrap endpoints.
base=p[p.epsilon==.1].set_index('job')
for record in table9.itertuples():
    g=d[d.case_id==record.case_id]
    if record.comparison.startswith('adaptive_KL_vs_'):
        method=record.comparison.removeprefix('adaptive_KL_vs_')
        extra=g[(g.family=='method_comparison')&(g.method==method)]
        old=extra.set_index('job').cost_upper
        new=base.loc[old.index,'cost_upper']
    else:
        budget=int(record.comparison.split('_')[1])
        extra=g[(g.family=='budget')&(g.budget==budget)]
        new=extra.set_index('job').cost_upper
        old=base.loc[new.index,'cost_upper']
    reduction=1-new/old
    assert len(reduction)==record.n==20
    assert np.isclose(reduction.mean(),record.mean_cost_reduction,rtol=0,atol=1e-14)
    assert int((new<old-1e-8).sum())==record.new_better
    assert int((new==old).sum())==record.new_equal
    assert int((new>old+1e-8).sum())==record.new_worse
save('table09_paired_components',table9[['comparison','case_id','n','mean_cost_reduction','reduction_ci_low','reduction_ci_high','new_better','new_equal','new_worse']])

t=pd.read_csv(ROOT/'results/bma_components/saved_trace_diagnostics.csv')
trace=t[t.study=='primary'].merge(p[p.epsilon==.1][['job','verified_feasible','verified_optimal','cost_lower','cost_upper','exact_cost']],on='job',validate='one_to_one')
assert len(trace)==500
assert np.allclose(trace.final_lower,trace.cost_lower)
assert np.allclose(trace.final_cost.fillna(np.inf),trace.cost_upper.fillna(np.inf))
rows=[]
for name,cases,terminals in [('Pairs',[10,11,12],59049),('Quartets',[23,24],243)]:
    g=trace[trace.case_index.isin(cases)]
    rows.append(dict(design=name,datasets=len(g),supplied_terminal_regions=terminals,
        median_exposed_terminal_regions=g.exposed_terminals.median(),
        stochastic_refinement=int((g.sampled_checkpoints>0).sum()),feasible=int(g.verified_feasible.sum()),certified=int(g.verified_optimal.sum())))
save('table10_search_diagnostics',rows)

save('table05_tecator_patterns',pd.read_csv(ROOT/'results/tecator/real_patterns.csv'))
tec=pd.read_csv(ROOT/'results/tecator/outcomes.csv')
tec['target']=1-tec.epsilon
tec['omitted_bound_rounded_up']=np.ceil(tec.omitted_bound*1e6)/1e6
save('table11_tecator',tec)
sc=pd.read_csv(ROOT/'results/sachs/comparison.csv')
save('table12_sachs',sc[(sc.method=='optimized')&sc.target.isin([.95,.99])])

sel=pd.read_csv(ROOT/'results/bma_selectors/outcomes.csv')
reference=sel[sel.method=='ars'].set_index('job')
selector=[]
for method,g in sel.groupby('method',sort=False):
    row=dict(method=method,n=len(g),feasible=int(g.feasible.sum()),certified=int(g.certified.sum()))
    if method!='ars':
        old=g.set_index('job').cost;new=reference.loc[old.index,'cost']
        row.update(mean_complexity_reduction=float((1-new/old).mean()),lower=int((new<old).sum()),equal=int((new==old).sum()),higher=int((new>old).sum()))
    selector.append(row)
save('text_selector_comparisons',selector)
summary={'primary_attempts':len(p),'primary_datasets':p.job.nunique(),
    'primary_feasible':int(p.verified_feasible.sum()),'primary_certified':int(p.verified_optimal.sum()),
    'primary_attained':int((p.verified_feasible&np.isclose(p.cost_upper,p.exact_cost,atol=1e-8)).sum()),
    'primary_nonreturns':int((~p.verified_feasible).sum()),'primary_statuses':p.status.value_counts().to_dict(),
    'bootstrap_intervals':'Copied from original saved analysis; point estimates and paired counts checked against retained outcomes',
    'new_simulations':False,'bootstrap_resampling_run':False}
assert (summary['primary_feasible'],summary['primary_attained'],summary['primary_certified'])==(8832,7241,7190)
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Rebuilt Tables 4, 5, 7-12 and selector summaries from saved results.')
