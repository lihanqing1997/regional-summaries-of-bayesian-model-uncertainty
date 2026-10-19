"""Opt-in recomputation of the original component bootstrap (5,000 resamples).
Not used by rebuild.py: the ordinary rebuild reads the saved intervals.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();OUT=args.output;OUT.mkdir(parents=True,exist_ok=False)
ROOT=Path(__file__).resolve().parents[1]
rng=np.random.default_rng(20261001)
d=pd.read_csv(ROOT/'results/bma_components/outcomes.csv')
d['cost_ratio']=d.cost_upper/d.exact_cost
p=d[d.family=='primary'].copy()
def boot(x):
    x=np.asarray(x,float);x=x[np.isfinite(x)]
    if not len(x):return [None,None]
    means=x[rng.integers(0,len(x),(5000,len(x)))].mean(axis=1)
    return np.quantile(means,[.025,.975]).tolist()
base=p[p.epsilon==.1].set_index('job')
pairrows=[];pairedraw=[]
for (family,method,budget,order),g in d[d.family!='primary'].groupby(['family','method','budget','tree_order']):
    m=g.merge(base[['cost_upper','verified_feasible','verified_optimal','seconds_used','calls','cost_ratio']],left_on='job',right_index=True,suffixes=('_extra','_primary'),validate='many_to_one')
    if family=='method_comparison':
        ref='extra';new='primary';label='adaptive_KL_vs_'+method
    else:ref='primary';new='extra';label=f'{family}_{budget}_{order}'
    for case,sub in [('ALL',m)]+list(m.groupby('case_id')):
        both=sub[f'verified_feasible_{ref}'] & sub[f'verified_feasible_{new}']
        z=sub[both].copy()
        reduction=1-z[f'cost_upper_{new}']/z[f'cost_upper_{ref}']
        ci=boot(reduction)
        lo=sub[f'cost_upper_{new}'].fillna(np.inf);hi=sub[f'cost_upper_{ref}'].fillna(np.inf)
        pairrows.append(dict(comparison=label,case_id=case,n=len(sub),reference_feasible=int(sub[f'verified_feasible_{ref}'].sum()),new_feasible=int(sub[f'verified_feasible_{new}'].sum()),
          reference_certified=int(sub[f'verified_optimal_{ref}'].sum()),new_certified=int(sub[f'verified_optimal_{new}'].sum()),both_feasible=int(both.sum()),
          new_better=int((lo<hi-1e-8).sum()),new_equal=int((lo==hi).sum()),new_worse=int((lo>hi+1e-8).sum()),
          mean_cost_reduction=reduction.mean(),median_cost_reduction=reduction.median(),reduction_ci_low=ci[0],reduction_ci_high=ci[1],
          reference_ratio_median=sub[f'cost_ratio_{ref}'].median(),new_ratio_median=sub[f'cost_ratio_{new}'].median(),
          median_seconds_change=(sub[f'seconds_used_{new}']-sub[f'seconds_used_{ref}']).median()))
    m['comparison']=label;pairedraw.append(m)
paired=pd.DataFrame(pairrows);paired.to_csv(OUT/'paired_comparisons.csv',index=False)
pd.concat(pairedraw).to_csv(OUT/'paired_observations.csv',index=False)
