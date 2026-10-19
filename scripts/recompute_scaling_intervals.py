"""Prespecified descriptive and stratified paired-bootstrap analysis.

Refuses to analyze an incomplete study. This script is separate from practical
workers and cannot affect method decisions, dataset generation, or run order.
"""
from pathlib import Path
from collections import defaultdict, Counter
import csv, hashlib, json, math
import numpy as np

import argparse
parser=argparse.ArgumentParser(description='Opt-in original scaling analysis, including 2,000 bootstrap resamples.')
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results/bma_scaling'
OUT=args.output;OUT.mkdir(parents=True,exist_ok=False)
with (RESULTS/'outcomes.csv').open(newline='') as stream:raw=list(csv.DictReader(stream))
assert len(raw)==1200
NUMERIC={'case_index','p','group_size','between','replication','budget','calls','seconds',
         'cost','lower','exact_cost','retained_mass','node_intervals','looks','cost_ratio',
         'omitted_bound','estimated_mass','preparation_seconds','sampling_interval_seconds',
         'selection_seconds','cache_setup_seconds','frontier_nodes','terminal_frontier',
         'unique_models','draws','mean_acceptance'}
BOOLEAN={'feasible','attained_optimal','reported_feasible','reported_certified','certified_optimal',
         'cost_bracket_covered','mass_bound_covered','all_node_intervals_covered','all_looks_covered'}
data=[]
for original in raw:
    row=original.copy()
    for key in NUMERIC: row[key]=float(row[key]) if row[key] not in {'NA',''} else math.nan
    for key in BOOLEAN: row[key]=None if row[key]=='NA' else row[key]=='TRUE'
    data.append(row)

def quantile(values,q):
    a=np.asarray(values,dtype=float);a=a[np.isfinite(a)]
    return float(np.quantile(a,q)) if len(a) else None
def wilson(k,n):
    z=1.959963984540054;p=k/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den
    radius=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0.0,center-radius),min(1.0,center+radius)
def summarize(rows,keys):
    grouped=defaultdict(list)
    for row in rows: grouped[tuple(row[k] for k in keys)].append(row)
    output=[]
    for group,rs in sorted(grouped.items()):
        n=len(rs);feasible=sum(r['feasible'] is True for r in rs)
        low,high=wilson(feasible,n)
        record=dict(zip(keys,group));record.update(n=n,feasible=feasible,
            feasible_rate=feasible/n,feasible_wilson_lower=low,feasible_wilson_upper=high,
            attained_optimal=sum(r['attained_optimal'] is True for r in rs),
            reported_certified=sum(r['reported_certified'] is True for r in rs),
            certified_optimal=sum(r['certified_optimal'] is True for r in rs),
            false_certificates=sum(r['reported_certified'] is True and r['certified_optimal'] is not True for r in rs),
            returned_summaries=sum(math.isfinite(r['cost']) for r in rs),
            no_summary=sum(not math.isfinite(r['cost']) for r in rs),
            runtime_available=sum(math.isfinite(r['seconds']) for r in rs),
            mass_median=quantile([r['retained_mass'] for r in rs],.5),
            mass_minimum=quantile([r['retained_mass'] for r in rs],0),
            mass_q05=quantile([r['retained_mass'] for r in rs],.05),
            mass_q25=quantile([r['retained_mass'] for r in rs],.25),
            mass_q75=quantile([r['retained_mass'] for r in rs],.75),
            undercovered_returned=sum(math.isfinite(r['cost']) and math.isfinite(r['retained_mass']) and r['retained_mass']<.9-1e-9 for r in rs),
            median_shortfall_undercovered=quantile([.9-r['retained_mass'] for r in rs
                if math.isfinite(r['cost']) and math.isfinite(r['retained_mass']) and r['retained_mass']<.9-1e-9],.5),
            maximum_shortfall=quantile([max(0,.9-r['retained_mass']) for r in rs
                if math.isfinite(r['cost']) and math.isfinite(r['retained_mass'])],1),
            mass_median_returned=quantile([r['retained_mass'] for r in rs if math.isfinite(r['cost'])],.5),
            cost_ratio_median_all_returned=quantile([r['cost_ratio'] for r in rs],.5),
            cost_ratio_q25_all_returned=quantile([r['cost_ratio'] for r in rs],.25),
            cost_ratio_q75_all_returned=quantile([r['cost_ratio'] for r in rs],.75),
            cost_ratio_median_feasible=quantile([r['cost_ratio'] for r in rs if r['feasible']],.5),
            finite_cost_brackets=sum(math.isfinite(r['cost']) and math.isfinite(r['lower']) for r in rs),
            cost_gap_median=quantile([r['cost']-r['lower'] for r in rs
                if math.isfinite(r['cost']) and math.isfinite(r['lower'])],.5),
            cost_gap_q25=quantile([r['cost']-r['lower'] for r in rs
                if math.isfinite(r['cost']) and math.isfinite(r['lower'])],.25),
            cost_gap_q75=quantile([r['cost']-r['lower'] for r in rs
                if math.isfinite(r['cost']) and math.isfinite(r['lower'])],.75),
            preparation_seconds_median=quantile([r['preparation_seconds'] for r in rs],.5),
            selection_seconds_median=quantile([r['selection_seconds'] for r in rs],.5),
            sampling_interval_seconds_median=quantile([r['sampling_interval_seconds'] for r in rs],.5),
            seconds_median=quantile([r['seconds'] for r in rs],.5),
            seconds_q25=quantile([r['seconds'] for r in rs],.25),
            seconds_q75=quantile([r['seconds'] for r in rs],.75),
            calls_median=quantile([r['calls'] for r in rs],.5),
            looks=sum(r['looks'] for r in rs),
            interval_miss_runs=sum(r['all_node_intervals_covered'] is False for r in rs),
            checkpoint_miss_runs=sum(r['all_looks_covered'] is False for r in rs),
            cost_bracket_miss_runs=sum(r['cost_bracket_covered'] is False for r in rs),
            omitted_mass_bound_miss_runs=sum(r['mass_bound_covered'] is False for r in rs))
        output.append(record)
    return output
def write_csv(name,rows):
    with (OUT/name).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
primary=data
write_csv('by_condition.csv',summarize(data,['family','p','group_size','between','budget','method']))
write_csv('by_dimension.csv',summarize(primary,['family','p','group_size','budget','method']))
pooled=summarize(data,['family','budget','method']);write_csv('pooled.csv',pooled)

rng=np.random.default_rng(1703000100)
index={(r['dataset'],r['method'],int(r['budget'])):r for r in data}
contrasts=[]
def contrast(label,family,method_a,budget_a,method_b,budget_b):
    datasets=sorted({r['dataset'] for r in data if r['family']==family})
    pairs=[(index[(ds,method_a,budget_a)],index[(ds,method_b,budget_b)]) for ds in datasets
           if (ds,method_a,budget_a) in index and (ds,method_b,budget_b) in index]
    strata=defaultdict(list)
    for i,(a,b) in enumerate(pairs): strata[a['case_index']].append(i)
    def measures(indices):
        selected=[pairs[int(i)] for i in indices]
        feasible_pairs=[(a,b) for a,b in selected if a['feasible'] and b['feasible']]
        return {
            'feasibility_difference':float(np.mean([int(a['feasible'])-int(b['feasible']) for a,b in selected])),
            'attained_optimality_difference':float(np.mean([int(a['attained_optimal'])-int(b['attained_optimal']) for a,b in selected])),
            'certification_difference':float(np.mean([int(a['certified_optimal'])-int(b['certified_optimal']) for a,b in selected])),
            'median_cost_ratio_both_feasible':quantile([a['cost']/b['cost'] for a,b in feasible_pairs],.5),
            'median_runtime_ratio_recorded_pairs':quantile([a['seconds']/b['seconds'] for a,b in selected
                if math.isfinite(a['seconds']) and math.isfinite(b['seconds']) and b['seconds']>0],.5),
            'median_score_request_ratio_recorded_pairs':quantile([a['calls']/b['calls'] for a,b in selected
                if math.isfinite(a['calls']) and math.isfinite(b['calls']) and b['calls']>0],.5)}
    point=measures(range(len(pairs)));boots={key:[] for key in point}
    for _ in range(2000):
        draw=np.concatenate([rng.choice(ids,len(ids),replace=True) for ids in strata.values()])
        values=measures(draw)
        for key,value in values.items():
            if value is not None: boots[key].append(value)
    record=dict(contrast=label,family=family,method_a=method_a,budget_a=budget_a,
        method_b=method_b,budget_b=budget_b,paired_datasets=len(pairs),
        both_feasible=sum(a['feasible'] and b['feasible'] for a,b in pairs),
        recorded_runtime_pairs=sum(math.isfinite(a['seconds']) and math.isfinite(b['seconds']) for a,b in pairs))
    for key,value in point.items():
        record[key]=value;record[key+'_lower']=quantile(boots[key],.025)
        record[key+'_upper']=quantile(boots[key],.975)
        record[key+'_bootstrap_valid']=len(boots[key])
    contrasts.append(record)
for family in ['main','stress']:
    for budget in [25000,100000]:
        for other in ['mc3_frequency','mc3_renormalized']:
            contrast('ARS versus independent MC3',family,'ars_newton',budget,other,budget)
    for method in ['ars_newton','mc3_frequency','mc3_renormalized']:
        contrast('100k versus 25k budget',family,method,100000,method,25000)
write_csv('paired_bootstrap.csv',contrasts)
statuses=Counter((r['family'],r['method'],int(r['budget']),r['status']) for r in data)
write_csv('statuses.csv',[dict(family=k[0],method=k[1],budget=k[2],status=k[3],n=v) for k,v in sorted(statuses.items())])
