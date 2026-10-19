"""Read-only release checks. Never imports an experiment or runs a search."""
from pathlib import Path
import argparse, ast, csv, hashlib, itertools, json, math

ROOT=Path(__file__).resolve().parents[1]
def read(name):
    with (ROOT/name).open(newline='') as stream:return list(csv.DictReader(stream))
def truth(x):return x.lower()=='true'
def number(x):return math.nan if x in ('NA','') else float(x)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def verify():
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    for name,expected in manifest['files'].items():
        assert sha(ROOT/name)==expected,'Changed or missing file: '+name
    for path in ROOT.rglob('*.py'):
        if 'build' not in path.relative_to(ROOT).parts:ast.parse(path.read_text(),filename=str(path))
    scales=read('results/bma_scaling/outcomes.csv');jobs=read('experiments/bma_scaling/scaling_config/jobs.csv')
    datasets=read('experiments/bma_scaling/scaling_config/datasets.csv')
    assert len(scales)==1200 and len(jobs)==800 and len(datasets)==200
    assert len({r['dataset'] for r in scales})==200
    expected={(j['job'],m) for j in jobs for m in (['mc3_frequency','mc3_renormalized'] if j['method']=='mc3' else [j['method']])}
    assert len({(r['job'],r['method']) for r in scales})==1200
    assert {(r['job'],r['method']) for r in scales}==expected
    for budget in ['25000','100000']:
        main=[r for r in scales if r['family']=='main' and r['method']=='ars_newton' and r['budget']==budget]
        assert len(main)==160 and sum(truth(r['feasible']) for r in main)==148
        assert sum(not math.isfinite(number(r['cost'])) for r in main)==12
    components=read('results/bma_components/outcomes.csv')
    conditions=read('experiments/bma_components/config/conditions.csv')
    assert len(components)==len(conditions)==10400
    assert len({r['condition'] for r in components})==10400
    assert {r['condition'] for r in components}=={r['condition'] for r in conditions}
    assert {r['family'] for r in components}=={'primary','budget','method_comparison','tree_order'}
    primary=[r for r in components if r['family']=='primary']
    assert len(primary)==9600 and len({r['job'] for r in primary})==2400
    assert sum(truth(r['verified_feasible']) for r in primary)==8832
    assert sum(truth(r['verified_optimal']) for r in primary)==7190
    source={r['condition']:r for r in components}
    fig=read('results/bma_components/figure3_data.csv')
    assert len(fig)==200 and len({r['job'] for r in fig})==40
    for r in fig:
        original=source[r['condition']]
        for key in ['job','case_id','method','budget','cost_upper','exact_cost']:
            assert r[key]==original[key],(r['condition'],key)
        assert truth(original['verified_feasible'])
        assert math.isclose(float(r['ratio']),float(original['cost_upper'])/float(original['exact_cost']),rel_tol=1e-14)
    selectors=read('results/bma_selectors/outcomes.csv');sj=read('experiments/bma_selectors/jobs.csv')
    assert len(selectors)==144 and len(sj)==48
    assert {(r['job'],r['method']) for r in selectors}==set(itertools.product([r['job'] for r in sj],['ars','mass','mass_cost']))
    assert all(r['kind']=='simulation' for r in selectors)
    net=read('results/bayesian_network/bayesian_network_exact.csv')
    key=lambda r:(int(r['sample_size']),float(r['signal']),int(r['replication']),float(r['epsilon']),r['tree_order'])
    assert len(net)==1600
    assert {key(r) for r in net}==set(itertools.product([100,300],[.7,1.2],range(1,51),[.01,.05,.1,.2],['lexicographic','reverse']))
    for r in net:
        assert r['status']=='bounds-met' and int(r['n_dags'])==29281
        assert float(r['selected_cost'])==float(r['exact_optimal_cost'])==float(r['cost_lower'])==float(r['cost_upper'])
        assert float(r['true_omitted_mass'])<=float(r['epsilon'])+1e-12
        assert float(r['maximum_skeleton_edge_error'])<=float(r['true_omitted_mass'])+1e-12
    lookup={key(r):r for r in net}
    assert all(r['selected_cost']==lookup[(*key(r)[:4],'reverse')]['selected_cost'] for r in net if r['tree_order']=='lexicographic')
    tec=read('results/tecator/outcomes.csv')
    assert len(tec)==6 and all(truth(r['feasible']) for r in tec)
    assert sum(truth(r['certified']) for r in tec)==4
    for r in tec:assert (ROOT/'results/tecator/records'/f"{r['dataset']}_eps_{r['epsilon']}.rds").is_file()
    sachs=read('results/sachs/comparison.csv')
    reported=[r for r in sachs if r['method']=='optimized' and float(r['target']) in [.95,.99]]
    assert len(reported)==6
    for spec in ['primary','ess_10','absence_060']:
        folder=ROOT/'results/sachs'/('sachs_'+spec)
        records=json.loads((folder/'comparison.json').read_text())
        for r in records:
            if r['method']=='optimized':
                assert r['selection_solver']['optimality_proven']
                assert r['retained_mass']>=r['target']-1e-8
                ranked=next(v for v in records if v['target']==r['target'] and v['method']=='ranked')
                assert ranked['cost']==r['cost']
    return {'status':'PASS','hashed_files':len(manifest['files']),
        'scaling_outcomes':1200,'component_outcomes':10400,'selector_outcomes':144,
        'network_outcomes':1600,'tecator_conditions':6,'sachs_specifications':3,
        'new_simulations':False,'bootstrap_resampling':False}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(verify(),indent=2))
