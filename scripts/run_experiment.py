"""Prepare isolated study workspaces; computation requires --execute explicitly."""
from pathlib import Path
import argparse, csv, hashlib, json, os, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parents[1]
STUDIES=('bma_scaling','bma_components','bma_selectors','bayesian_network','tecator','sachs')

def read_csv(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))

def tree(source,dest):
    shutil.copytree(source,dest,dirs_exist_ok=True)

def stage_core(dest,variant):
    tree(ROOT/'src/bma/common',dest)
    tree(ROOT/'src/bma'/variant,dest)
    tree(ROOT/'experiments/bma_components/config',dest/'config')
    # This is a release integrity manifest, not the original preregistration.
    with (dest/'config/frozen_manifest.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['path','md5'])
        for p in sorted(dest.rglob('*')):
            if p.is_file() and p.name!='frozen_manifest.csv':
                writer.writerow([str(p.relative_to(dest)),hashlib.md5(p.read_bytes()).hexdigest()])
    for name in ['results','report','status','build']:(dest/name).mkdir(exist_ok=True)

def prepare(study,work):
    work=Path(work).resolve()
    if work==ROOT or ROOT in work.parents:
        raise ValueError('Use an output directory outside the release package.')
    work.mkdir(parents=True,exist_ok=False)
    tree(ROOT/'experiments'/study,work)
    if study in ('bma_scaling','bma_selectors','tecator'):
        stage_core(work/'core','staged' if study=='bma_selectors' else 'newton')
    if study=='bma_components':stage_core(work,'staged')
    if study=='bma_scaling':
        (work/'scaling_results').mkdir()
        files={str(p.relative_to(work)):hashlib.sha256(p.read_bytes()).hexdigest()
               for p in work.rglob('*') if p.is_file()}
        (work/'scaling_freeze.json').write_text(json.dumps({'scope':'Prepared release workspace; original IDs and seeds preserved','files':files},indent=2))
    if study=='tecator':
        (work/'data').mkdir();shutil.copy2(ROOT/'data/tecator/tecator.csv',work/'data/tecator.csv')
    if study=='bayesian_network':
        (work/'sim/src').mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/'src/bayesian_network/exact_bayesian_network_helpers.R',work/'sim/src')
        shutil.copy2(ROOT/'src/bma/staged/vendor/region_restricted_smc.R',work/'sim/src')
    if study=='sachs':
        tree(ROOT/'src/sachs',work)
        tree(ROOT/'data/sachs',work/'data')
        (work/'config').mkdir();shutil.copy2(work/'protocol.json',work/'config/protocol.json')
        (work/'results').mkdir()
    (work/'release_workspace.json').write_text(json.dumps({'study':study,'execution_requested':False,'source_manifest_sha256':hashlib.sha256((ROOT/'MANIFEST.json').read_bytes()).hexdigest() if (ROOT/'MANIFEST.json').exists() else None},indent=2))
    return work

def execute(study,work,rscript):
    env=os.environ.copy()
    if study=='bma_scaling':env.update(LANG='C',LC_ALL='C',TZ='UTC')
    def call(script,*args,timeout=None):
        subprocess.run([rscript,str(script),*map(str,args)],cwd=work,env=env,check=True,timeout=timeout)
    if study=='bma_components':
        for job in read_csv(work/'config/jobs.csv'):
            call('worker.R','--job='+job['job'],timeout=14400)
        subprocess.run([sys.executable,'collect.py'],cwd=work,check=True)
    elif study=='bma_selectors':
        for job in read_csv(work/'jobs.csv'):call('run_job.R',job['job'],timeout=1800)
        combine(work.glob('results/*/outcomes.csv'),work/'outcomes.csv')
    elif study=='tecator':call('run.R')
    elif study=='bayesian_network':call('sim/run_bayesian_network_study.R')
    elif study=='sachs':subprocess.run([sys.executable,'network.py'],cwd=work,check=True)
    elif study=='bma_scaling':
        jobs=read_csv(work/'scaling_config/jobs.csv')
        for index,dataset in enumerate(read_csv(work/'scaling_config/datasets.csv'),1):
            call('make_scaling_data.R',index,timeout=120)
            # Pass row positions, retaining original job/index columns as provenance.
            for position,job in enumerate(jobs,1):
                if job['dataset']!=dataset['dataset']:continue
                try:call('run_scaling_job.R',position,timeout=420)
                except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as error:
                    saved=work/'scaling_results'/job['dataset']/(job['job']+'.rds')
                    if saved.exists() or getattr(error,'returncode',None)==9:
                        raise RuntimeError('Worker reported a correctness error or failed after saving a result; inspect the preserved record.') from error
                    target=work/'scaling_results'/job['dataset']/(job['job']+'.failure.json')
                    target.write_text(json.dumps({'job':job['job'],'status':'hard-worker-timeout' if isinstance(error,subprocess.TimeoutExpired) else 'worker-error','error':str(error)}))
            call('audit_scaling_dataset.R',index,timeout=1200)
        combine(work.glob('scaling_results/*/outcomes.csv'),work/'scaling_results/outcomes.csv')

def combine(paths,target):
    rows=[row for path in sorted(paths) for row in read_csv(path)]
    if not rows:raise ValueError('No completed outcome tables')
    with target.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',choices=STUDIES)
    parser.add_argument('--output',type=Path,required=True,help='New directory outside this package')
    parser.add_argument('--execute',action='store_true',help='Run NEW computations after preparation (can be expensive)')
    parser.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'))
    args=parser.parse_args();work=prepare(args.study,args.output)
    if args.execute:
        marker=json.loads((work/'release_workspace.json').read_text());marker['execution_requested']=True
        (work/'release_workspace.json').write_text(json.dumps(marker,indent=2))
        execute(args.study,work,args.rscript)
    print(('Executed' if args.execute else 'Prepared without executing')+': '+str(work))

if __name__=='__main__':main()
