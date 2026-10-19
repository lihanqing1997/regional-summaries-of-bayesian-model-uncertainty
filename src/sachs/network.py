"""Structure-modular DAG sums via unique longest-path layers, all parent sets."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import itertools, json, math, time, heapq
import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.optimize import milp, Bounds, LinearConstraint
from numba import njit
from prepare import ROOT, dump

@njit(cache=True)
def layer_log_sum(weights, required, forbidden):
    n=weights.shape[0]; N=1<<n
    # A[i,B] sums legal parent sets inside B, including required neighbours in B.
    A=np.full((n,N),-np.inf)
    for i in range(n):
        for B in range(N):
            if B & (1<<i): continue
            must=required[i]&B
            allowed=B & ~forbidden[i]
            if must & ~allowed: continue
            free=allowed & ~must
            sub=free
            total=-np.inf
            while True:
                total=np.logaddexp(total,weights[i,sub|must])
                if sub==0: break
                sub=(sub-1)&free
            A[i,B]=total
    F=np.full((N,N),-np.inf)
    ternary=np.zeros(N,np.int64)
    for B in range(N):
        for i in range(n):
            if B&(1<<i):ternary[B]+=3**i
    # Many layer transitions reuse the same (node, earlier set, previous layer).
    # A compact ternary key caches the positive parent sum, including fallbacks.
    hit_cache=np.full((n,3**n),np.nan)
    for S in range(1,N):
        T=S
        while T:
            valid=True; base=0.
            for i in range(n):
                if T & (1<<i):
                    if required[i]&T:
                        valid=False; break
                    base+=weights[i,0]
            if valid:
                B=S^T
                if B==0:
                    F[S,T]=base
                else:
                    U=B; total=-np.inf
                    while U:
                        prev=F[B,U]
                        if np.isfinite(prev):
                            product=prev
                            for i in range(n):
                                if not (T & (1<<i)): continue
                                key=ternary[B]+ternary[U]
                                val=hit_cache[i,key]
                                if np.isnan(val):
                                    val=A[i,B]
                                    if not (required[i]&U):
                                        difference=A[i,B^U]-val
                                        # Avoid subtractive cancellation: positive direct sum.
                                        if difference>=-1e-10:
                                            must=required[i]&B
                                            free=B & ~forbidden[i] & ~must
                                            sub=free; val=-np.inf
                                            while True:
                                                P=sub|must
                                                if P&U: val=np.logaddexp(val,weights[i,P])
                                                if sub==0: break
                                                sub=(sub-1)&free
                                        else:
                                            val+=np.log(-np.expm1(difference))
                                    hit_cache[i,key]=val
                                product+=val
                                if not np.isfinite(product): break
                            total=np.logaddexp(total,product)
                        U=(U-1)&B
                    F[S,T]=total
            T=(T-1)&S
    total=-np.inf
    for T in range(1,N):total=np.logaddexp(total,F[N-1,T])
    return total

def layer_sum(weights,required,forbidden):
    return np.exp(layer_log_sum(np.log(weights),required,forbidden))

def local_scores(data,ess,absence):
    n=data.shape[1]; out=np.full((n,1<<n),-np.inf)
    ratio=((1-absence)/2)/absence
    for i in range(n):
        for mask in range(1<<n):
            if mask&(1<<i): continue
            pars=[j for j in range(n) if mask&(1<<j)]
            q=3**len(pars)
            code=data[:,pars]@(3**np.arange(len(pars),dtype=np.int64)) if pars else np.zeros(len(data),int)
            counts=np.bincount(code*3+data[:,i],minlength=q*3).reshape(q,3)
            # Zero-count terms cancel exactly; accumulate only occupied configurations.
            counts=counts[counts.sum(1)>0]
            value=np.sum(gammaln(ess/q)-gammaln(ess/q+counts.sum(1)))
            value+=np.sum(gammaln(ess/(3*q)+counts)-gammaln(ess/(3*q)))
            out[i,mask]=value+len(pars)*np.log(ratio)
    return out

@njit(cache=True)
def acyclic(parents,n):
    done=0
    for k in range(n):
        old=done
        for i in range(n):
            if not(done&(1<<i)) and (parents[i]&~done)==0: done|=1<<i
        if done==old: break
    return done==(1<<n)-1

@njit(cache=True)
def enumerate_skeleton(edges,scores):
    n=scores.shape[0]; E=len(edges); size=1<<E
    logs=np.full(size,-np.inf)
    for orient in range(size):
        pars=np.zeros(n,np.int64)
        for e in range(E):
            i,j=edges[e]
            if orient&(1<<e): pars[j]|=1<<i
            else: pars[i]|=1<<j
        if acyclic(pars,n):
            value=0.
            for i in range(n): value+=scores[i,pars[i]]
            logs[orient]=value
    return logs

def select_mass(mass,cost,target,mode):
    mass=np.asarray(mass);cost=np.asarray(cost)
    select_mass.last_status={'solver':'ranking','optimality_proven':False}
    if mass.sum()<target-1e-10: return None
    order=np.argsort(-(mass if mode=='ranked' else mass/cost),kind='stable')
    if mode in ('ranked','ratio'):
        return order[:np.searchsorted(np.cumsum(mass[order]),target)+1]
    res=milp(cost.astype(float),integrality=np.ones(len(cost)),bounds=Bounds(0,1),
       constraints=LinearConstraint(mass,target,np.inf),options={'time_limit':30,'mip_rel_gap':0})
    select_mass.last_status={'solver':'HiGHS MILP','status':int(res.status),'message':str(res.message),
        'optimality_proven':bool(res.status==0),'mip_gap':float(res.mip_gap) if getattr(res,'mip_gap',None) is not None else None}
    if res.x is None:return None
    chosen=np.flatnonzero(res.x>.5)
    if mass[chosen].sum()<target-1e-8: raise RuntimeError('MILP feasibility check failed')
    return chosen

def fractional_bound(mass,cost,target):
    order=np.argsort(-np.asarray(mass)/np.asarray(cost))
    remaining=target;value=0.
    for j in order:
        if mass[j]<=0:continue
        use=min(remaining,mass[j]); value+=use/mass[j]*cost[j]; remaining-=use
        if remaining<=1e-12: break
    return value if remaining<=1e-10 else None

class RegionalSum:
    def __init__(self,scores,folder):
        self.scores=scores; self.n=len(scores)
        self.pairs=list(itertools.combinations(range(self.n),2))
        self.shift=np.max(scores,axis=1)
        self.weights=scores-self.shift[:,None]
        self.folder=folder;self.cache={};self.calls=0
        self.path=folder/'mass_cache.json'
        if self.path.exists():self.cache=json.loads(self.path.read_text())
    def mass(self,present,absent):
        key=f'{present}:{absent}'
        if key not in self.cache:
            req=np.zeros(self.n,np.int64);forbid=np.zeros(self.n,np.int64)
            for e,(i,j) in enumerate(self.pairs):
                for bits,arr in ((present,req),(absent,forbid)):
                    if bits&(1<<e):arr[i]|=1<<j;arr[j]|=1<<i
            z=float(layer_log_sum(self.weights,req,forbid))
            if not np.isfinite(z): raise ArithmeticError(f'Invalid log DAG sum {z} at {key}')
            self.cache[key]=z; self.calls+=1
            if self.calls%20==0:dump(self.path,self.cache)
        return self.cache[key]

def run_variant(config,targets,maxcalls,maxseconds,orientation_cap):
    folder=ROOT/'results'/f"sachs_{config['name']}";folder.mkdir(exist_ok=True)
    original=json.loads((folder/'summary.json').read_text()) if (folder/'summary.json').exists() else None
    started=time.perf_counter()
    data=np.load(ROOT/'data/processed/sachs.npz')
    scores=local_scores(data['discrete'],config['ess'],config['absence'])
    np.save(folder/'local_scores.npy',scores)
    oracle=RegionalSum(scores,folder)
    z=oracle.mass(0,0); logz=z+oracle.shift.sum()
    prep=time.perf_counter()-started
    pairprobs=np.array([-np.expm1(oracle.mass(0,1<<e)-z) for e in range(len(oracle.pairs))])
    np.savetxt(folder/'edge_probabilities.csv',np.column_stack([np.array(oracle.pairs)+1,pairprobs]),delimiter=',',header='node1,node2,posterior_probability',comments='')
    queue=[(-1.,0,0,0)]; leaves=[]; checks=[]; discovered=0.;expanded=0
    while queue and discovered<max(targets)+1e-8 and len(oracle.cache)+2<=maxcalls and time.perf_counter()-started<maxseconds:
        negmass,depth,present,absent=heapq.heappop(queue);pm=-negmass
        if depth==len(oracle.pairs):
            leaves.append((present,pm));discovered+=pm;continue
        bit=1<<depth
        m0=np.exp(oracle.mass(present,absent|bit)-z)
        m1=np.exp(oracle.mass(present|bit,absent)-z)
        error=abs(m0+m1-pm)
        checks.append(error)
        if error>1e-8*max(pm,1e-8):raise ArithmeticError(f'Partition check failed: {error}')
        heapq.heappush(queue,(-m0,depth+1,present,absent|bit))
        heapq.heappush(queue,(-m1,depth+1,present|bit,absent))
        expanded+=1
        if expanded%25==0:
            dump(folder/'progress.json',dict(expanded=expanded,leaves=len(leaves),discovered_mass=discovered,elapsed=time.perf_counter()-started))
            print(config['name'],expanded,len(leaves),round(discovered,6),round(time.perf_counter()-started,1),flush=True)
    # Incorporate any complete skeletons already in the frontier.
    newqueue=[]
    for entry in queue:
        nm,d,p,a=entry
        if d==len(oracle.pairs):leaves.append((p,-nm));discovered-=nm
        else:newqueue.append(entry)
    queue=newqueue
    leaves.sort(key=lambda a:-a[1])
    mass=np.array([v for _,v in leaves]);cost=np.array([1+p.bit_count() for p,_ in leaves])
    table=[dict(skeleton_hex=hex(p),mass=float(m),cost=int(c)) for (p,m),c in zip(leaves,cost)]
    dump(folder/'skeletons.json',table)
    rows=[];selections={}
    fm=[-v[0] for v in queue];fc=[1+v[2].bit_count() for v in queue]
    for target in targets:
        lb=fractional_bound(list(mass)+fm,list(cost)+fc,target)
        # Every region with mass exceeding epsilon must be retained by any feasible union.
        mandatory=int(cost[mass>1-target+1e-10].sum())
        lb=max(lb,mandatory)
        for mode in ['ranked','ratio','optimized']:
            ix=select_mass(mass,cost,target,mode)
            row=dict(target=target,method=mode,status='inconclusive' if ix is None else 'mass_target_met_numerically',
                retained_mass=None,cost=None,cost_lower_bound=lb,mandatory_region_cost=mandatory,max_edge_error=None,
                selection_solver=select_mass.last_status.copy())
            if ix is not None:
                retained=float(mass[ix].sum());c=int(cost[ix].sum())
                marginal=sum(mass[k]*np.array([bool(leaves[k][0]&(1<<e)) for e in range(len(oracle.pairs))]) for k in ix)/retained
                row.update(retained_mass=retained,cost=c,max_edge_error=float(np.max(abs(marginal-pairprobs))),
                    cost_gap_upper=float(c-lb),global_cost_optimal_numerically=bool(c-np.ceil(lb-1e-9)<1e-9),entries=len(ix))
                selections[f'{target}_{mode}']=[hex(leaves[k][0]) for k in ix]
            rows.append(row)
    # Individual-DAG comparator, explicit enumeration only within discovered skeletons.
    dags=[];skipped=[];orientation_checks=[]
    for p,m in leaves:
        edges=np.array([pair for e,pair in enumerate(oracle.pairs) if p&(1<<e)],dtype=np.int64).reshape(-1,2)
        if 2**len(edges)>orientation_cap:skipped.append(hex(p));continue
        logs=enumerate_skeleton(edges,scores)
        ids=np.flatnonzero(np.isfinite(logs));probs=np.exp(logs[ids]-logz)
        orientation_checks.append(abs(probs.sum()-m))
        if abs(probs.sum()-m)>1e-8:raise ArithmeticError('Skeleton orientation sum mismatch')
        dags.extend((float(prob),len(edges)+1,hex(p),int(o)) for prob,o in zip(probs,ids))
    dags.sort(reverse=True)
    for target in targets:
        probs=np.array([a[0] for a in dags]);acc=np.cumsum(probs)
        k=int(np.searchsorted(acc,target))+1
        ok=k<=len(dags)
        # A ranked prefix is globally ranked only if every omitted DAG is bounded below cutoff.
        omitted_upper=max([*fm,*[m for p,m in leaves if hex(p) in skipped],0.])
        global_rank=bool(ok and dags[k-1][0]>=omitted_upper)
        rows.append(dict(target=target,method='ranked_dags',status=('mass_target_met_numerically' if ok else 'inconclusive'),
            retained_mass=float(acc[k-1]) if ok else None,cost=sum(d[1] for d in dags[:k]) if ok else None,
            entries=k if ok else None,global_ranking_verified=global_rank))
    dump(folder/'selections.json',selections)
    dump(folder/'comparison.json',rows)
    result=dict(config=config,normalizer_log=float(logz),normalizer_shifted_log=z,preparation_seconds=prep,
        total_seconds=time.perf_counter()-started,mass_calls=oracle.calls,cached_mass_regions=len(oracle.cache),
        expanded=expanded,discovered_skeletons=len(leaves),discovered_mass=float(discovered),
        remaining_frontier_mass=float(sum(fm)),partition_max_absolute_error=max(checks,default=0.),
        orientation_max_absolute_error=max(orientation_checks,default=0.),enumerated_dags=len(dags),
        skipped_skeletons=skipped,full_state_space='all DAGs on 11 labelled nodes',
        numerical_status='positive dynamic-program summation; floating-point validation, not interval arithmetic',
        stop_reason='target_mass_discovered' if discovered>=max(targets) else 'registered_resource_limit')
    if original:
        result['summary_replay_seconds']=result['total_seconds']
        result['total_seconds']=original['total_seconds']
        result['preparation_seconds']=original['preparation_seconds']
        result['mass_calls']=original['mass_calls']
        result['timing_note']='Original completed execution retained; replay uses mass cache. Primary execution reused 20 pre-repair cached regions; timings are not cold-start comparisons.'
    dump(folder/'summary.json',result);dump(oracle.path,oracle.cache)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    cfg=json.loads((ROOT/'config/protocol.json').read_text())
    for v in cfg['sachs']['variants']:
        if (ROOT/'results'/f"sachs_{v['name']}"/'summary.json').exists():continue
        run_variant(v,cfg['targets'],cfg['sachs']['max_mass_calls'],cfg['sachs']['max_seconds'],cfg['sachs']['orientation_cap'])
