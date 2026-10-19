"""Recreate processed Sachs data from the supplied raw copies."""
from pathlib import Path
import gzip,io,json
import numpy as np
ROOT=Path(__file__).resolve().parent
def dump(path,obj):Path(path).write_text(json.dumps(obj,indent=2),encoding='utf8')
def main():
    (ROOT/'data/processed').mkdir(parents=True,exist_ok=True)
    rawtext=gzip.decompress((ROOT/'data/raw/sachs.data.txt.gz').read_bytes()).decode()
    disctext=gzip.decompress((ROOT/'data/raw/sachs.discretised.txt.gz').read_bytes()).decode()
    header=[v.strip('"') for v in disctext.splitlines()[0].split()]
    raw=np.loadtxt(io.StringIO(rawtext),skiprows=1)
    labels=np.loadtxt(io.StringIO(disctext),skiprows=1,dtype=str)
    disc=np.empty(labels.shape,dtype=int)
    bin_labels={}
    for j,name in enumerate(header):
        levels=sorted(set(labels[:,j]),key=lambda s:float(s.strip('"')[1:].split(',')[0]))
        assert len(levels)==3
        bin_labels[name]=levels
        for k,level in enumerate(levels): disc[labels[:,j]==level,j]=k+1
    assert raw.shape==disc.shape==(853,11), (raw.shape,disc.shape)
    assert np.isfinite(raw).all() and np.all(np.isin(disc,[1,2,3]))
    assert [v.strip('"') for v in rawtext.splitlines()[0].split()]==header
    np.savez_compressed(ROOT/'data/processed/sachs.npz',raw=raw,discrete=disc-1,names=np.array(header))
    np.savetxt(ROOT/'data/processed/sachs_discrete.csv',disc,delimiter=',',fmt='%d',header=','.join(header),comments='')
    ranges={}
    for j,name in enumerate(header):
        ranges[name]={str(k):[float(raw[disc[:,j]==k,j].min()),float(raw[disc[:,j]==k,j].max()),int((disc[:,j]==k).sum())] for k in (1,2,3)}
    dump(ROOT/'data/processed/sachs_observed_bin_ranges.json',dict(labels=bin_labels,observed_ranges=ranges))

if __name__=='__main__':main()
