"""Small, independently testable helpers for saved-result replication."""
from __future__ import annotations
import hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

class ReplicationError(RuntimeError):
    """A missing input or failed reproducibility contract."""

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda:f.read(1024*1024),b''):
            h.update(part)
    return h.hexdigest()

def verify_integrity(root: Path) -> int:
    manifest=root/'documentation/SHA256SUMS.json'
    if not manifest.is_file():
        raise ReplicationError('Missing package integrity manifest.')
    records=json.loads(manifest.read_text(encoding='utf-8'))
    failures=[]
    for name,digest in records.items():
        p=root/name
        if not p.is_file():failures.append(name+' [missing]')
        elif sha256(p)!=digest:failures.append(name+' [changed]')
    if failures:
        raise ReplicationError('Package integrity failure:\n'+'\n'.join(failures[:20]))
    return len(records)

def load_csv(path: Path, required=(), n: int|None=None, keys=()) -> pd.DataFrame:
    if not path.is_file():raise ReplicationError(f'Missing required input: {path}')
    frame=pd.read_csv(path)
    missing=set(required)-set(frame.columns)
    if missing:raise ReplicationError(f'{path.name}: missing columns {sorted(missing)}')
    if n is not None and len(frame)!=n:
        raise ReplicationError(f'{path.name}: expected {n} rows, found {len(frame)}')
    if keys:
        if frame[list(keys)].isna().any().any():raise ReplicationError(f'{path.name}: null key')
        if frame.duplicated(list(keys)).any():raise ReplicationError(f'{path.name}: duplicate key')
    return frame

def ari_from_contingency(counts) -> float:
    a=np.asarray(counts,dtype=float)
    if a.ndim!=2 or not np.all(np.isfinite(a)) or np.any(a<0) or np.any(a!=np.floor(a)):
        raise ReplicationError('A contingency matrix must contain finite nonnegative integer counts.')
    n=a.sum()
    if n<2:raise ReplicationError('ARI requires at least two observations.')
    choose2=lambda x: x*(x-1)/2
    actual=choose2(a).sum(); left=choose2(a.sum(axis=1)).sum();right=choose2(a.sum(axis=0)).sum()
    expected=left*right/choose2(n);maximum=(left+right)/2;den=maximum-expected
    if abs(den)<=np.finfo(float).eps*max(1,abs(maximum)):
        return float(abs(actual-maximum)<=np.finfo(float).eps*max(1,abs(maximum)))
    return float((actual-expected)/den)

def auc_rank(y,score) -> float:
    y=np.asarray(y);score=np.asarray(score,dtype=float)
    if not np.all(np.isin(y,[0,1])) or not np.isfinite(score).all():
        raise ReplicationError('AUC needs binary labels and finite scores.')
    n1=int(y.sum());n0=len(y)-n1
    if min(n0,n1)==0:raise ReplicationError('AUC needs both outcome classes.')
    ranks=pd.Series(score).rank(method='average').to_numpy()
    return float((ranks[y==1].sum()-n1*(n1+1)/2)/(n1*n0))

class Checks:
    def __init__(self):self.rows=[]
    def add(self,name,ok,detail=''):
        self.rows.append({'check':name,'status':'PASS' if bool(ok) else 'FAIL','detail':str(detail)})
        if not bool(ok):raise ReplicationError(f'{name}: {detail}')
    def close(self,name,value,expected,atol=1e-10):
        self.add(name,np.isfinite(value) and abs(float(value)-expected)<=atol,
                 f'computed={value!r}; expected={expected!r}; atol={atol:g}')
    def published(self,name,value,printed,digits):
        tol=0.5*10**(-digits)
        self.add(name,np.isfinite(value) and abs(float(value)-printed)<=tol+1e-12,
                 f'computed={value:.15g}; manuscript={printed}; display tolerance={tol:g}')
    def frame(self):return pd.DataFrame(self.rows)

def equal_values(a: pd.DataFrame,b: pd.DataFrame,columns) -> bool:
    for c in columns:
        x=a[c].reset_index(drop=True);y=b[c].reset_index(drop=True)
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            if not np.allclose(x,y,rtol=0,atol=1e-12,equal_nan=True):return False
        else:
            if not x.fillna('<NA>').astype(str).equals(y.fillna('<NA>').astype(str)):return False
    return True

def latex_table(frame: pd.DataFrame, digits: int=6, formats: dict|None=None) -> str:
    """No jinja2 dependency; round for display only. Full precision stays in CSV."""
    def text(v,col=None):
        if pd.isna(v):return '--'
        if isinstance(v,(float,np.floating)):
            if np.isinf(v):return r'$\infty$' if v>0 else r'$-\infty$'
            places=(formats or {}).get(col,digits)
            return f'{v:.{places}f}'
        s=str(v)
        for a,b in [('\\',r'\textbackslash{}'),('&',r'\&'),('%',r'\%'),('_',r'\_'),('#',r'\#')]:s=s.replace(a,b)
        return s
    cols='l'*len(frame.columns)
    rows=[r'\begin{tabular}{'+cols+'}',r'\toprule',' & '.join(map(text,frame.columns))+r' \\',r'\midrule']
    rows.extend(' & '.join(text(v,col) for v,col in zip(row,frame.columns))+r' \\' for row in frame.itertuples(index=False,name=None))
    return '\n'.join(rows+[r'\bottomrule',r'\end{tabular}',''])
