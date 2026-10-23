"""Read-only audit of saved fits and new banking outputs; no estimator fitting.
Raw matrices remain local. Output tables contain aggregates and variable codes only.
"""
from pathlib import Path
import argparse, hashlib, json, tempfile
import numpy as np
import pandas as pd
import scipy.io as sio
import h5py
from scipy.linalg import orth, svd

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,required=True,help='Directory containing uploaded banking ZIP/MAT files')
parser.add_argument('--sources',type=Path,required=True,help='Directory containing data_ECB_ANOR.mat and ECB_cov_input.mat or the R2 ECB ZIP')
parser.add_argument('--output',type=Path,required=True,help='New aggregate-output directory; must not already exist')
parser.add_argument('--previous-labels',type=Path,help='Optional historical anonymous label export; not needed for V15 replication')
args=parser.parse_args()
ROOT=args.inputs.resolve(); BASE=args.output.resolve();BASE.mkdir(parents=True,exist_ok=False)
_scratch=tempfile.TemporaryDirectory(prefix='csda_private_banking_')
EX=Path(_scratch.name)/'extracted';EX.mkdir()
import zipfile
for name in ['ANOR_POSTSELECTION_ANALYSIS_V1_output.zip','analysis_descriptive.zip','ANOR_TRKM_final_V1.zip','ECB_TFKM_final_V1(1).zip']:
 p=ROOT/name
 if not p.exists():raise FileNotFoundError(p)
 target=EX/p.stem;target.mkdir(exist_ok=True)
 with zipfile.ZipFile(p) as zz:
  if zz.testzip() is not None:raise ValueError('Corrupt ZIP: '+str(p))
  for zi in zz.infolist():
   dest=(target/zi.filename).resolve()
   if not dest.is_relative_to(target.resolve()):raise ValueError('Unsafe archive path')
   zz.extract(zi,target)
def locate_source(name):
 hits=list(args.sources.resolve().rglob(name))
 if len(hits)!=1:raise ValueError(f'Expected one {name}, found {len(hits)} in --sources')
 return hits[0]
try:LOG_SOURCE=locate_source('ECB_cov_input.mat')
except ValueError:
 zips=list(args.sources.resolve().glob('R2_ECB_APPLICATION_V12_CEIL_package*.zip'))
 if len(zips)!=1:raise FileNotFoundError('Provide ECB_cov_input.mat or one R2 ECB ZIP')
 with zipfile.ZipFile(zips[0]) as zz:
  hits=[n for n in zz.namelist() if Path(n).name=='ECB_cov_input.mat']
  if len(hits)!=1:raise ValueError('ECB_cov_input.mat missing or duplicated in ZIP')
  private=Path(_scratch.name)/'private_inputs';private.mkdir(exist_ok=True)
  LOG_SOURCE=private/'ECB_cov_input.mat';LOG_SOURCE.write_bytes(zz.read(hits[0]))
OUT=BASE/'tables';OUT.mkdir(exist_ok=True)
RATIO=next((EX/'ANOR_POSTSELECTION_ANALYSIS_V1_output').rglob('SUMMARY.csv')).parent
LOG=next((EX/'analysis_descriptive').rglob('ECB_TRKM_DESCRIPTIVE_SUMMARY.csv')).parent
checks=[];metrics=[]

def valmat(path,key):
 with h5py.File(path,'r') as f:
  d=f[key]
  if d.attrs.get('MATLAB_empty',0):return np.empty(0)
  a=np.array(d)
  return a.T if a.ndim==2 else a

def check(name,got,expected,source,tol=1e-10,relative=False):
 a=np.asarray(got);b=np.asarray(expected)
 if a.shape!=b.shape:passed=False;err=float('inf')
 else:
  err=float(np.nanmax(np.abs(a-b))) if a.size else 0.
  scale=max(1.,float(np.nanmax(np.abs(b)))) if b.size else 1.
  passed=bool(np.allclose(a,b,atol=tol,rtol=tol if relative else 0,equal_nan=True))
 checks.append(dict(check=name,pass_check=passed,max_absolute_error=err,tolerance=tol,relative=relative,source=source))
 if not passed:print('FAIL',name,err)

def metric(name,v,src):metrics.append(dict(metric=name,value=float(v),source=src))
def read(p):return pd.read_csv(p)
def orient(U):
 U=U.copy()
 for j in range(U.shape[1]):
  if U[np.argmax(np.abs(U[:,j])),j]<0:U[:,j]*=-1
 return U

def aricount(tab):
 tab=np.asarray(tab,dtype=float);n=tab.sum();c=lambda x:x*(x-1)/2
 a=c(tab).sum();b=c(tab.sum(1)).sum();d=c(tab.sum(0)).sum();e=b*d/c(n);m=(b+d)/2
 return (a-e)/(m-e) if m!=e else 1.
def corr(X,Z):
 X=np.asarray(X);Z=np.asarray(Z)
 if X.ndim==1:X=X[:,None]
 if Z.ndim==1:Z=Z[:,None]
 X=X-X.mean(0);Z=Z-Z.mean(0)
 den=np.linalg.norm(X,axis=0)[:,None]*np.linalg.norm(Z,axis=0)[None,:]
 return np.divide(X.T@Z,den,out=np.full((X.shape[1],Z.shape[1]),np.nan),where=den>0)
def rowsupport(X,U,ret,ratio=True):
 Y=X[ret];Y=Y-Y.mean(0);_,s,Vt=svd(Y,full_matrices=False)
 tol=max(Y.shape)*np.spacing(max(float(s.max()),1.) if ratio else float(s.max()))
 rk=int((s>tol).sum());Q=Vt[:rk].T; mass=np.clip(np.sum((Q.T@U)**2,axis=0),0,1)
 return rk,mass

# Ratios: exact source and prescribed row removal.
raw=sio.loadmat(locate_source('data_ECB_ANOR.mat'))['DataECB']
check('ratios raw shape',raw.shape,[365,383],'data_ECB_ANOR.mat',tol=0)
check('ratios duplicate 93=99',raw[92],raw[98],'data_ECB_ANOR.mat',tol=0)
check('ratios duplicate 108=126',raw[107],raw[125],'data_ECB_ANOR.mat',tol=0)
keep=np.ones(365,bool);keep[[92,107]]=False
X=raw[keep,1:];assets=raw[keep,0];ids=np.arange(1,366)[keep]
check('ratios unique clean rows',np.unique(X,axis=0).shape[0],363,'data_ECB_ANOR.mat',tol=0)
ref=ROOT/'ANOR_TFKM_DESCRIPTIVE_REFIT.mat'
yF=valmat(ref,'Labels').ravel().astype(int);tF=valmat(ref,'Trim').ravel().astype(bool);UF=orient(valmat(ref,'U'))
tr=EX/'ANOR_TRKM_final_V1/ANOR_TRKM_PRIMARY_FINAL.mat'
yR=valmat(tr,'FinalLabels').ravel().astype(int);tR=valmat(tr,'FinalTrim').ravel().astype(bool);UR=orient(valmat(tr,'FinalU'))
check('ratios TFKM no descriptive trimming',tF.sum(),0,ref.name,tol=0)
check('ratios TFKM Labels top vs fit',yF,valmat(ref,'fit/labels').ravel(),ref.name,tol=0)
check('ratios TFKM converged',valmat(ref,'fit/converged').item(),1,ref.name,tol=0)
check('ratios TFKM U orthonormal',UF.T@UF,np.eye(2),ref.name)
check('ratios TRKM U orthonormal',UR.T@UR,np.eye(1),tr.name)
part=read(EX/'ANOR_TRKM_final_V1/ANOR_TRKM_PRIMARY_PARTITION_ANONYMOUS.csv')
check('ratios row mapping',part.original_365_row.to_numpy(),ids,part.columns.tolist().__str__(),tol=0)
if args.previous_labels is not None:
 old=args.previous_labels
 o=read(old)
 labelcol='cluster' if 'cluster' in o else 'tfkm_cluster'
 check('ratios saved MAT labels vs optional historical labels',yF,o[labelcol].to_numpy(),old.name,tol=0)
ZF=X@UF;ZR=X@UR
T=pd.crosstab(yF[~tR],yR[~tR]);ari=aricount(T.to_numpy());over=float(np.sum((orth(UF).T@orth(UR))**2));angle=float(np.degrees(np.arccos(np.sqrt(over))))
check('ratios ARI vs archived',ari,read(RATIO/'SUBSPACE_AND_PARTITION_COMPARISON.csv').iloc[0,0],'SUBSPACE_AND_PARTITION_COMPARISON.csv')
check('ratios overlap vs archived',over,read(RATIO/'SUBSPACE_AND_PARTITION_COMPARISON.csv').iloc[0,1],'SUBSPACE_AND_PARTITION_COMPARISON.csv')
check('ratios angle vs archived',angle,read(RATIO/'SUBSPACE_AND_PARTITION_COMPARISON.csv').iloc[0,2],'SUBSPACE_AND_PARTITION_COMPARISON.csv')
for n,v in [('ratios_ARI',ari),('ratios_overlap_squared',over),('ratios_angle_degrees',angle)]:metric(n,v,'ANOR_TFKM_DESCRIPTIVE_REFIT.mat + ANOR_TRKM_PRIMARY_FINAL.mat')
T.to_csv(OUT/'RATIOS_CROSSTAB_RECALCULATED.csv',index_label='TFKM_cluster')
trimtab=pd.DataFrame([dict(TFKM_cluster=k,total=int((yF==k).sum()),TRKM_trimmed=int((tR&(yF==k)).sum()),TRKM_retained=int((~tR&(yF==k)).sum())) for k in np.unique(yF)])
trimtab.to_csv(OUT/'RATIOS_TRIMMING_RECALCULATED.csv',index=False)
check('ratios retained crosstab',read(RATIO/'TFKM_TRKM_CROSSTAB_ON_TRKM_RETAINED.csv').iloc[:,2].to_numpy(),T.to_numpy().ravel(),'TFKM_TRKM_CROSSTAB_ON_TRKM_RETAINED.csv',tol=0)
for m,U,Z,y,ret in [('TFKM',UF,ZF,yF,~tF),('TRKM',UR,ZR,yR,~tR)]:
 S=read(RATIO/f'{m}_SCORE_ASSET_ASSOCIATION.csv')
 c=corr(assets,Z).ravel();check(f'ratios {m} size correlations',c,S.corr_latent_score_log_total_assets.to_numpy(),f'{m}_SCORE_ASSET_ASSOCIATION.csv')
 for j,v in enumerate(c):metric(f'ratios_{m}_size_corr_{j+1}',v,f'{m}_SCORE_ASSET_ASSOCIATION.csv')
 rk,mass=rowsupport(X,U,ret)
 S=read(RATIO/f'{m}_ROWSPACE_SUPPORT.csv')
 check(f'ratios {m} row rank',rk,S.rank_centered_X.iloc[0],f'{m}_ROWSPACE_SUPPORT.csv',tol=0)
 check(f'ratios {m} row mass',mass,S.rowspace_mass.to_numpy(),f'{m}_ROWSPACE_SUPPORT.csv',tol=2e-10)
 for j,v in enumerate(1-mass):metric(f'ratios_{m}_null_mass_{j+1}',v,f'{m}_ROWSPACE_SUPPORT.csv')
 S=read(RATIO/f'{m}_TOP_LOADINGS.csv')
 got=np.array([U[int(r.variable_index)-1,int(r.component)-1] for r in S.itertuples()]);check(f'ratios {m} all listed loadings',got,S.loading.to_numpy(),f'{m}_TOP_LOADINGS.csv')
 C=corr(X,Z);S=read(RATIO/f'{m}_TOP_SCORE_CORRELATIONS.csv');got=np.array([C[int(r.variable_index)-1,int(r.component)-1] for r in S.itertuples()]);check(f'ratios {m} all listed variable-score correlations',got,S.corr_with_latent_score.to_numpy(),f'{m}_TOP_SCORE_CORRELATIONS.csv')
 S=read(RATIO/('TFKM_DESCRIPTIVE_CLUSTER_SUMMARY.csv' if m=='TFKM' else 'TRKM_PRIMARY_CLUSTER_SUMMARY.csv'))
 for r in S.itertuples():
  ix=(y==r.cluster)&ret
  check(f'ratios {m} cluster {r.cluster} summary',[ix.sum(),np.mean(assets[ix]),np.median(assets[ix])], [r.n,r.mean_log_total_assets,r.median_log_total_assets],m+'_CLUSTER_SUMMARY')
  check(f'ratios {m} cluster {r.cluster} centroid',np.mean(Z[ix],0),np.array([getattr(r,f'latent_centroid_{j+1}') for j in range(Z.shape[1])]),m+'_CLUSTER_SUMMARY')
S=read(RATIO/'TRKM_RETAINED_TRIMMED_ASSETS.csv')
for i,mask in enumerate([~tR,tR]):
 check('ratios assets '+S.status.iloc[i],[mask.sum(),assets[mask].mean(),np.median(assets[mask])],S.iloc[i,1:].to_numpy(float),'TRKM_RETAINED_TRIMMED_ASSETS.csv')
# All published contrast rows, including their full numeric decompositions.
for fname in ['TFKM_TOP_MULTIGROUP_CONTRASTS.csv','TRKM_TOP_CLUSTER_CONTRASTS.csv','TRKM_TOP_TRIMMED_VS_RETAINED_CONTRASTS.csv']:
 S=read(RATIO/fname);diffs=[]
 for row in S.itertuples():
  j=int(row.variable_index)-1
  if fname.startswith('TFKM'):
   a=X[yF==row.cluster_1,j];b=X[yF==row.cluster_2,j];truth=[row.standardized_difference_2_minus_1,row.mean_cluster_1,row.mean_cluster_2,row.nonzero_rate_cluster_1,row.nonzero_rate_cluster_2]
  else:
   if fname=='TRKM_TOP_CLUSTER_CONTRASTS.csv':a=X[(yR==row.cluster_1)&~tR,j];b=X[(yR==row.cluster_2)&~tR,j]
   else:a=X[~tR,j];b=X[tR,j]
   truth=[row.standardized_difference_B_minus_A,row.mean_A,row.mean_B,row.nonzero_rate_A,row.nonzero_rate_B]
  sd=np.sqrt((np.var(a)+np.var(b))/2);d=(b.mean()-a.mean())/sd if sd>0 else 0.
  diffs.append(np.array([d,a.mean(),b.mean(),np.mean(a!=0),np.mean(b!=0)])-truth)
 check('ratios all numeric rows '+fname,np.array(diffs),np.zeros_like(diffs),fname)

# Log-level input from previously supplied ECB covariance package.
rawfile=LOG_SOURCE
D=sio.loadmat(rawfile);rawL=D['X'];keepL=np.ones(364,bool);keepL[92]=False;XL=rawL[keepL];idsL=np.arange(1,365)[keepL]
check('log raw covariance relative error',np.linalg.norm(np.cov(rawL,rowvar=False)-D['C'])/np.linalg.norm(D['C']),0.,'ECB_cov_input.mat',tol=1e-10)
check('log duplicate 93=99',rawL[92],rawL[98],'ECB_cov_input.mat',tol=0)
check('log clean unique rows',len(np.unique(XL,axis=0)),363,'ECB_cov_input.mat',tol=0)
metric('log_zero_fraction',np.mean(XL==0),'ECB_cov_input.mat / prescribed row removal')
U=valmat(ROOT/'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat','Uor');mu=valmat(ROOT/'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat','mu').ravel();cent=valmat(ROOT/'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat','Centroids').ravel()
check('log saved center vs raw clean mean',XL.mean(0),mu,'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat + ECB_cov_input.mat')
Z=(XL-mu)@U;dist=(Z-cent[None,:])**2
# Reconstruct assignments using archived centers; this is not a refit.
y=np.argmin(dist,axis=1)+1;resid=dist.min(1)
check('log saved residual vector',resid,valmat(ROOT/'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat','R2all').ravel(),'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat',tol=1e-9)
check('log descriptive counts from nearest saved centers',np.bincount(y)[1:],[70,293],'ECB_TRKM_DESCRIPTIVE_CLUSTER_SIZES.csv',tol=0)
for k in [1,2]:check(f'log center {k} reconstructed',float(Z[y==k].mean()),cent[k-1],'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat')
F=EX/'ECB_TFKM_final_V1(1)/ECB_TFKM_PRIMARY_FINAL.mat';ylF=valmat(F,'FinalLabels').ravel().astype(int);tlF=valmat(F,'FinalTrim').ravel().astype(bool);UlF=valmat(F,'FinalU')
part=read(EX/'ECB_TFKM_final_V1(1)/ECB_TFKM_PRIMARY_PARTITION_ANONYMOUS.csv')
mapcol='original_364_row' if 'original_364_row' in part.columns else 'original_row_index'
check('log row mapping',idsL,part[mapcol].to_numpy(),part.columns.tolist().__str__(),tol=0)
TL=pd.crosstab(ylF[~tlF],y[~tlF]);ariL=aricount(TL.to_numpy());ovL=float(np.sum((orth(UlF).T@orth(U))**2));angL=float(np.degrees(np.arccos(np.sqrt(ovL))))
S=read(LOG/'ECB_TFKM_TRKM_SECONDARY_COMPARISON.csv')
check('log crossmethod full row',[int((~tlF).sum()),int(tlF.sum()),ariL,ovL,angL],S.iloc[0].to_numpy(float),'ECB_TFKM_TRKM_SECONDARY_COMPARISON.csv')
for n,v in [('log_ARI',ariL),('log_overlap_squared',ovL),('log_angle_degrees',angL)]:metric(n,v,'ECB_TFKM_PRIMARY_FINAL.mat + saved TRKM U/centers + raw X')
TL.to_csv(OUT/'LOGLEVEL_CROSSTAB_RECALCULATED.csv',index_label='TFKM_cluster')
trimtabL=pd.DataFrame([dict(TRKM_cluster=k,total=int((y==k).sum()),TFKM_trimmed=int((tlF&(y==k)).sum()),TFKM_retained=int((~tlF&(y==k)).sum())) for k in [1,2]])
trimtabL.to_csv(OUT/'LOGLEVEL_TRIMMING_RECALCULATED.csv',index=False)
S=read(LOG/'ECB_TFKM_TRIM_RATE_BY_TRKM_DESCRIPTIVE_CLUSTER.csv')
check('log trimming counts',trimtabL[['total','TFKM_trimmed']].to_numpy(),S[['n_total','n_TFKM_trimmed']].to_numpy(),'ECB_TFKM_TRIM_RATE_BY_TRKM_DESCRIPTIVE_CLUSTER.csv',tol=0)
sizecorr=float(corr(XL[:,0],Z).item());S=read(LOG/'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv')
check('log TRKM size corr',sizecorr,S.corr_latent_score_log_total_assets.iloc[0],'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv')
metric('log_TRKM_size_corr',sizecorr,'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv')
rk,mass=rowsupport(XL,U,np.ones(363,bool),False)
check('log TRKM row rank',rk,S.rank_centered_X.iloc[0],'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv',tol=0)
check('log TRKM row mass',mass[0],S.mean_rowspace_mass.iloc[0],'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv')
metric('log_TRKM_null_mass',1-mass[0],'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv')
# TFKM size result not included in new CSVs: compare to manuscript rounded values.
cf=float(corr(XL[~tlF,0],XL[~tlF]@UlF).item());metric('log_TFKM_retained_size_abs_corr',abs(cf),'raw X + ECB_TFKM_PRIMARY_FINAL.mat')
check('log TFKM retained abs size corr vs manuscript',round(abs(cf),3),.475,'V15 S9.2',tol=0)
for k,expected in [(1,21.64),(2,23.29)]:
 v=XL[(ylF==k)&~tlF,0].mean();metric(f'log_TFKM_cluster_{k}_mean_assets',v,'raw X + ECB_TFKM_PRIMARY_FINAL.mat');check(f'log TFKM group {k} mean assets vs manuscript',round(v,2),expected,'V15 S9.2',tol=0)
S=read(LOG/'ECB_TRKM_DESCRIPTIVE_SIZE_BY_CLUSTER.csv')
for row in S.itertuples():
 ix=y==row.cluster;check('log TRKM group '+str(row.cluster)+' assets',[ix.sum(),XL[ix,0].mean(),np.median(XL[ix,0]),Z[ix].mean()],[row.n,row.mean_log_total_assets,row.median_log_total_assets,row.mean_latent_score_1],S.columns.tolist().__str__())
S=read(LOG/'ECB_TRKM_DESCRIPTIVE_SCORE_CORRELATIONS.csv');CC=corr(XL,Z).ravel()
check('log all listed score correlations',CC[S.variable_index.to_numpy(int)-1],S.corr_with_latent_score.to_numpy(),'ECB_TRKM_DESCRIPTIVE_SCORE_CORRELATIONS.csv')
S=read(LOG/'ECB_TRKM_DESCRIPTIVE_RAW_LOADINGS.csv');check('log all listed raw loadings',U[S.variable_index.to_numpy(int)-1,0],S.loading.to_numpy(),'ECB_TRKM_DESCRIPTIVE_RAW_LOADINGS.csv')
S=read(LOG/'ECB_TRKM_DESCRIPTIVE_DATA_SUPPORTED_LOADINGS.csv');Ur=valmat(ROOT/'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat','Urow');check('log listed data-supported loadings',Ur[S.variable_index.to_numpy(int)-1,0],S.loading.to_numpy(),'ECB_TRKM_DESCRIPTIVE_DATA_SUPPORTED_LOADINGS.csv')
S=read(LOG/'ECB_TRKM_DESCRIPTIVE_CONTRAST_ZERO_DECOMPOSITION.csv');a=XL[y==1];b=XL[y==2];sd=XL.std(0,ddof=1);diff=[]
for r in S.itertuples():
 j=int(r.variable_index)-1;av=a[:,j];bv=b[:,j];got=[(bv.mean()-av.mean())/sd[j],np.mean(av!=0),np.mean(bv!=0),np.median(av[av!=0]),np.median(bv[bv!=0])];exp=[r.standardized_difference_2_minus_1,r.nonzero_rate_cluster_1,r.nonzero_rate_cluster_2,r.median_nonzero_cluster_1,r.median_nonzero_cluster_2];diff.append(np.array(got)-exp)
check('log all listed zero/contrast decompositions',np.array(diff),np.zeros_like(diff),'ECB_TRKM_DESCRIPTIVE_CONTRAST_ZERO_DECOMPOSITION.csv')
# Save aggregates only.
pd.DataFrame(checks).to_csv(OUT/'BANKING_PRIVATE_CHECKS.csv',index=False)
pd.DataFrame(metrics).to_csv(OUT/'BANKING_PRIVATE_METRICS.csv',index=False)
print(pd.DataFrame(metrics).to_string(index=False))
print('CHECKS',len(checks),'PASSED',sum(c['pass_check'] for c in checks))
assert all(c['pass_check'] for c in checks), 'Inspect banking check failures before reporting.'

_scratch.cleanup()
