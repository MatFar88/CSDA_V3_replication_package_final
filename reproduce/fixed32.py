"""FINAL32 fixed-parameter comparison: 4 scenarios x 8 replicates x 6 methods.

The supplied merged records are the sole reporting input. Worker files and
published summaries are independent cross-checks, never additional observations.
The method identities correspond to the final HIER4/RV-CV procedures.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from .common import load_csv, Checks, equal_values, ReplicationError

KEY=['scenario_id','replicate','method']
METHODS={
 'TKM':'TKM', 'PCA_TKM':'PCA-TKM', 'RPCA_TKM':'ROBPCA-TKM',
 'TCLUST':'TCLUST', 'TFKM_V2_HIER4':'REG-HD-TFKM', 'TRKM_RVCV':'REG-HD-TRKM'}
SCENARIOS={33001,33004,33005,33007}

def aggregate(frame, groups):
    """Rebuild the supplied summaries with their original names/denominators."""
    spec={
      'n':('ari_inliers_full','size'),
      'mean_ari_inliers_full':('ari_inliers_full','mean'),
      'min_ari_inliers_full':('ari_inliers_full','min'),
      'mean_ari_inliers_retained':('ari_inliers_retained','mean'),
      'mean_trim_f1':('trim_f1','mean'),
      'outlier_set_exact_rate':('outlier_set_exact','mean'),
      'mean_false_trim_rate':('false_trim_rate','mean'),
      'mean_masking_rate':('masking_rate','mean'),
      'mean_subspace_projection_distance':('subspace_projection_distance','mean'),
      'mean_subspace_overlap':('subspace_overlap','mean'),
      'mean_principal_angle_degrees':('mean_principal_angle_degrees','mean'),
      'mean_maximum_principal_angle_degrees':('maximum_principal_angle_degrees','mean'),
      'mean_centroid_nrmse_observed':('centroid_nrmse_observed','mean'),
      'mean_centroid_nrmse_latent_aligned':('centroid_nrmse_latent_aligned','mean'),
      'mean_runtime_seconds':('runtime_seconds','mean')}
    return frame.groupby(groups,sort=True).agg(**spec).reset_index()

def fixed_tables(root:Path,path:Path,checks:Checks,expected:dict):
    required=KEY+['session','valid','runtime_seconds','ari_inliers_full',
      'ari_inliers_retained','N_starts','internal_selector','r0','g0',
      'alpha0','contamination','truth_regime']
    t=load_csv(path,required,192,KEY)
    checks.add('FINAL32: four declared scenarios',set(t.scenario_id)==SCENARIOS)
    checks.add('FINAL32: six original method identities',set(t.method)==set(METHODS))
    checks.add('FINAL32: 32 unique data sets',len(t[['scenario_id','replicate']].drop_duplicates())==32)
    checks.add('FINAL32: eight replicates in every scenario/method',
      t.groupby(['scenario_id','method']).replicate.apply(lambda x:set(x)==set(range(1,9))).all())
    checks.add('FINAL32: six methods for every data set',t.groupby(KEY[:2]).method.nunique().eq(6).all())
    checks.add('FINAL32: all fits valid',t.valid.eq(1).all())
    checks.add('FINAL32: 100 starts recorded',t.N_starts.eq(100).all())
    checks.add('FINAL32: session letters match replicates',
      t.session.eq(t.replicate.map(dict(enumerate('ABCDEFGH',1)))).all())
    for col in ['ari_inliers_full','ari_inliers_retained']:
        checks.add('FINAL32: bounded finite '+col,np.isfinite(t[col]).all() and t[col].between(-1,1).all())
    checks.add('FINAL32: finite nonnegative runtime',np.isfinite(t.runtime_seconds).all() and t.runtime_seconds.ge(0).all())
    checks.add('FINAL32: HIER4 arbitration',t.loc[t.method.eq('TFKM_V2_HIER4'),'internal_selector'].eq('HIER_STABLE4').all())
    checks.add('FINAL32: RV_CV arbitration',t.loc[t.method.eq('TRKM_RVCV'),'internal_selector'].eq('RV_CV').all())
    for col in ['r0','g0','alpha0','contamination','truth_regime']:
        checks.add('FINAL32: common generating '+col,t.groupby('scenario_id')[col].nunique().eq(1).all())
    checks.add('FINAL32: recorded nominal contamination levels',np.array_equal(t.alpha0,t.scenario_id.map({33001:0.,33004:.1,33005:.1,33007:.05})))
    workers=[]
    for letter in 'ABCDEFGH':
        f=root/'results/reference/fixed32/workers'/f'fixed_{letter}.csv'
        w=load_csv(f,t.columns,24,KEY);workers.append(w)
        checks.add('FINAL32 worker '+letter+' matches merged records',
          equal_values(w.sort_values(KEY),t[t.session.eq(letter)].sort_values(KEY),list(t.columns)))
    all_workers=pd.concat(workers,ignore_index=True)
    checks.add('FINAL32: workers are exact coverage, not extra observations',
      not all_workers.duplicated(KEY).any() and equal_values(all_workers.sort_values(KEY),t.sort_values(KEY),list(t.columns)))
    tables={}
    for suffix,groups,n in [('summary',['method'],6),('by_scenario',['method','scenario_id'],24),('by_regime',['method','truth_regime'],18)]:
        result=aggregate(t,groups)
        original=load_csv(root/'results/reference/fixed32'/f'FIXED_{suffix}.csv',groups,n,groups)
        checks.add('FINAL32 '+suffix+': original summary columns',set(original.columns)==set(result.columns))
        a=result.sort_values(groups);b=original.sort_values(groups)
        for col in original.columns:
            checks.add('FINAL32 '+suffix+': recomputed '+col,equal_values(a,b,[col]))
        tables['FIXED32_'+suffix]=result
    out=t.groupby('method').agg(N=('ari_inliers_full','size'),Mean_ARI=('ari_inliers_full','mean'),
      Median_ARI=('ari_inliers_full','median'),Mean_seconds=('runtime_seconds','mean')).reindex(METHODS).reset_index()
    out['method']=out.method.map(METHODS)
    for row in out.itertuples(index=False):
        for metric,value,target,digits in zip(['mean','median','runtime'],[row.Mean_ARI,row.Median_ARI,row.Mean_seconds],expected['table_s4'][row.method],[4,4,3]):
            checks.published('TableS4 FINAL32 '+row.method+' '+metric,value,target,digits)
    tables={'Supplement_TableS4':out,**tables}
    return tables
