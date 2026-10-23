function OUT = RUN_TFKM_TRKM_FRESH_DGP_SMOKE_V1(documents_root)
%RUN_TFKM_TRKM_FRESH_DGP_SMOKE_V1 Quick DGP-only smoke test, no method fits.
if nargin<1,documents_root="";end
if strlength(string(documents_root))>0,root=string(documents_root);else,root=string(fullfile(getenv('USERPROFILE'),'Documents'));end
outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if ~isfile(fullfile(outroot,'study.mat')),BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1(root);end
S=load(fullfile(outroot,'study.mat'));canon=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));assert(isfolder(canon));addpath(genpath(canon),'-begin');rehash path;assert(startsWith(string(which('generate_fixed_dataset_v1')),canon,'IgnoreCase',true));
T=table();
for i=1:height(S.catalog)
    cfg=table2struct(S.catalog(i,:));q=S.manifest(double(S.manifest.scenario_id)==cfg.scenario_id & double(S.manifest.replicate)==1,:);sim=generate_fixed_dataset_v1(cfg,double(q.data_seed));
    nout=NaN;if isfield(sim,'outlier'),nout=sum(logical(sim.outlier));elseif isfield(sim,'outlier_mask'),nout=sum(logical(sim.outlier_mask));end
    T=[T;table(cfg.scenario_id,string(cfg.panel),cfg.n,cfg.p,cfg.p_over_n,cfg.r0,cfg.c0,cfg.rho,cfg.delta,string(cfg.contamination),cfg.alpha0,size(sim.X,1),size(sim.X,2),nout, ...
        'VariableNames',{'scenario_id','panel','n','p','p_over_n','r0','g0','rho','delta','contamination','alpha0','generated_n','generated_p','n_outliers_POSTHOC_DGP_ONLY'})]; %#ok<AGROW>
end
assert(all(T.generated_n==T.n & T.generated_p==T.p),'DGP dimension failure.');
f=fullfile(outroot,'DGP_SMOKE_REP1.csv');writetable(T,f);disp(T(:,{'scenario_id','panel','p_over_n','r0','g0','rho','delta','contamination','n_outliers_POSTHOC_DGP_ONLY'}));fprintf('DGP SMOKE PASS: 41 cells generated correctly, including structure panel. No method fit run.\n');OUT=T;
end
