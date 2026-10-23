function OUT=smoke_fixed_candidate(root)
%SMOKE_FIXED_CANDIDATE Exactly one candidate (r,g,alpha)=(1,2,0) per method.
% Retains the ORIGINAL 100-start allocations. This is an environment test,
% not an additional published experiment or a reduced-budget replacement.
S=load(fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output','study.mat'));
sid=double(S.manifest.scenario_id(1));rep=double(S.manifest.replicate(1));
cfg=table2struct(S.catalog(double(S.catalog.scenario_id)==sid,:));
tf=fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor');
tr=fullfile(root,'TRKM_TARGETED_MAX_BACKWARD_V1_1','vendor','TFKM_TRKM_UNSEEN_FINAL_v1');
addpath(genpath(tf),'-begin');rehash path;
assert(startsWith(which('run_hd_multistart_top_hier_stable4_unseen_v1'),tf,'IgnoreCase',true));
sim=generate_fixed_dataset_v1(cfg,double(S.manifest.data_seed(1)));X=double(sim.X);options=S.options;
TF=run_hd_multistart_top_hier_stable4_unseen_v1(X,1,2,0,options,double(S.manifest.tfkm_pair_seed(1)));
check_fit(TF,size(X,1));
rmpath(genpath(tf));addpath(genpath(tr),'-begin');rehash path;
assert(startsWith(which('run_hd_multistart_top_v1'),tr,'IgnoreCase',true));
TR=run_hd_multistart_top_v1(X,1,2,0,"REG_TRKM_HD",options,double(S.manifest.trkm_pair_seed(1)));
check_fit(TR,size(X,1));
OUT=struct('status','PASS','purpose','ENVIRONMENT_SMOKE_NOT_PAPER_RESULT','scenario_id',sid,'replicate',rep,'r',1,'g',2,'alpha',0,'starts_per_method',100);
d=fullfile(root,'replication_smoke');if ~isfolder(d),mkdir(d);end
f=fullfile(d,['candidate_smoke_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')) '.mat']);
save(f,'TF','TR','OUT','-v7.3');
fprintf('Candidate smoke passed: two candidates, original100-start allocations. Saved outside production fits: %s\n',f);
end
function check_fit(F,n)
assert(isfield(F,'labels')&&numel(F.labels)==n,'Invalid labels in smoke fit.');
assert(isfield(F,'trim')&&numel(F.trim)==n,'Invalid trimming mask in smoke fit.');
assert(isfield(F,'U')&&all(isfinite(F.U(:))),'Invalid subspace in smoke fit.');
assert(isfield(F,'objective')&&isfinite(F.objective),'Invalid fit objective.');
end
