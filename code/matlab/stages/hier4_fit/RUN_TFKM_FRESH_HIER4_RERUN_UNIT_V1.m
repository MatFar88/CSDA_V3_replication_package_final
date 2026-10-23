function OUT = RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1(scenario_id,replicate,documents_root)
%RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1
% Corrected pre-truth TFKM fresh rerun using the already-frozen HIER_STABLE4
% multistart arbitration. Same DGP, data seeds, pair seeds, alpha/RG grid,
% initialization families, allocations, ridge and ALS core as the original
% fresh study. Writes to a NEW output tree; never overwrites the aborted
% SIGN_ALL_STABLE TFKM cubes.
%
% NO truth labels, r0/g0/alpha0, contamination labels, ARI, recovery or
% outlier masks are read, computed or printed here.
%
% Usage:
%   RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1(51001,1);
%   RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1("SELFTEST",[]);
%   RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1("SOURCECHECK",[]);

if nargin>=1 && strcmpi(string(scenario_id),"SELFTEST")
    RG=rg_grid_local();
    assert(size(RG,1)==15 && isequal(RG(1,:),[1 2]) && isequal(RG(end,:),[5 6]));
    assert(abs(hier_identity_local(.4,.7,.4+.2*.7))<=1e-12);
    fprintf('SELFTEST PASS: HIER4 rerun grid + score identity contract.\n');
    OUT=struct('selftest_pass',true);return
end
if nargin<3,documents_root="";end
if nargin>=1 && strcmpi(string(scenario_id),"SOURCECHECK")
    [root,~,~]=load_study_local(documents_root);
    P=lock_tfkm_hier4_vendor_local(root);
    fprintf('\n============================================================\n');
    fprintf('TFKM HIER4 RERUN SOURCE CHECK - TRUTH SEALED\n');
    fprintf('============================================================\n');
    disp(P);
    fprintf('STATUS: HIER4_SOURCE_LOCK_CONFIRMED\n');
    OUT=struct('source_check_pass',true,'paths',P);return
end
if nargin<2,error('Provide scenario_id and replicate.');end

[root,srcroot,study]=load_study_local(documents_root);
options=study.options;catalog=study.catalog;manifest=study.manifest;CFG=study.CFG;
validate_options_local(options);
assert(replicate>=1 && replicate<=CFG.n_replicates && replicate==round(replicate),'replicate out of range.');
ic=find(double(catalog.scenario_id)==double(scenario_id),1);assert(~isempty(ic),'Unknown scenario_id.');
cfg=table2struct(catalog(ic,:));
im=find(double(manifest.scenario_id)==double(scenario_id) & double(manifest.replicate)==double(replicate),1);assert(~isempty(im),'Manifest row missing.');
data_seed=double(manifest.data_seed(im));
pair_seed=double(manifest.tfkm_pair_seed(im));

paths=lock_tfkm_hier4_vendor_local(root);
sim=generate_fixed_dataset_v1(cfg,data_seed);X=sim.X;n=size(X,1);p=size(X,2);
assert(n==cfg.n && p==cfg.p,'Generated dimension drift.');
alpha_grid=double(CFG.alpha_grid);RG=rg_grid_local();N=150;

outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
unitdir=fullfile(outroot,'fits','TFKM',sprintf('s%d',scenario_id),sprintf('rep%02d',replicate));
if ~isfolder(unitdir),mkdir(unitdir);end
matfile=fullfile(unitdir,'cube.mat');csvfile=fullfile(unitdir,'candidates.csv');completefile=fullfile(unitdir,'COMPLETE.txt');
if isfile(completefile) && isfile(matfile)
    Q=load(matfile,'Candidates','unit_info');
    fprintf('TFKM-HIER4 s%d rep%02d already COMPLETE (%d/150).\n',scenario_id,replicate,height(Q.Candidates));
    OUT=struct('complete',true,'unit_dir',string(unitdir),'candidates',Q.Candidates,'unit_info',Q.unit_info);return
end
if isfile(matfile)
    Q=load(matfile,'Candidates','Labels','Trim','unit_info');Candidates=Q.Candidates;Labels=Q.Labels;Trim=Q.Trim;unit_info=Q.unit_info;
else
    Candidates=empty_candidate_table_local();Labels=zeros(n,N,'uint16');Trim=false(n,N);
    hier_row = paths.function_name == "run_hd_multistart_top_hier_stable4_unseen_v1";
    assert(sum(hier_row)==1,'HIER4 wrapper provenance row is not unique.');
    hier_wrapper_path = string(paths.resolved_path(hier_row));
    unit_info=struct('study',CFG.study_name,'rerun_name',"TFKM_HIER4_RERUN_V1_PATCH1",'scenario_id',double(scenario_id), ...
        'replicate',double(replicate),'method',"REG_TFKM_HD",'data_seed',data_seed,'pair_seed',pair_seed,'n',n,'p',p, ...
        'source_study_root',string(srcroot),'hier_wrapper_path',hier_wrapper_path,'complete',false,'started_at',datetime('now'));
end

fprintf('\n============================================================\n');
fprintf('TFKM HIER4 FRESH RERUN | sid=%d rep=%02d | n=%d p=%d | existing=%d/150\n',scenario_id,replicate,n,p,height(Candidates));
fprintf('IMPLEMENTATION-CORRECTION RERUN ONLY | TRUTH PERFORMANCE SEALED.\n');
fprintf('Selector: HIER_STABLE4 | same fixed data/pair seeds and 50/35/10/5 starts.\n');
fprintf('============================================================\n');
row_index=0;
for ia=1:numel(alpha_grid)
    a=alpha_grid(ia);
    for irg=1:size(RG,1)
        row_index=row_index+1;r=RG(irg,1);g=RG(irg,2);
        if ~isempty(Candidates) && any(double(Candidates.row_index)==row_index),continue,end
        t0=tic;
        fit=run_hd_multistart_top_hier_stable4_unseen_v1(X,r,g,a,options,pair_seed);
        runtime=toc(t0);
        assert(isfield(fit,'valid') && logical(fit.valid),'Invalid HIER4 fit sid=%d rep=%d (%d,%d,%.3f).',scenario_id,replicate,r,g,a);
        stats=compute_external_rg_fit_stats_v1(X,fit,"REG_TFKM_HD");
        A=scalar_field_local(fit,'hier_cross_ari_score',NaN);
        J=scalar_field_local(fit,'hier_cross_j_score',NaN);
        H=scalar_field_local(fit,'hier_cross_H_score',NaN);
        S=scalar_field_local(fit,'selector_score',NaN);
        assert(all(isfinite([A J H S])),'Missing finite HIER4 evidence sid=%d rep=%d (%d,%d,%.3f).',scenario_id,replicate,r,g,a);
        assert(abs(H-(A+double(options.trim_stability_weight)*J))<=1e-10,'HIER4 A+0.2J identity failure.');
        assert(abs(S-H)<=1e-10,'selector_score != hier_cross_H_score.');
        Labels(:,row_index)=uint16(fit.labels(:));Trim(:,row_index)=logical(fit.trim(:));
        row=candidate_row_local(row_index,scenario_id,replicate,data_seed,pair_seed,n,p,r,g,a,stats,fit,runtime);
        Candidates=[Candidates;row]; %#ok<AGROW>
        Candidates=sortrows(Candidates,'row_index');
        unit_info.completed_candidates=height(Candidates);unit_info.last_saved_at=datetime('now');
        save(matfile,'Candidates','Labels','Trim','unit_info','-v7.3');writetable(Candidates,csvfile);
        fprintf('%3d/150 | a=%5.3f | (%d,%d) | H=% .6f | A=% .6f | J=% .6f | %s | %.1fs\n', ...
            row_index,a,r,g,S,A,J,string(row.selected_init_family),runtime);
    end
end
assert(height(Candidates)==150 && all(logical(Candidates.valid)),'Unit cube incomplete/invalid.');
assert(all(isfinite(double(Candidates.hier_cross_ari_score))) && all(isfinite(double(Candidates.hier_cross_j_score))),'Unit HIER evidence incomplete.');
err=abs(double(Candidates.selector_score)-(double(Candidates.hier_cross_ari_score)+double(options.trim_stability_weight)*double(Candidates.hier_cross_j_score)));
assert(max(err)<=1e-10,'Unit HIER identity audit failed.');
unit_info.complete=true;unit_info.completed_at=datetime('now');unit_info.completed_candidates=150;unit_info.max_hier_identity_error=max(err);
save(matfile,'Candidates','Labels','Trim','unit_info','-v7.3');writetable(Candidates,csvfile);
fid=fopen(completefile,'w');assert(fid>0);fprintf(fid,'COMPLETE\nrerun=TFKM_HIER4_RERUN_V1\nscenario_id=%d\nreplicate=%d\nmethod=REG_TFKM_HD\ntruth_opened=NO\n',scenario_id,replicate);fclose(fid);
fprintf('TFKM-HIER4 s%d rep%02d COMPLETE.\n',scenario_id,replicate);
OUT=struct('complete',true,'unit_dir',string(unitdir),'candidates',Candidates,'unit_info',unit_info);
end

function row=candidate_row_local(idx,sid,rep,ds,ps,n,p,r,g,a,stats,fit,runtime)
row=table(idx,double(sid),double(rep),double(ds),double(ps),double(n),double(p),"REG_TFKM_HD",double(r),double(g),double(a), ...
    scalar_field_local(stats,'R2',NaN),scalar_field_local(stats,'W_ratio',NaN),scalar_field_local(stats,'n_keep',NaN), ...
    scalar_field_local(fit,'selector_score',NaN),scalar_field_local(fit,'selector_cv_score',NaN),scalar_field_local(fit,'objective',NaN), ...
    string_field_local(fit,'selected_init_family'),scalar_field_local(fit,'selected_start_id',NaN),scalar_field_local(fit,'selected_objective_rank',NaN), ...
    scalar_field_local(fit,'valid_starts',NaN),scalar_field_local(fit,'N_starts',NaN),string_field_local(fit,'internal_selector'), ...
    scalar_field_local(fit,'hier_within_family_score',NaN),scalar_field_local(fit,'hier_cross_ari_score',NaN), ...
    scalar_field_local(fit,'hier_cross_j_score',NaN),scalar_field_local(fit,'hier_cross_H_score',NaN),true,double(runtime), ...
    'VariableNames',{'row_index','scenario_id','replicate','data_seed','pair_seed','n','p','method','r','g','alpha','R2','W_ratio','n_keep', ...
    'selector_score','selector_cv_score','objective','selected_init_family','selected_start_id','selected_objective_rank','valid_starts','N_starts','internal_selector', ...
    'hier_within_family_score','hier_cross_ari_score','hier_cross_j_score','hier_cross_H_score','valid','runtime_seconds'});
end
function T=empty_candidate_table_local()
T=table('Size',[0 29],'VariableTypes',{'double','double','double','double','double','double','double','string','double','double','double','double','double','double','double','double','double','string','double','double','double','double','string','double','double','double','double','logical','double'}, ...
    'VariableNames',{'row_index','scenario_id','replicate','data_seed','pair_seed','n','p','method','r','g','alpha','R2','W_ratio','n_keep','selector_score','selector_cv_score','objective','selected_init_family','selected_start_id','selected_objective_rank','valid_starts','N_starts','internal_selector','hier_within_family_score','hier_cross_ari_score','hier_cross_j_score','hier_cross_H_score','valid','runtime_seconds'});
end
function x=scalar_field_local(s,name,default)
x=default;if isstruct(s)&&isfield(s,name),v=s.(name);if ~isempty(v)&&isscalar(v)&&isnumeric(v),x=double(v);end,end
end
function x=string_field_local(s,name)
x="";if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),x=string(s.(name));if numel(x)>1,x=x(1);end,end
end
function e=hier_identity_local(A,J,H),e=H-(A+.2*J);end
function RG=rg_grid_local()
RG=zeros(0,2);for r=1:5,for g=2:6,if r<g,RG(end+1,:)=[r g];end,end,end
end
function [root,srcroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";srcroot="";
for c=cands(:)'
    o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
    if isfile(fullfile(o,'study.mat')),root=c;srcroot=string(o);break;end
end
assert(strlength(root)>0,'Original fresh study not found. The rerun reuses its immutable study.mat and manifest.');
S=load(fullfile(srcroot,'study.mat'));
end
function validate_options_local(options)
req={'hd_init_names','hd_tfkm_allocation','trim_stability_weight','tauW','tauT','max_iterations'};
for j=1:numel(req),assert(isfield(options,req{j}),'options missing %s.',req{j});end
assert(isequal(string(options.hd_init_names(:)'),["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"]),'Initialization family drift.');
assert(isequal(double(options.hd_tfkm_allocation(:)'),[50 35 10 5]),'TFKM allocation drift.');
assert(abs(double(options.trim_stability_weight)-.2)<=1e-12,'trim_stability_weight drift.');
end
function P=lock_tfkm_hier4_vendor_local(root)
tf=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));
tr=string(fullfile(root,'TRKM_TARGETED_MAX_BACKWARD_V1_1','vendor','TFKM_TRKM_UNSEEN_FINAL_v1'));
assert(isfolder(tf),'TFKM canonical vendor missing: %s',tf);
if isfolder(tr),safe_rmpath_tree_local(tr);end
addpath(genpath(char(tf)),'-begin');rehash path;
req=["generate_fixed_dataset_v1","regularized_trimmed_fit_corrected_v1","make_robust_initialization_v1", ...
    "compute_external_rg_fit_stats_v1","adjusted_rand_index_v1","run_hd_multistart_top_hier_stable4_unseen_v1","select_tfkm_hier_stable4_unseen_v1"];
paths=strings(numel(req),1);
for q=1:numel(req)
    f=string(which(char(req(q))));assert(strlength(f)>0,'Missing %s.',req(q));
    assert(startsWith(f,tf,'IgnoreCase',true),'TFKM HIER4 PATH LOCK FAIL %s -> %s',req(q),f);paths(q)=f;
end
assert(contains(paths(req=="run_hd_multistart_top_hier_stable4_unseen_v1"),[filesep 'hier4' filesep],'IgnoreCase',true),'HIER4 wrapper did not resolve in vendor/hier4.');
P=table(req(:),paths,'VariableNames',{'function_name','resolved_path'});
end
function safe_rmpath_tree_local(tree)
tree=char(tree);dirs=string(strsplit(genpath(tree),pathsep));dirs=dirs(strlength(dirs)>0);cur=string(strsplit(path,pathsep));
for i=1:numel(dirs),if any(strcmpi(cur,dirs(i))),rmpath(char(dirs(i)));end,end
end
