function OUT = RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1(scenario_id,replicate,method,documents_root)
%RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1 Fit one fresh scenario x replicate x method cube.
%
% FIT PHASE ONLY. Full 15 (r,g) x 10 alpha = 150 selected candidates.
% No ARI, no truth labels, no recovery metrics are computed or printed.
% Checkpoint/resume occurs after every candidate.
%
% method: "REG_TFKM_HD" or "REG_TRKM_HD".
%
% Usage:
%   RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1(51001,1,"REG_TFKM_HD");
%   RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1(51001,1,"REG_TRKM_HD");
%   RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1("SELFTEST",[],[]);

if nargin>=1 && strcmpi(string(scenario_id),"SELFTEST")
    RG=rg_grid_local();assert(size(RG,1)==15&&isequal(RG(1,:),[1 2])&&isequal(RG(end,:),[5 6]));
    x=struct('selector_score',1.2);assert(abs(scalar_field_local(x,'selector_score',NaN)-1.2)<1e-12);
    fprintf('SELFTEST PASS: fresh unit fit contracts + 15-pair grid.\n');OUT=struct('selftest_pass',true);return
end
if nargin<3,error('Provide scenario_id, replicate, method.');end
if nargin<4,documents_root="";end
method=string(method);assert(any(method==["REG_TFKM_HD","REG_TRKM_HD"]),'Unsupported method.');
[root,outroot,study]=load_study_local(documents_root);
options=study.options;catalog=study.catalog;manifest=study.manifest;CFG=study.CFG;
assert(replicate>=1&&replicate<=CFG.n_replicates&&replicate==round(replicate),'replicate out of range.');
ic=find(double(catalog.scenario_id)==double(scenario_id),1);assert(~isempty(ic),'Unknown scenario_id.');
cfg=table2struct(catalog(ic,:));
im=find(double(manifest.scenario_id)==double(scenario_id)&double(manifest.replicate)==double(replicate),1);assert(~isempty(im),'Manifest row missing.');
data_seed=double(manifest.data_seed(im));
if method=="REG_TFKM_HD",pair_seed=double(manifest.tfkm_pair_seed(im));else,pair_seed=double(manifest.trkm_pair_seed(im));end

canon=lock_vendor_local(root,method);
req=["generate_fixed_dataset_v1","run_hd_multistart_top_v1","compute_external_rg_fit_stats_v1"];
for q=req
    f=string(which(char(q)));assert(strlength(f)>0,'Missing %s.',q);assert(startsWith(f,canon,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,f);
end

sim=generate_fixed_dataset_v1(cfg,data_seed);X=sim.X;n=size(X,1);p=size(X,2);
assert(n==cfg.n&&p==cfg.p,'Generated dimension drift.');
alpha_grid=double(CFG.alpha_grid);RG=rg_grid_local();N=150;

mshort="TFKM";if method=="REG_TRKM_HD",mshort="TRKM";end
unitdir=fullfile(outroot,'fits',char(mshort),sprintf('s%d',scenario_id),sprintf('rep%02d',replicate));
if ~isfolder(unitdir),mkdir(unitdir);end
matfile=fullfile(unitdir,'cube.mat');csvfile=fullfile(unitdir,'candidates.csv');completefile=fullfile(unitdir,'COMPLETE.txt');
if isfile(completefile)&&isfile(matfile)
    Q=load(matfile,'Candidates','unit_info');
    fprintf('%s s%d rep%02d already COMPLETE (%d/150).\n',mshort,scenario_id,replicate,height(Q.Candidates));
    OUT=struct('complete',true,'unit_dir',string(unitdir),'candidates',Q.Candidates,'unit_info',Q.unit_info);return
end
if isfile(matfile)
    Q=load(matfile,'Candidates','Labels','Trim','unit_info');Candidates=Q.Candidates;Labels=Q.Labels;Trim=Q.Trim;unit_info=Q.unit_info;
else
    Candidates=empty_candidate_table_local();Labels=zeros(n,N,'uint16');Trim=false(n,N);
    unit_info=struct('study',CFG.study_name,'scenario_id',double(scenario_id),'replicate',double(replicate), ...
        'method',method,'data_seed',data_seed,'pair_seed',pair_seed,'n',n,'p',p,'complete',false,'started_at',datetime('now'));
end
fprintf('\n============================================================\n');
fprintf('%s FRESH FIT | sid=%d rep=%02d | n=%d p=%d | existing=%d/150\n',mshort,scenario_id,replicate,n,p,height(Candidates));
fprintf('FIT PHASE ONLY: truth performance is sealed.\n');
fprintf('============================================================\n');
row_index=0;
for ia=1:numel(alpha_grid)
    a=alpha_grid(ia);
    for irg=1:size(RG,1)
        row_index=row_index+1;r=RG(irg,1);g=RG(irg,2);
        if ~isempty(Candidates)&&any(double(Candidates.row_index)==row_index),continue,end
        t0=tic;fit=run_hd_multistart_top_v1(X,r,g,a,method,options,pair_seed);runtime=toc(t0);
        assert(isfield(fit,'valid')&&logical(fit.valid),'Invalid fit sid=%d rep=%d %s (%d,%d,%.3f).',scenario_id,replicate,method,r,g,a);
        stats=compute_external_rg_fit_stats_v1(X,fit,method);
        Labels(:,row_index)=uint16(fit.labels(:));Trim(:,row_index)=logical(fit.trim(:));
        row=candidate_row_local(row_index,scenario_id,replicate,data_seed,pair_seed,n,p,method,r,g,a,stats,fit,runtime);
        Candidates=[Candidates;row]; %#ok<AGROW>
        Candidates=sortrows(Candidates,'row_index');
        unit_info.completed_candidates=height(Candidates);unit_info.last_saved_at=datetime('now');
        save(matfile,'Candidates','Labels','Trim','unit_info','-v7.3');writetable(Candidates,csvfile);
        fprintf('%3d/150 | a=%5.3f | (%d,%d) | score=% .6f | R2=% .6f | %s | %.1fs\n', ...
            row_index,a,r,g,double(row.selector_score),double(row.R2),string(row.selected_init_family),runtime);
    end
end
assert(height(Candidates)==150&&all(logical(Candidates.valid)),'Unit cube incomplete/invalid.');
unit_info.complete=true;unit_info.completed_at=datetime('now');unit_info.completed_candidates=150;
save(matfile,'Candidates','Labels','Trim','unit_info','-v7.3');writetable(Candidates,csvfile);
fid=fopen(completefile,'w');fprintf(fid,'COMPLETE\nstudy=%s\nscenario_id=%d\nreplicate=%d\nmethod=%s\n',CFG.study_name,scenario_id,replicate,method);fclose(fid);
fprintf('%s s%d rep%02d COMPLETE.\n',mshort,scenario_id,replicate);
OUT=struct('complete',true,'unit_dir',string(unitdir),'candidates',Candidates,'unit_info',unit_info);
end

function row=candidate_row_local(idx,sid,rep,ds,ps,n,p,method,r,g,a,stats,fit,runtime)
row=table(idx,double(sid),double(rep),double(ds),double(ps),double(n),double(p),method,double(r),double(g),double(a), ...
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
function RG=rg_grid_local()
RG=zeros(0,2);for r=1:5,for g=2:6,if r<g,RG(end+1,:)=[r g];end,end,end
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";
for c=cands(:)'
    o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break;end
end
assert(strlength(root)>0,'Fresh study not built. Run BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1 first.');S=load(fullfile(outroot,'study.mat'));
end
function canon=lock_vendor_local(root,method)
% Quietly enforce method-specific canonical vendor precedence.
% rmpath(genpath(tree)) warns for every directory that is not currently on
% the MATLAB path.  Remove only directories that are actually present.
tf=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));
tr=string(fullfile(root,'TRKM_TARGETED_MAX_BACKWARD_V1_1','vendor','TFKM_TRKM_UNSEEN_FINAL_v1'));
assert(isfolder(tf),'TFKM canonical vendor missing: %s',tf);assert(isfolder(tr),'TRKM canonical vendor missing: %s',tr);
if method=="REG_TFKM_HD"
    safe_rmpath_tree_local(tr);
    addpath(genpath(char(tf)),'-begin');
    canon=tf;
else
    safe_rmpath_tree_local(tf);
    addpath(genpath(char(tr)),'-begin');
    canon=tr;
end
rehash path;
end

function safe_rmpath_tree_local(tree)
% Remove only exact path entries belonging to TREE; emit no benign warnings.
tree=char(tree);
dirs=string(strsplit(genpath(tree),pathsep));
dirs=dirs(strlength(dirs)>0);
cur=string(strsplit(path,pathsep));
for i=1:numel(dirs)
    if any(strcmpi(cur,dirs(i)))
        rmpath(char(dirs(i)));
    end
end
end
