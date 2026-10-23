function report = run_multistart_basin_audit_v1(output_root,mode)
%RUN_MULTISTART_BASIN_AUDIT_V1 Exact deterministic refit audit of v11 basins.
%
% Default mode PRIMARY:
%   - TFKM current-selected pools, clean + bridge, reps 1:8
%   - TRKM current-selected pools, clean + bridge, reps 1:8
%
% Other modes:
%   TFKM_SELECTED       : 16 TFKM current-selected clean/bridge pools
%   TFKM_SELECTED_TRUE  : TFKM selected + true-fixed controls (32 pools)
%   TRKM_SELECTED       : 16 TRKM current-selected clean/bridge pools
%   FULL_CLEAN_BRIDGE   : TFKM selected + true-fixed, TRKM selected (48 pools)
%   FAILURES            : current-selected clean/bridge pools with stored ARI < .95
%
% The audit never changes the frozen v11 selection. Truth is used only after
% each terminal fit has been produced and after truth-blind challenger scores
% have been computed.

if nargin<2 || strlength(string(mode))==0, mode="PRIMARY"; end
mode=upper(string(mode));
pkg=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(pkg),'-begin');
output_root=char(string(output_root));
assert(isfolder(output_root),'Output root does not exist: %s',output_root);
S=load(fullfile(output_root,'study.mat'));
assert(isfield(S,'study')&&isfield(S,'options')&&isfield(S,'catalog'),'study.mat is incomplete.');
study=S.study; options=S.options; catalog=S.catalog;

summary_file=fullfile(output_root,'summary','all_joint_candidates.csv');
if ~isfile(summary_file)
    fprintf('Merged candidate cube not found; attempting frozen merge...\n');
    merge_final_fresh_validation_v1(output_root);
end
assert(isfile(summary_file),'Missing merged candidate cube: %s',summary_file);
C=readtable(summary_file,'TextType','string');

scenario_keep=[23001 23004]; % clean and bridge: primary basin stress test.
reps=double(study.validation_reps(:));
targets=build_targets(C,catalog,study,scenario_keep,reps,mode);
assert(~isempty(targets),'No audit targets generated for mode %s.',mode);

audit_root=fullfile(output_root,'multistart_basin_audit_v1');
pool_root=fullfile(audit_root,'pools'); summary_root=fullfile(audit_root,'summary');
if ~isfolder(pool_root),mkdir(pool_root);end
if ~isfolder(summary_root),mkdir(summary_root);end
writetable(targets,fullfile(summary_root,'audit_targets.csv'));

fprintf('\n============================================================\n');
fprintf('MULTISTART BASIN AUDIT v1 | mode=%s | pools=%d\n',mode,height(targets));
fprintf('Primary scenarios: 23001 clean, 23004 bridge.\n');
fprintf('Frozen estimator is NOT modified.\n');
fprintf('============================================================\n\n');

last_sid=NaN; last_rep=NaN; sim=[];
for i=1:height(targets)
    t=targets(i,:); key=pool_key(t);
    starts_file=fullfile(pool_root,key+"__starts.csv");
    rules_file=fullfile(pool_root,key+"__rules.csv");
    fam_file=fullfile(pool_root,key+"__families.csv");
    mat_file=fullfile(pool_root,key+"__compact.mat");
    if isfile(rules_file) && isfile(starts_file) && isfile(fam_file)
        fprintf('[%2d/%2d] SKIP completed %s\n',i,height(targets),key);
        continue;
    end

    sid=double(t.scenario_id); rep=double(t.replicate); method=string(t.method);
    cfg=table2struct(catalog(double(catalog.scenario_id)==sid,:));
    if sid~=last_sid || rep~=last_rep || isempty(sim)
        data_seed=double(study.validation_base_seed)+sid*100000+rep;
        sim=generate_fixed_dataset_v1(cfg,data_seed);
        last_sid=sid;last_rep=rep;
    else
        data_seed=double(study.validation_base_seed)+sid*100000+rep;
    end
    if method=="REG_TFKM_HD", pair_seed=data_seed+110000; else, pair_seed=data_seed+220000; end

    fprintf('[%2d/%2d] %s | sid=%d rep=%d %s | (%d,%d,%.3f) ...\n', ...
        i,height(targets),string(t.target_type),sid,rep,method,double(t.r),double(t.g),double(t.alpha));
    A=run_hd_multistart_pool_audit_v1(sim,double(t.r),double(t.g),double(t.alpha),method,options,pair_seed);

    % Reproduction gate: compare frozen-current refit against stored v11 row.
    cur=A.rules(A.rules.rule==current_rule_name(method),:);
    assert(height(cur)==1,'Current rule row missing.');
    repcheck=reproduction_check(t,cur);

    % Add pool metadata and truth-open pool references to all outputs.
    starts=add_start_metadata(A.starts,t);
    families=add_family_metadata(A.family_summary,t);
    rules=add_rule_metadata(A.rules,t,A,repcheck);

    writetable(starts,starts_file);
    writetable(families,fam_file);
    writetable(rules,rules_file);
    pool=A.pool; pairwise_similarity=A.pairwise_similarity; %#ok<NASGU>
    save(mat_file,'pool','pairwise_similarity','repcheck','-v7.3');

    fprintf('          reproduction=%d | current ARI=%.3f | best ALL=%.3f | best SIGN=%.3f | oracle family=%s\n', ...
        repcheck.pass,cur.selected_ari,A.best_all_ari,A.best_sign_ari,A.oracle_best_family);
    if ~repcheck.pass
        error('basinAudit:ReproductionFailed', ...
            'Frozen current winner did not reproduce for %s. Audit stopped before interpretation.',key);
    end
end

report=summarize_multistart_basin_audit_v1(audit_root);
fprintf('\nAudit summary written to:\n  %s\n',summary_root);
end

function T=build_targets(C,catalog,study,sids,reps,mode)
T=table();
for sid=sids
    cfg=catalog(double(catalog.scenario_id)==sid,:); assert(height(cfg)==1,'Catalog row missing.');
    for rep=reps(:)'
        for method=["REG_TFKM_HD","REG_TRKM_HD"]
            if ~include_method(mode,method),continue;end
            U=C(double(C.scenario_id)==sid & double(C.replicate)==rep & string(C.method)==method,:);
            assert(height(U)==150,'Expected 150 cube rows for sid=%d rep=%d %s; got %d.',sid,rep,method,height(U));
            sel=select_final_pipeline_v1(U,method); sr=sel.selected_row;
            if mode=="FAILURES" && double(sr.ari_inlier)>=.95,continue;end
            T=[T; target_row(sid,rep,method,"CURRENT_SELECTED",sel.r,sel.g,sel.alpha,sr)]; %#ok<AGROW>
            if include_true_control(mode,method)
                r0=double(cfg.r0);g0=double(cfg.c0);a0=double(cfg.alpha0);
                q=U(double(U.r)==r0 & double(U.g)==g0 & abs(double(U.alpha)-a0)<1e-12,:);
                assert(height(q)==1,'True-fixed row absent.');
                T=[T; target_row(sid,rep,method,"TRUE_FIXED",r0,g0,a0,q)]; %#ok<AGROW>
            end
        end
    end
end
T=sortrows(T,{'scenario_id','replicate','method','target_type'});
end

function tf=include_method(mode,method)
switch mode
    case {"PRIMARY","FULL_CLEAN_BRIDGE","FAILURES"},tf=true;
    case {"TFKM_SELECTED","TFKM_SELECTED_TRUE"},tf=(method=="REG_TFKM_HD");
    case "TRKM_SELECTED",tf=(method=="REG_TRKM_HD");
    otherwise,error('basinAudit:Mode','Unknown audit mode %s.',mode);
end
end
function tf=include_true_control(mode,method)
tf=(method=="REG_TFKM_HD") && ismember(mode,["TFKM_SELECTED_TRUE","FULL_CLEAN_BRIDGE"]);
end

function row=target_row(sid,rep,method,target_type,r,g,a,q)
row=table(double(sid),double(rep),string(method),string(target_type),double(r),double(g),double(a), ...
    double(q.ari_inlier),string(q.selected_init_family),double(q.selected_objective_rank), ...
    double(q.selector_score),double(q.selector_cv_score),double(q.objective), ...
    'VariableNames',{'scenario_id','replicate','method','target_type','r','g','alpha', ...
    'stored_current_ari','stored_selected_family','stored_objective_rank','stored_selector_score', ...
    'stored_selector_cv_score','stored_objective'});
end

function key=pool_key(t)
key=sprintf('sid%d_rep%d_%s_%s_r%d_g%d_a%0.3f',double(t.scenario_id),double(t.replicate), ...
    char(erase(string(t.method),"REG_")),char(string(t.target_type)),double(t.r),double(t.g),double(t.alpha));
key=string(strrep(key,'.','p'));
end

function name=current_rule_name(method)
if method=="REG_TFKM_HD",name="CURRENT_SIGN_ALL_STABLE";else,name="CURRENT_RV_CV";end
end

function R=reproduction_check(t,cur)
ari_ok=close_num(double(t.stored_current_ari),double(cur.selected_ari),1e-10);
fam_ok=string(t.stored_selected_family)==string(cur.chosen_family);
rank_ok=close_num(double(t.stored_objective_rank),double(cur.objective_rank),1e-12);
obj_ok=close_num(double(t.stored_objective),double(cur.objective),1e-8);
if string(t.method)=="REG_TFKM_HD"
    score_ok=close_num(double(t.stored_selector_score),double(cur.rule_score),1e-10);
else
    score_ok=close_num(double(t.stored_selector_cv_score),double(cur.rule_score),1e-10);
end
pass=ari_ok&&fam_ok&&rank_ok&&obj_ok&&score_ok;
R=struct('pass',pass,'ari_ok',ari_ok,'family_ok',fam_ok,'objective_rank_ok',rank_ok,'objective_ok',obj_ok,'score_ok',score_ok);
end
function tf=close_num(a,b,tol)
if isnan(a)&&isnan(b),tf=true;elseif ~isfinite(a)||~isfinite(b),tf=(a==b);else,tf=abs(a-b)<=tol*max(1,max(abs([a b])));end
end

function T=add_start_metadata(T,t)
n=height(T);T.scenario_id=repmat(double(t.scenario_id),n,1);T.replicate=repmat(double(t.replicate),n,1);
T.method=repmat(string(t.method),n,1);T.target_type=repmat(string(t.target_type),n,1);
T.r=repmat(double(t.r),n,1);T.g=repmat(double(t.g),n,1);T.alpha=repmat(double(t.alpha),n,1);
T=movevars(T,{'scenario_id','replicate','method','target_type','r','g','alpha'},'Before',1);
end
function T=add_family_metadata(T,t)
n=height(T);T.scenario_id=repmat(double(t.scenario_id),n,1);T.replicate=repmat(double(t.replicate),n,1);
T.method=repmat(string(t.method),n,1);T.target_type=repmat(string(t.target_type),n,1);
T.r=repmat(double(t.r),n,1);T.g=repmat(double(t.g),n,1);T.alpha=repmat(double(t.alpha),n,1);
T=movevars(T,{'scenario_id','replicate','method','target_type','r','g','alpha'},'Before',1);
end
function T=add_rule_metadata(T,t,A,R)
n=height(T);T.scenario_id=repmat(double(t.scenario_id),n,1);T.replicate=repmat(double(t.replicate),n,1);
T.method=repmat(string(t.method),n,1);T.target_type=repmat(string(t.target_type),n,1);
T.r=repmat(double(t.r),n,1);T.g=repmat(double(t.g),n,1);T.alpha=repmat(double(t.alpha),n,1);
T.stored_current_ari=repmat(double(t.stored_current_ari),n,1);
T.best_all_ari=repmat(A.best_all_ari,n,1);T.best_sign_ari=repmat(A.best_sign_ari,n,1);T.best_non_sign_ari=repmat(A.best_non_sign_ari,n,1);
T.basin_gap_all_vs_sign=T.best_all_ari-T.best_sign_ari;T.selector_regret=T.best_all_ari-T.selected_ari;
T.oracle_best_family=repmat(string(A.oracle_best_family),n,1);T.reproduction_pass=repmat(logical(R.pass),n,1);
T=movevars(T,{'scenario_id','replicate','method','target_type','r','g','alpha'},'Before',1);
end
