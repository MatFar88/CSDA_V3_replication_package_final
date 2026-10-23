function OUT = RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1(documents_root)
%RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1
% Broad truth-open forensic diagnostic for the frozen fresh study.
%
% IMPORTANT SCIENTIFIC STATUS
%   * Truth has already been opened by fresh_truth_open_v1.
%   * This runner NEVER fits or refits an estimator.
%   * It NEVER changes selector/action/alpha/anchor/backward rules.
%   * It reads only frozen predictions, frozen audit objects, frozen candidate
%     cubes, and the fixed DGP truth for POSTHOC diagnosis.
%   * All oracle quantities are explicitly labelled POSTHOC and must not be
%     interpreted as validation performance or used to retune V1.
%
% Main goals
%   1. Decompose action, alpha, anchor, backward, ambiguity and local-fit failures.
%   2. Diagnose TFKM leverage non-emission by the frozen J/Smed/DS gates.
%   3. Diagnose TRKM false positives and separate radial vs archived sources.
%   4. Diagnose structural failure by p/n, rho, delta and (r0,g0).
%   5. Ask whether the frozen selected fit AT THE TRUE STATE was itself viable.
%   6. Compute a no-refit POSTHOC oracle over the 15 frozen states at alpha0.
%
% Usage
%   RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1;
%   RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1("C:\Users\name\Documents");
%   RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1("SELFTEST");

if nargin>=1 && strcmpi(string(documents_root),"SELFTEST")
    assert(reachable_local([4 6],[3 4]));
    assert(reachable_local([3 5],[3 4]));
    assert(~reachable_local([2 3],[3 4]));
    assert(isequal(predecessor_local([3 6]),[3 5]));
    assert(isequal(predecessor_local([3 4]),[2 3]));
    assert(abs(f1_local(1,1,1)-1)<=eps);
    assert(abs(f1_local(1,2,1)-2/3)<=10*eps);
    M=[3 4;3 5;2 4];
    Mp=componentwise_maximal_local(M);
    assert(isequal(Mp,[3 5]));
    fprintf('SELFTEST PASS: posthoc diagnostic contracts.\n');
    OUT=struct('selftest_pass',true);return
end
if nargin<1,documents_root="";end
[root,outroot,S]=load_study_local(documents_root);
catalog=S.catalog;manifest=S.manifest;
assert(height(catalog)==41 && height(manifest)==1230,'Fresh study dimensions drifted.');

truthroot=fullfile(outroot,'fresh_truth_open_v1');
assert(isfile(fullfile(truthroot,'TRUTH_OPENED_AT.txt')),'Truth-open marker missing. Run first truth-open evaluation first.');
unitTruth=fullfile(truthroot,'FRESH_TRUTH_UNIT_RESULTS.csv');
assert(isfile(unitTruth),'Missing FRESH_TRUTH_UNIT_RESULTS.csv.');
R=readtable(unitTruth,'TextType','string');
assert(height(R)==2460,'Expected 2460 truth-open unit rows.');

% Frozen pretruth seal is still used as the authoritative prediction table.
seal=fullfile(outroot,'FROZEN_PRETRUTH_20261010');
assert(isfolder(seal),'Pretruth frozen snapshot missing.');
PT=readtable(fullfile(seal,'TFKM_HIER4_BLIND_FINAL_PREDICTIONS_FROZEN.csv'),'TextType','string');
PR=readtable(fullfile(seal,'TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv'),'TextType','string');
assert(height(PT)==1230 && height(PR)==1230,'Frozen prediction count mismatch.');

lock_truth_vendor_local(root);
assert(exist('generate_fixed_dataset_v1','file')==2,'generate_fixed_dataset_v1 missing.');
assert(exist('adjusted_rand_index_v1','file')==2,'adjusted_rand_index_v1 missing.');

hierRoot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
assert(isfolder(hierRoot),'Corrected HIER4 rerun root missing.');
tfBlind=fullfile(hierRoot,'blind_finalization_hier4_v1','TFKM');
trBlind=fullfile(outroot,'blind_finalization_v1','TRKM');
tfCube=fullfile(hierRoot,'fits','TFKM');
trCube=fullfile(outroot,'fits','TRKM');
for p=[string(tfBlind),string(trBlind),string(tfCube),string(trCube)]
    assert(isfolder(p),'Required frozen diagnostic source missing: %s',p);
end

fprintf('\n============================================================\n');
fprintf('TFKM/TRKM FRESH POSTHOC DIAGNOSTIC V1\n');
fprintf('TRUTH OPEN | NO REFITS | NO RETUNING | FROZEN V1 ONLY\n');
fprintf('============================================================\n');

D=repmat(empty_diag_struct_local(),2460,1);
pos=0;
for im=1:2
    if im==1
        method="TFKM";P=sortrows(PT,{'scenario_id','replicate'});blindBase=tfBlind;cubeBase=tfCube;
    else
        method="TRKM";P=sortrows(PR,{'scenario_id','replicate'});blindBase=trBlind;cubeBase=trCube;
    end
    for i=1:height(P)
        pos=pos+1;
        sid=double(P.scenario_id(i));rep=double(P.replicate(i));
        ic=find(double(catalog.scenario_id)==sid,1);assert(~isempty(ic),'Catalog sid missing.');
        imf=find(double(manifest.scenario_id)==sid & double(manifest.replicate)==rep,1);assert(~isempty(imf),'Manifest key missing.');
        c=catalog(ic,:);cfg=table2struct(c);
        sim=generate_fixed_dataset_v1(cfg,double(manifest.data_seed(imf)));
        [ytrue,out]=truth_local(sim);in=~out;

        rt=R(string(R.method)==method & double(R.scenario_id)==sid & double(R.replicate)==rep,:);
        assert(height(rt)==1,'Truth result key mismatch.');
        pr=P(i,:);
        bf=fullfile(blindBase,sprintf('s%d',sid),sprintf('rep%02d',rep),'blind_final.mat');
        assert(isfile(bf),'Missing blind_final.mat: %s',bf);
        Z=load(bf,'Prediction','Edges','Audit');
        assert(height(Z.Prediction)==1,'Prediction row missing in blind_final.');
        E=Z.Edges;A=Z.Audit;

        cf=fullfile(cubeBase,sprintf('s%d',sid),sprintf('rep%02d',rep),'cube.mat');
        assert(isfile(cf),'Missing cube.mat: %s',cf);
        Q=load(cf,'Candidates','Labels','Trim');C=Q.Candidates;L=Q.Labels;T=Q.Trim;
        assert(height(C)==150 && size(L,2)==150 && size(T,2)==150,'Candidate cube shape drift.');

        d=empty_diag_struct_local();
        d.method=method;d.scenario_id=sid;d.replicate=rep;d.panel=string(c.panel);
        d.n=double(c.n);d.p=double(c.p);d.p_over_n=double(c.p_over_n);d.rho=double(c.rho);d.delta=double(c.delta);
        d.r0_POSTHOC=double(c.r0);d.g0_POSTHOC=double(c.g0);d.contamination_POSTHOC=string(c.contamination);d.alpha0_POSTHOC=double(c.alpha0);
        d.leverage_POSTHOC=(d.contamination_POSTHOC=="leverage");

        % First truth-open outcomes (copied, not recomputed).
        d.alpha_emitted=logical(rt.alpha_emitted);d.action_correct_POSTHOC=logical(rt.action_correct_POSTHOC);
        d.alpha_hat=double(rt.alpha_hat);d.alpha_exact_POSTHOC=logical(rt.alpha_exact_POSTHOC);
        d.terminal_unique=logical(rt.terminal_unique);d.terminal_r=double(rt.terminal_r);d.terminal_g=double(rt.terminal_g);
        d.structure_exact_POSTHOC=logical(rt.structure_exact_POSTHOC);d.ambiguous_output_POSTHOC=logical(rt.ambiguous_output_POSTHOC);
        d.joint_correct_POSTHOC=logical(rt.joint_correct_POSTHOC);d.final_fit_materialized=logical(rt.final_fit_materialized);
        d.final_inlier_ARI_POSTHOC=double(rt.final_inlier_ARI_POSTHOC);d.final_trim_F1_POSTHOC=double(rt.final_trim_F1_POSTHOC);
        d.action_status=string(pr.action_status);d.terminal_status=string(pr.terminal_status);
        d.n_anchors=double(pr.n_anchors);d.N_edges=double(pr.N_edges);d.N_accepts=double(pr.N_accepts);
        d.dr_terminal_POSTHOC=d.terminal_r-d.r0_POSTHOC;d.dg_terminal_POSTHOC=d.terminal_g-d.g0_POSTHOC;
        if d.terminal_unique && isfinite(d.terminal_r) && isfinite(d.terminal_g)
            d.complexity_error_POSTHOC=(d.terminal_r+d.terminal_g)-(d.r0_POSTHOC+d.g0_POSTHOC);
        end

        % Frozen selected fit at the exact truth state; no refit.
        [d.truth_state_inlier_ARI_POSTHOC,d.truth_state_trim_F1_POSTHOC,d.truth_state_family_POSTHOC,d.truth_state_score_POSTHOC] = ...
            state_metrics_local(C,L,T,d.r0_POSTHOC,d.g0_POSTHOC,d.alpha0_POSTHOC,ytrue,out);
        d.truth_state_ARI_ge_090_POSTHOC=isfinite(d.truth_state_inlier_ARI_POSTHOC) && d.truth_state_inlier_ARI_POSTHOC>=.90;

        % POSTHOC no-refit oracle among the 15 already-frozen selected states at alpha0.
        [d.best_alpha0_inlier_ARI_POSTHOC,d.best_alpha0_r_POSTHOC,d.best_alpha0_g_POSTHOC,d.best_alpha0_family_POSTHOC] = ...
            best_alpha0_state_local(C,L,d.alpha0_POSTHOC,ytrue,in);
        d.best_alpha0_ARI_ge_090_POSTHOC=isfinite(d.best_alpha0_inlier_ARI_POSTHOC) && d.best_alpha0_inlier_ARI_POSTHOC>=.90;
        d.truth_state_ARI_regret_POSTHOC=d.best_alpha0_inlier_ARI_POSTHOC-d.truth_state_inlier_ARI_POSTHOC;
        d.best_alpha0_structure_is_truth_POSTHOC=(d.best_alpha0_r_POSTHOC==d.r0_POSTHOC && d.best_alpha0_g_POSTHOC==d.g0_POSTHOC);

        % Method-specific anchor/action evidence.
        if method=="TFKM"
            W=double(A.per_alpha);Mplus=double(A.Mplus);AA=A.anchor_action;
            d.per_alpha_structure_exact_count_POSTHOC=sum(W(:,2)==d.r0_POSTHOC & W(:,3)==d.g0_POSTHOC);
            d.anchor_truth_in_set_POSTHOC=contains_row_local(Mplus,[d.r0_POSTHOC d.g0_POSTHOC]);
            [d.any_anchor_reaches_truth_POSTHOC,d.all_anchors_reach_truth_POSTHOC]=reachability_set_local(Mplus,[d.r0_POSTHOC d.g0_POSTHOC]);
            d.tfkm_Jmax_min_POSTHOC=safe_min_local(double(AA.Jmax));d.tfkm_Jmax_max_POSTHOC=safe_max_local(double(AA.Jmax));
            d.tfkm_Smed_min_POSTHOC=safe_min_local(double(AA.Smed));d.tfkm_Smed_max_POSTHOC=safe_max_local(double(AA.Smed));
            d.tfkm_DS_min_POSTHOC=safe_min_local(double(AA.DS));d.tfkm_DS_max_POSTHOC=safe_max_local(double(AA.DS));
            d.tfkm_action_block_reason_POSTHOC=tfkm_action_block_local(AA,logical(pr.action_ambiguous));
        else
            AW=double(A.archived_per_alpha);d.per_alpha_structure_exact_count_POSTHOC=sum(AW(:,2)==d.r0_POSTHOC & AW(:,3)==d.g0_POSTHOC);
            d.trkm_archived_modal_exact_POSTHOC=isequal(double(A.archived_modal(:)'),[d.r0_POSTHOC d.g0_POSTHOC]);
            d.trkm_archived_emit=logical(pr.archived_emit);d.trkm_archived_alpha_hat=double(pr.archived_alpha_hat);d.trkm_archived_delta_RV=double(pr.archived_delta_RV);
            d.trkm_radial_emit=logical(pr.radial_emit);d.trkm_radial_alpha_hat=double(pr.radial_alpha_hat);d.trkm_radial_high_weight=double(pr.radial_high_weight);
            d.trkm_radial_delta_BIC=double(pr.radial_delta_BIC);d.trkm_radial_separation_sd=double(pr.radial_separation_sd);
            if d.alpha_emitted
                RW=double(A.raw_per_alpha);anchors=componentwise_maximal_local(RW(RW(:,1)>=d.alpha_hat-1e-12,2:3));
                d.trkm_raw_tail_structure_exact_count_POSTHOC=sum(RW(RW(:,1)>=d.alpha_hat-1e-12,2)==d.r0_POSTHOC & RW(RW(:,1)>=d.alpha_hat-1e-12,3)==d.g0_POSTHOC);
            else
                anchors=double(A.archived_modal(:)');
            end
            d.anchor_truth_in_set_POSTHOC=contains_row_local(anchors,[d.r0_POSTHOC d.g0_POSTHOC]);
            [d.any_anchor_reaches_truth_POSTHOC,d.all_anchors_reach_truth_POSTHOC]=reachability_set_local(anchors,[d.r0_POSTHOC d.g0_POSTHOC]);
            if ~d.leverage_POSTHOC && d.alpha_emitted
                if d.trkm_radial_emit,d.trkm_false_positive_source_POSTHOC="RADIAL_GMM";elseif d.trkm_archived_emit,d.trkm_false_positive_source_POSTHOC="ARCHIVED_FALLBACK";else,d.trkm_false_positive_source_POSTHOC="UNKNOWN";end
            end
        end

        % Backward-path evidence already frozen before truth opening.
        [d.n_stop_edges_POSTHOC,d.min_selected_edge_ARI_POSTHOC,d.rejected_edge_into_truth_POSTHOC, ...
            d.accepted_edge_below_truth_POSTHOC,d.premature_stop_above_truth_POSTHOC] = ...
            edge_truth_diagnostics_local(E,[d.r0_POSTHOC d.g0_POSTHOC]);

        d.structure_failure_stage_POSTHOC=structure_failure_stage_local(d);
        d.primary_failure_stage_POSTHOC=primary_failure_stage_local(d);
        D(pos)=d;

        if mod(pos,100)==0,fprintf('POSTHOC progress %d/2460\n',pos);end
    end
end
Tdiag=struct2table(D);
Tdiag=sortrows(Tdiag,{'method','scenario_id','replicate'});
assert(height(Tdiag)==2460,'Diagnostic table size mismatch.');

outdiag=fullfile(outroot,'fresh_posthoc_diagnostic_v1');
if ~isfolder(outdiag),mkdir(outdiag);end
writetable(Tdiag,fullfile(outdiag,'POSTHOC_UNIT_DIAGNOSTIC.csv'));

S_overall=summarize_local(Tdiag,{'method'});
S_panel_cont=summarize_local(Tdiag,{'method','panel','contamination_POSTHOC'});
S_main=summarize_local(Tdiag(Tdiag.panel=="MAIN",:),{'method','contamination_POSTHOC','p_over_n','rho'});
S_weak=summarize_local(Tdiag(Tdiag.panel=="WEAK_DELTA7",:),{'method','p_over_n','rho'});
S_structure=summarize_local(Tdiag(Tdiag.panel=="STRUCTURE",:),{'method','contamination_POSTHOC','r0_POSTHOC','g0_POSTHOC'});
S_primary=count_strings_local(Tdiag,{'method'},'primary_failure_stage_POSTHOC');
S_structfail=count_strings_local(Tdiag,{'method','panel'},'structure_failure_stage_POSTHOC');
S_terminal=terminal_confusion_local(Tdiag);
S_fit=summarize_fit_local(Tdiag,{'method','panel','contamination_POSTHOC'});
S_tfkm_block=count_strings_local(Tdiag(Tdiag.method=="TFKM" & Tdiag.leverage_POSTHOC,:),{'panel','p_over_n','rho'},'tfkm_action_block_reason_POSTHOC');
S_trkm_fp=count_strings_local(Tdiag(Tdiag.method=="TRKM" & ~Tdiag.leverage_POSTHOC & Tdiag.alpha_emitted,:),{'panel','contamination_POSTHOC'},'trkm_false_positive_source_POSTHOC');
S_action_quant=action_quantiles_local(Tdiag);

writetable(S_overall,fullfile(outdiag,'SUMMARY_OVERALL_FORENSIC.csv'));
writetable(S_panel_cont,fullfile(outdiag,'SUMMARY_PANEL_CONTAMINATION.csv'));
writetable(S_main,fullfile(outdiag,'SUMMARY_MAIN_GEOMETRY.csv'));
writetable(S_weak,fullfile(outdiag,'SUMMARY_WEAK_GEOMETRY.csv'));
writetable(S_structure,fullfile(outdiag,'SUMMARY_STRUCTURE_PANEL.csv'));
writetable(S_primary,fullfile(outdiag,'FAILURE_DECOMPOSITION_PRIMARY.csv'));
writetable(S_structfail,fullfile(outdiag,'STRUCTURE_FAILURE_DECOMPOSITION.csv'));
writetable(S_terminal,fullfile(outdiag,'TERMINAL_CONFUSION.csv'));
writetable(S_fit,fullfile(outdiag,'TRUE_STATE_FIT_DIAGNOSTIC.csv'));
writetable(S_tfkm_block,fullfile(outdiag,'TFKM_LEVERAGE_ACTION_BLOCKS.csv'));
writetable(S_trkm_fp,fullfile(outdiag,'TRKM_NONLEVERAGE_FALSE_POSITIVE_SOURCE.csv'));
writetable(S_action_quant,fullfile(outdiag,'ACTION_FEATURE_QUANTILES.csv'));

write_report_local(fullfile(outdiag,'POSTHOC_DIAGNOSTIC_REPORT.txt'),S_overall,S_primary,S_structfail,S_fit,S_tfkm_block,S_trkm_fp,S_main,S_weak,S_structure,outdiag);

fprintf('\n================ FORENSIC OVERALL ================\n');disp(S_overall);
fprintf('\n========== PRIMARY FAILURE DECOMPOSITION ==========\n');disp(S_primary);
fprintf('\n========== TRUE-STATE FIT DIAGNOSTIC ==============\n');disp(S_fit);
fprintf('\n========== TFKM LEVERAGE ACTION BLOCKS ============\n');disp(S_tfkm_block);
fprintf('\n========== TRKM NONLEVERAGE FP SOURCE =============\n');disp(S_trkm_fp);
fprintf('\nPOSTHOC DIAGNOSTIC COMPLETE.\nOutput: %s\n',outdiag);
fprintf('NO V1 rule has been changed or refit. These outputs are diagnostic only.\n');

OUT=struct('unit',Tdiag,'overall',S_overall,'panel_contamination',S_panel_cont,'main_geometry',S_main, ...
    'weak_geometry',S_weak,'structure_panel',S_structure,'primary_failure',S_primary, ...
    'structure_failure',S_structfail,'terminal_confusion',S_terminal,'fit_diagnostic',S_fit, ...
    'tfkm_action_blocks',S_tfkm_block,'trkm_false_positive_source',S_trkm_fp,'action_quantiles',S_action_quant, ...
    'output_root',string(outdiag));
end

% =========================================================================
% UNIT DIAGNOSTICS
% =========================================================================
function d=empty_diag_struct_local()
d=struct( ...
'method',"",'scenario_id',NaN,'replicate',NaN,'panel',"",'n',NaN,'p',NaN,'p_over_n',NaN,'rho',NaN,'delta',NaN, ...
'r0_POSTHOC',NaN,'g0_POSTHOC',NaN,'contamination_POSTHOC',"",'alpha0_POSTHOC',NaN,'leverage_POSTHOC',false, ...
'alpha_emitted',false,'action_correct_POSTHOC',false,'alpha_hat',NaN,'alpha_exact_POSTHOC',false, ...
'terminal_unique',false,'terminal_r',NaN,'terminal_g',NaN,'structure_exact_POSTHOC',false,'ambiguous_output_POSTHOC',false, ...
'joint_correct_POSTHOC',false,'final_fit_materialized',false,'final_inlier_ARI_POSTHOC',NaN,'final_trim_F1_POSTHOC',NaN, ...
'action_status',"",'terminal_status',"",'n_anchors',NaN,'N_edges',NaN,'N_accepts',NaN,'dr_terminal_POSTHOC',NaN,'dg_terminal_POSTHOC',NaN,'complexity_error_POSTHOC',NaN, ...
'truth_state_inlier_ARI_POSTHOC',NaN,'truth_state_trim_F1_POSTHOC',NaN,'truth_state_family_POSTHOC',"",'truth_state_score_POSTHOC',NaN,'truth_state_ARI_ge_090_POSTHOC',false, ...
'best_alpha0_inlier_ARI_POSTHOC',NaN,'best_alpha0_r_POSTHOC',NaN,'best_alpha0_g_POSTHOC',NaN,'best_alpha0_family_POSTHOC',"",'best_alpha0_ARI_ge_090_POSTHOC',false, ...
'truth_state_ARI_regret_POSTHOC',NaN,'best_alpha0_structure_is_truth_POSTHOC',false, ...
'per_alpha_structure_exact_count_POSTHOC',NaN,'anchor_truth_in_set_POSTHOC',false,'any_anchor_reaches_truth_POSTHOC',false,'all_anchors_reach_truth_POSTHOC',false, ...
'tfkm_Jmax_min_POSTHOC',NaN,'tfkm_Jmax_max_POSTHOC',NaN,'tfkm_Smed_min_POSTHOC',NaN,'tfkm_Smed_max_POSTHOC',NaN,'tfkm_DS_min_POSTHOC',NaN,'tfkm_DS_max_POSTHOC',NaN,'tfkm_action_block_reason_POSTHOC',"", ...
'trkm_archived_modal_exact_POSTHOC',false,'trkm_archived_emit',false,'trkm_archived_alpha_hat',NaN,'trkm_archived_delta_RV',NaN, ...
'trkm_radial_emit',false,'trkm_radial_alpha_hat',NaN,'trkm_radial_high_weight',NaN,'trkm_radial_delta_BIC',NaN,'trkm_radial_separation_sd',NaN, ...
'trkm_raw_tail_structure_exact_count_POSTHOC',NaN,'trkm_false_positive_source_POSTHOC',"", ...
'n_stop_edges_POSTHOC',0,'min_selected_edge_ARI_POSTHOC',NaN,'rejected_edge_into_truth_POSTHOC',false,'accepted_edge_below_truth_POSTHOC',false,'premature_stop_above_truth_POSTHOC',false, ...
'structure_failure_stage_POSTHOC',"",'primary_failure_stage_POSTHOC',"");
end

function [ari,f1,fam,score]=state_metrics_local(C,L,T,r,g,a,y,out)
q=find(double(C.r)==r & double(C.g)==g & abs(double(C.alpha)-a)<=1e-12,1);assert(~isempty(q),'Truth state absent from cube.');
idx=q;if ismember('row_index',C.Properties.VariableNames),idx=double(C.row_index(q));end
lab=double(L(:,idx));tr=logical(T(:,idx));in=~out;ari=NaN;if sum(in)>1,ari=adjusted_rand_index_v1(y(in),lab(in));end
f1=NaN;if any(out),tp=sum(tr&out);f1=f1_local(tp,sum(tr),sum(out));end
fam="";if ismember('selected_init_family',C.Properties.VariableNames),fam=string(C.selected_init_family(q));end
score=NaN;if ismember('selector_score',C.Properties.VariableNames),score=double(C.selector_score(q));end
end

function [best,r,g,fam]=best_alpha0_state_local(C,L,a,y,in)
ids=find(abs(double(C.alpha)-a)<=1e-12);assert(numel(ids)==15,'Expected 15 states at alpha0.');best=-Inf;r=NaN;g=NaN;fam="";
for z=ids(:)'
    idx=z;if ismember('row_index',C.Properties.VariableNames),idx=double(C.row_index(z));end
    lab=double(L(:,idx));v=adjusted_rand_index_v1(y(in),lab(in));if ~isfinite(v),continue,end
    if v>best+1e-12
        best=v;r=double(C.r(z));g=double(C.g(z));if ismember('selected_init_family',C.Properties.VariableNames),fam=string(C.selected_init_family(z));end
    elseif abs(v-best)<=1e-12
        keyNew=[double(C.r(z))+double(C.g(z)),double(C.g(z)),double(C.r(z))];keyOld=[r+g,g,r];
        if lex_less_local(keyNew,keyOld),r=double(C.r(z));g=double(C.g(z));if ismember('selected_init_family',C.Properties.VariableNames),fam=string(C.selected_init_family(z));end,end
    end
end
if best==-Inf,best=NaN;end
end

function tf=lex_less_local(a,b)
tf=false;for k=1:numel(a),if a(k)<b(k),tf=true;return;elseif a(k)>b(k),return;end,end
end

function s=tfkm_action_block_local(A,amb)
if amb,s="ANCHOR_OR_ALPHA_DISAGREEMENT";return,end
if all(logical(A.emit)),s="EMIT";return,end
jq=double(A.Jmax)>=.50-1e-12;sig=(double(A.Smed)<.44-1e-12)|(double(A.DS)>=.12-1e-12);
if all(~jq),s="J_BELOW_050_ALL_ANCHORS";
elseif all(jq & ~sig),s="J_OK_SIGNAL_BELOW_CUTOFF_ALL";
elseif all(~logical(A.emit)),s="MIXED_BLOCK_NO_EMIT";
else,s="OTHER";
end
end

function [nstop,minari,rejectInto,acceptBelow,premature]=edge_truth_diagnostics_local(E,truth)
nstop=0;minari=NaN;rejectInto=false;acceptBelow=false;premature=false;if isempty(E)||height(E)==0,return,end
if ismember('selected_edge_ARI',E.Properties.VariableNames),v=double(E.selected_edge_ARI);v=v(isfinite(v));if ~isempty(v),minari=min(v);end,end
for i=1:height(E)
    act=string(E.selected_action(i));fr=[double(E.r_from(i)) double(E.g_from(i))];to=[double(E.r_to(i)) double(E.g_to(i))];
    if act=="STOP",nstop=nstop+1;end
    if isequal(to,truth)&&act=="STOP",rejectInto=true;end
    if isequal(fr,truth)&&act=="ACCEPT",acceptBelow=true;end
    if act=="STOP" && ~isequal(fr,truth) && reachable_local(fr,truth),premature=true;end
end
end

function s=structure_failure_stage_local(d)
if d.structure_exact_POSTHOC,s="EXACT";return,end
if ~d.terminal_unique,s="MULTI_ANCHOR_OR_TERMINAL_AMBIGUITY";return,end
if ~d.alpha_emitted
    if d.method=="TRKM",s="NO_ACTION_ARCHIVED_MODAL_MISS";else,s="NO_ACTION_MPLUS_MISS";end
    return
end
if ~d.any_anchor_reaches_truth_POSTHOC,s="ANCHOR_UNREACHABLE_FROM_TRUTH_PATH";return,end
if d.rejected_edge_into_truth_POSTHOC,s="REJECTED_EDGE_INTO_TRUTH";return,end
if d.premature_stop_above_truth_POSTHOC,s="PREMATURE_STOP_ABOVE_TRUTH";return,end
if d.accepted_edge_below_truth_POSTHOC,s="ACCEPTED_EDGE_BELOW_TRUTH";return,end
if d.anchor_truth_in_set_POSTHOC,s="TRUTH_ANCHOR_PRESENT_BUT_TERMINAL_MISS";return,end
s="OTHER_REACHABLE_PATH_MISS";
end

function s=primary_failure_stage_local(d)
if d.joint_correct_POSTHOC,s="JOINT_CORRECT";return,end
if d.leverage_POSTHOC
    if ~d.action_correct_POSTHOC,s="ACTION_NONEMISSION";return,end
    if ~d.alpha_exact_POSTHOC,s="ALPHA_ERROR";return,end
    if ~d.structure_exact_POSTHOC,s="STRUCTURE_"+d.structure_failure_stage_POSTHOC;return,end
else
    if ~d.action_correct_POSTHOC,s="FALSE_POSITIVE_ACTION";return,end
    if ~d.structure_exact_POSTHOC,s="STRUCTURE_"+d.structure_failure_stage_POSTHOC;return,end
end
s="OTHER";
end

% =========================================================================
% SUMMARIES
% =========================================================================
function S=summarize_local(T,vars)
if isempty(T),S=table();return,end
G=findgroups(T(:,vars));ng=max(G);S=table();
for g=1:ng
    U=T(G==g,:);head=U(1,vars);
    vals={height(U),mean(double(U.action_correct_POSTHOC)),mean(double(U.structure_exact_POSTHOC)),mean(double(U.joint_correct_POSTHOC)), ...
        mean(double(U.ambiguous_output_POSTHOC)),safe_mean_local(U.truth_state_inlier_ARI_POSTHOC),mean(double(U.truth_state_ARI_ge_090_POSTHOC)), ...
        safe_mean_local(U.best_alpha0_inlier_ARI_POSTHOC),mean(double(U.best_alpha0_ARI_ge_090_POSTHOC)),safe_mean_local(U.truth_state_ARI_regret_POSTHOC), ...
        mean(double(U.anchor_truth_in_set_POSTHOC)),mean(double(U.any_anchor_reaches_truth_POSTHOC)),safe_mean_local(U.final_inlier_ARI_POSTHOC)};
    tail=cell2table(vals,'VariableNames',{'N','action_accuracy','structure_exact_rate','joint_correct_rate','ambiguous_rate', ...
        'mean_truth_state_ARI_POSTHOC','truth_state_ARI_ge090_rate_POSTHOC','mean_best_alpha0_ARI_POSTHOC','best_alpha0_ARI_ge090_rate_POSTHOC', ...
        'mean_truth_state_regret_POSTHOC','anchor_truth_in_set_rate_POSTHOC','any_anchor_reaches_truth_rate_POSTHOC','mean_final_ARI_materialized_POSTHOC'});
    S=[S;[head,tail]]; %#ok<AGROW>
end
end

function S=summarize_fit_local(T,vars)
if isempty(T),S=table();return,end
G=findgroups(T(:,vars));ng=max(G);S=table();
for g=1:ng
    U=T(G==g,:);head=U(1,vars);
    vals={height(U),safe_mean_local(U.truth_state_inlier_ARI_POSTHOC),q_local(U.truth_state_inlier_ARI_POSTHOC,.10),q_local(U.truth_state_inlier_ARI_POSTHOC,.50), ...
        mean(double(U.truth_state_ARI_ge_090_POSTHOC)),safe_mean_local(U.best_alpha0_inlier_ARI_POSTHOC),q_local(U.best_alpha0_inlier_ARI_POSTHOC,.50), ...
        mean(double(U.best_alpha0_ARI_ge_090_POSTHOC)),safe_mean_local(U.truth_state_ARI_regret_POSTHOC),mean(double(U.best_alpha0_structure_is_truth_POSTHOC))};
    tail=cell2table(vals,'VariableNames',{'N','mean_truth_state_ARI','q10_truth_state_ARI','median_truth_state_ARI','truth_state_ge090_rate', ...
        'mean_best_alpha0_ARI_POSTHOC','median_best_alpha0_ARI_POSTHOC','best_alpha0_ge090_rate_POSTHOC','mean_best_minus_truth_state_ARI_POSTHOC', ...
        'oracle_best_alpha0_structure_is_truth_rate_POSTHOC'});
    S=[S;[head,tail]]; %#ok<AGROW>
end
end

function S=count_strings_local(T,groupvars,stringvar)
if isempty(T),S=table();return,end
V=T(:,groupvars);V.category=string(T.(stringvar));
% Version-robust manual counting over unique group/category rows.
K=unique(V,'rows','stable');cnt=zeros(height(K),1);
for i=1:height(K)
    m=true(height(V),1);
    for j=1:numel(groupvars),m=m & equal_col_local(V.(groupvars{j}),K.(groupvars{j})(i));end
    m=m & string(V.category)==string(K.category(i));cnt(i)=sum(m);
end
S=K;S.N=cnt;
% Add within-group rate.
S.rate_within_group=nan(height(S),1);
for i=1:height(S)
    m=true(height(T),1);for j=1:numel(groupvars),m=m & equal_col_local(T.(groupvars{j}),S.(groupvars{j})(i));end
    S.rate_within_group(i)=S.N(i)/sum(m);
end
S=sortrows(S,[groupvars,{'category'}]);
end

function m=equal_col_local(v,x)
if isstring(v)||iscellstr(v)||iscategorical(v),m=string(v)==string(x);else,m=double(v)==double(x);end
end

function S=terminal_confusion_local(T)
V=T(:,{'method','r0_POSTHOC','g0_POSTHOC','terminal_unique','terminal_r','terminal_g'});K=unique(V,'rows','stable');cnt=zeros(height(K),1);
for i=1:height(K)
    m=string(V.method)==string(K.method(i)) & double(V.r0_POSTHOC)==double(K.r0_POSTHOC(i)) & double(V.g0_POSTHOC)==double(K.g0_POSTHOC(i)) & logical(V.terminal_unique)==logical(K.terminal_unique(i));
    if isfinite(double(K.terminal_r(i))),m=m & double(V.terminal_r)==double(K.terminal_r(i));else,m=m & ~isfinite(double(V.terminal_r));end
    if isfinite(double(K.terminal_g(i))),m=m & double(V.terminal_g)==double(K.terminal_g(i));else,m=m & ~isfinite(double(V.terminal_g));end
    cnt(i)=sum(m);
end
S=K;S.N=cnt;S=sortrows(S,{'method','r0_POSTHOC','g0_POSTHOC','N'},{'ascend','ascend','ascend','descend'});
end

function S=action_quantiles_local(T)
S=table();
% TFKM leverage/nonleverage distributions of frozen action evidence.
U=T(T.method=="TFKM",:);groups=unique(U(:,{'panel','contamination_POSTHOC'}),'rows','stable');
for i=1:height(groups)
    V=U(U.panel==groups.panel(i)&U.contamination_POSTHOC==groups.contamination_POSTHOC(i),:);
    row=table("TFKM",string(groups.panel(i)),string(groups.contamination_POSTHOC(i)),height(V), ...
        q_local(V.tfkm_Jmax_min_POSTHOC,.10),q_local(V.tfkm_Jmax_min_POSTHOC,.50),q_local(V.tfkm_Jmax_min_POSTHOC,.90), ...
        q_local(V.tfkm_Smed_max_POSTHOC,.10),q_local(V.tfkm_Smed_max_POSTHOC,.50),q_local(V.tfkm_Smed_max_POSTHOC,.90), ...
        q_local(V.tfkm_DS_min_POSTHOC,.10),q_local(V.tfkm_DS_min_POSTHOC,.50),q_local(V.tfkm_DS_min_POSTHOC,.90), ...
        NaN,NaN,NaN,NaN,NaN,NaN, ...
        'VariableNames',{'method','panel','contamination','N','J_q10','J_q50','J_q90','Smed_q10','Smed_q50','Smed_q90','DS_q10','DS_q50','DS_q90','dBIC_q10','dBIC_q50','dBIC_q90','sep_q10','sep_q50','sep_q90'});
    S=[S;row]; %#ok<AGROW>
end
U=T(T.method=="TRKM",:);groups=unique(U(:,{'panel','contamination_POSTHOC'}),'rows','stable');
for i=1:height(groups)
    V=U(U.panel==groups.panel(i)&U.contamination_POSTHOC==groups.contamination_POSTHOC(i),:);
    row=table("TRKM",string(groups.panel(i)),string(groups.contamination_POSTHOC(i)),height(V), ...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
        q_local(V.trkm_radial_delta_BIC,.10),q_local(V.trkm_radial_delta_BIC,.50),q_local(V.trkm_radial_delta_BIC,.90), ...
        q_local(V.trkm_radial_separation_sd,.10),q_local(V.trkm_radial_separation_sd,.50),q_local(V.trkm_radial_separation_sd,.90), ...
        'VariableNames',{'method','panel','contamination','N','J_q10','J_q50','J_q90','Smed_q10','Smed_q50','Smed_q90','DS_q10','DS_q50','DS_q90','dBIC_q10','dBIC_q50','dBIC_q90','sep_q10','sep_q50','sep_q90'});
    S=[S;row]; %#ok<AGROW>
end
end

% =========================================================================
% GENERIC HELPERS
% =========================================================================
function [anyR,allR]=reachability_set_local(M,truth)
if isempty(M),anyR=false;allR=false;return,end
v=false(size(M,1),1);for i=1:size(M,1),v(i)=reachable_local(M(i,:),truth);end;anyR=any(v);allR=all(v);
end
function tf=contains_row_local(M,row),tf=~isempty(M)&&any(all(M==row,2));end
function tf=reachable_local(start,target)
cur=double(start(:)');target=double(target(:)');tf=isequal(cur,target);guard=0;
while ~tf && has_pred_local(cur) && guard<20,cur=predecessor_local(cur);tf=isequal(cur,target);guard=guard+1;end
end
function tf=has_pred_local(rg),r=rg(1);g=rg(2);tf=(g>r+1)||(g==r+1&&r>1);end
function q=predecessor_local(rg),r=rg(1);g=rg(2);if g>r+1,q=[r,g-1];elseif g==r+1&&r>1,q=[r-1,g-1];else,error('No predecessor.');end,end
function M=componentwise_maximal_local(U)
if isempty(U),M=zeros(0,2);return,end;U=unique(double(U),'rows');keep=true(size(U,1),1);for i=1:size(U,1),for j=1:size(U,1),if i~=j&&all(U(j,:)>=U(i,:))&&any(U(j,:)>U(i,:)),keep(i)=false;break,end,end,end;M=sortrows(U(keep,:),[1 2]);
end
function m=safe_mean_local(x),x=double(x(:));x=x(isfinite(x));if isempty(x),m=NaN;else,m=mean(x);end,end
function m=safe_min_local(x),x=double(x(:));x=x(isfinite(x));if isempty(x),m=NaN;else,m=min(x);end,end
function m=safe_max_local(x),x=double(x(:));x=x(isfinite(x));if isempty(x),m=NaN;else,m=max(x);end,end
function q=q_local(x,p)
x=sort(double(x(isfinite(x))));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x;return,end;t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end
end
function f=f1_local(tp,pred,truth),den=pred+truth;if den<=0,f=NaN;else,f=2*tp/den;end,end
function [y,out]=truth_local(sim)
if isfield(sim,'labels'),y=double(sim.labels(:));elseif isfield(sim,'gtrue'),y=double(sim.gtrue(:));else,error('Truth labels missing.');end
if isfield(sim,'outlier'),out=logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),out=logical(sim.outlier_mask(:));else,error('Outlier mask missing.');end
assert(numel(y)==numel(out),'Truth length mismatch.');
end

function write_report_local(f,O,PF,SF,FD,TB,FP,MG,WG,SG,root)
fid=fopen(f,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'TFKM/TRKM FRESH POSTHOC DIAGNOSTIC V1\n');fprintf(fid,'========================================\n');
fprintf(fid,'TRUTH OPEN. NO REFITS. NO RETUNING. Frozen V1 is not modified.\n');
fprintf(fid,'POSTHOC oracle columns are diagnostic only.\nOutput root: %s\n\n',root);
fprintf(fid,'OVERALL FORENSIC\n%s\n',evalc('disp(O)'));
fprintf(fid,'PRIMARY FAILURE DECOMPOSITION\n%s\n',evalc('disp(PF)'));
fprintf(fid,'STRUCTURE FAILURE DECOMPOSITION\n%s\n',evalc('disp(SF)'));
fprintf(fid,'TRUE-STATE FIT DIAGNOSTIC\n%s\n',evalc('disp(FD)'));
fprintf(fid,'TFKM LEVERAGE ACTION BLOCKS\n%s\n',evalc('disp(TB)'));
fprintf(fid,'TRKM NONLEVERAGE FALSE-POSITIVE SOURCE\n%s\n',evalc('disp(FP)'));
fprintf(fid,'MAIN GEOMETRY\n%s\n',evalc('disp(MG)'));
fprintf(fid,'WEAK GEOMETRY\n%s\n',evalc('disp(WG)'));
fprintf(fid,'STRUCTURE PANEL\n%s\n',evalc('disp(SG)'));
end

function lock_truth_vendor_local(root)
tf=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));assert(isfolder(tf),'Canonical vendor missing.');addpath(genpath(char(tf)),'-begin');rehash path;
for q=["generate_fixed_dataset_v1","adjusted_rand_index_v1"],p=string(which(char(q)));assert(strlength(p)>0,'Missing %s.',q);assert(startsWith(p,tf,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,p);end
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";
for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
