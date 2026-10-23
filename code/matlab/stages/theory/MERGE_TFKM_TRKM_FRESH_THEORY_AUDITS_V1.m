function OUT = MERGE_TFKM_TRKM_FRESH_THEORY_AUDITS_V1(documents_root)
%MERGE_TFKM_TRKM_FRESH_THEORY_AUDITS_V1 Merge the four direct theory audits.
if nargin<1,documents_root="";end
[~,outroot,S]=load_study_local(documents_root); %#ok<ASGLU>
outdir=fullfile(outroot,'fresh_theory_audits_v1');assert(isfolder(outdir),'Audit output directory missing.');
D=dir(fullfile(outdir,'worker_*_of_*_theory_audits.csv'));assert(~isempty(D),'No worker files found.');
T=table();for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');T=[T;Q];end %#ok<AGROW>
T=unique(T,'rows','stable');T=sortrows(T,{'scenario_id','replicate','method'});
keys=string(T.method)+"|"+string(T.scenario_id)+"|"+string(T.replicate);[~,ia]=unique(keys,'stable');T=T(ia,:);keys=string(T.method)+"|"+string(T.scenario_id)+"|"+string(T.replicate);
assert(height(T)==2460 && numel(unique(keys))==2460,'Expected 2460 unique method-unit rows, found %d.',height(T));
writetable(T,fullfile(outdir,'THEORY_AUDITS_ALL.csv'));save(fullfile(outdir,'THEORY_AUDITS_ALL.mat'),'T','-v7.3');

% Audit 1: strict true-state prerequisites and one-step return.
A1=group_summary_local(T,{'method','panel','contamination'}, ...
    {'strict_prerequisite_quotient','one_step_truth_state_recovered','assignment_inliers_exact','one_step_trim_exact'}, ...
    {'sample_generalized_gap','assignment_margin_inliers','trim_residual_margin_est'});
writetable(A1,fullfile(outdir,'AUDIT1_STRICTNESS_SUMMARY.csv'));

% Audit 2: method-specific spectral accessibility.
A2=group_summary_local(T,{'method','panel','p_over_n','rho','delta'}, ...
    {'population_target_correct','tfkm_population_condition_pass'}, ...
    {'population_generalized_gap','operator_error','operator_error_over_gap','G_spec_gap_over_error','sample_target_cstar','sample_target_overlap_mean'});
writetable(A2,fullfile(outdir,'AUDIT2_SPECTRAL_SUMMARY.csv'));

% Audit 3: projected effective margin.
A3=group_summary_local(T,{'method','panel','p_over_n','rho','delta'}, ...
    {'effective_projected_margin_gt1','projected_oracle_exact'}, ...
    {'intrinsic_partition_margin_8logn','sample_target_cstar','effective_projected_margin','projected_oracle_ARI'});
writetable(A3,fullfile(outdir,'AUDIT3_PROJECTED_MARGIN_SUMMARY.csv'));

% Audit 4: residual trimming order. Keep all contamination types as controls,
% but the theorem-facing positive target is diffuse leverage.
A4=group_summary_local(T,{'method','panel','contamination','p_over_n','rho','delta'}, ...
    {'oracle_trim_mask_exact','estimated_trim_mask_exact'}, ...
    {'oracle_trim_gap','estimated_trim_gap','max_residual_perturbation','G_trim_gap_over_2error'});
writetable(A4,fullfile(outdir,'AUDIT4_TRIMMING_GAP_SUMMARY.csv'));

% WEAK geometry transition: direct theorem quantities + old A4 basin audit.
W=T(string(T.panel)=="WEAK_DELTA7",:);
Weak=group_summary_local(W,{'method','p_over_n','rho'}, ...
    {'strict_prerequisite_quotient','one_step_truth_state_recovered','projected_oracle_exact','effective_projected_margin_gt1'}, ...
    {'G_spec_gap_over_error','sample_target_cstar','sample_target_overlap_mean','effective_projected_margin'});
wf=fullfile(outroot,'fresh_v2_diagnostic_suite_v1','weak_true_state_pool_audit_v1','WEAK_TRUESTATE_POOL_ALL.csv');
WeakJoin=table();WeakAUC=table();
if isfile(wf)
    P=readtable(wf,'TextType','string');assert(height(P)==300,'A4 pool table should have 300 rows.');
    K=T(:,{'method','scenario_id','replicate','G_spec_gap_over_error','sample_target_cstar','sample_target_overlap_mean','effective_projected_margin','projected_oracle_exact','strict_prerequisite_quotient','one_step_truth_state_recovered'});
    WeakJoin=innerjoin(P,K,'Keys',{'method','scenario_id','replicate'});assert(height(WeakJoin)==300,'A4 direct-theory join mismatch.');
    WeakJoin.basin_accessible_POSTHOC=double(WeakJoin.best_truth_ARI_POSTHOC)>=.90-1e-12;
    WeakJoin.selected_good_POSTHOC=double(WeakJoin.selected_truth_ARI_POSTHOC)>=.90-1e-12;
    writetable(WeakJoin,fullfile(outdir,'WEAK_A4_DIRECT_THEORY_JOIN.csv'));
    metrics=["G_spec_gap_over_error","sample_target_cstar","sample_target_overlap_mean","effective_projected_margin"];
    methods=unique(string(WeakJoin.method));
    for mm=methods(:)'
        Z=WeakJoin(string(WeakJoin.method)==mm,:);y=double(Z.basin_accessible_POSTHOC);
        for q=metrics
            auc=rank_auc_local(double(Z.(q)),y);
            WeakAUC=[WeakAUC;table(mm,q,auc,'VariableNames',{'method','metric','AUC_basin_accessible'})]; %#ok<AGROW>
        end
    end
    writetable(WeakAUC,fullfile(outdir,'WEAK_A4_DIRECT_THEORY_AUC.csv'));
    % Append basin coverage/arbitration by geometry to a separate table.
    [G,m,pn,rho]=findgroups(string(WeakJoin.method),double(WeakJoin.p_over_n),double(WeakJoin.rho));
    N=splitapply(@numel,WeakJoin.scenario_id,G);
    basin=splitapply(@mean,double(WeakJoin.basin_accessible_POSTHOC),G);
    sel=splitapply(@mean,double(WeakJoin.selected_good_POSTHOC),G);
    arb=splitapply(@(a,b)conditional_mean_local(a,b),double(WeakJoin.selected_good_POSTHOC),double(WeakJoin.basin_accessible_POSTHOC),G);
    gs=splitapply(@(x)median(double(x),'omitnan'),WeakJoin.G_spec_gap_over_error,G);
    cs=splitapply(@(x)median(double(x),'omitnan'),WeakJoin.sample_target_cstar,G);
    ge=splitapply(@(x)median(double(x),'omitnan'),WeakJoin.effective_projected_margin,G);
    WeakTransition=table(m,pn,rho,N,gs,cs,ge,basin,sel,arb,'VariableNames',{'method','p_over_n','rho','N','median_Gspec','median_cstar','median_GPeff','basin_accessible_rate','selected_good_rate','arbitration_success_given_basin'});
    writetable(WeakTransition,fullfile(outdir,'WEAK_THEORY_TRANSITION.csv'));
else
    WeakTransition=table();
end
writetable(Weak,fullfile(outdir,'WEAK_DIRECT_THEORY_SUMMARY.csv'));

% Join frozen true-state candidate diagnostics if present.
pf=fullfile(outroot,'fresh_posthoc_diagnostic_v1','POSTHOC_UNIT_DIAGNOSTIC.csv');
FrozenJoin=table();
if isfile(pf)
    P=readtable(pf,'TextType','string');K=T(:,{'method','scenario_id','replicate','G_spec_gap_over_error','sample_target_cstar','effective_projected_margin','one_step_truth_state_recovered','G_trim_gap_over_2error'});
    FrozenJoin=innerjoin(P(:,{'method','scenario_id','replicate','truth_state_inlier_ARI_POSTHOC','truth_state_trim_F1_POSTHOC'}),K,'Keys',{'method','scenario_id','replicate'});
    writetable(FrozenJoin,fullfile(outdir,'FROZEN_TRUESTATE_THEORY_JOIN.csv'));
end

write_report_local(fullfile(outdir,'FOUR_THEORY_AUDITS_REPORT.txt'),A1,A2,A3,A4,Weak,WeakTransition,WeakAUC);

fprintf('\n============================================================\n');
fprintf('FOUR THEORY AUDITS MERGED: 2460/2460 method-units\n');
fprintf('============================================================\n');
fprintf('\nAUDIT 1 - STRICTNESS / ONE-STEP\n');disp(A1);
fprintf('\nAUDIT 2 - SPECTRAL ACCESSIBILITY\n');disp(A2);
fprintf('\nAUDIT 3 - PROJECTED EFFECTIVE MARGIN\n');disp(A3);
fprintf('\nAUDIT 4 - TRIMMING GAP\n');disp(A4);
if ~isempty(WeakTransition),fprintf('\nWEAK THEORY TRANSITION + A4\n');disp(WeakTransition);end
if ~isempty(WeakAUC),fprintf('\nWEAK DIRECT-THEORY AUC FOR BASIN ACCESSIBILITY\n');disp(WeakAUC);end
fprintf('Output: %s\n',outdir);
OUT=struct('units',T,'audit1',A1,'audit2',A2,'audit3',A3,'audit4',A4,'weak',Weak,'weak_transition',WeakTransition,'weak_auc',WeakAUC,'frozen_join',FrozenJoin,'output_root',string(outdir));
end

function S=group_summary_local(T,groupVars,rateVars,metricVars)
if isempty(T),S=table();return,end
[Gu,keyTable]=findgroups(T(:,groupVars));
S=keyTable;N=splitapply(@numel,T.scenario_id,Gu);S.N=N;
for q=string(rateVars)
    S.("rate_"+q)=splitapply(@(x)mean(double(x),'omitnan'),T.(q),Gu);
end
for q=string(metricVars)
    S.("median_"+q)=splitapply(@(x)median(double(x),'omitnan'),T.(q),Gu);
    S.("q10_"+q)=splitapply(@(x)qlinear_local(x,.10),T.(q),Gu);
    S.("q90_"+q)=splitapply(@(x)qlinear_local(x,.90),T.(q),Gu);
end
end

function v=conditional_mean_local(sel,basin)
q=logical(basin);if ~any(q),v=NaN;else,v=mean(double(sel(q)));end
end

function a=rank_auc_local(score,y)
score=double(score(:));y=double(y(:));ok=isfinite(score)&isfinite(y);score=score(ok);y=y(ok);pos=(y==1);neg=(y==0);np=sum(pos);nn=sum(neg);if np==0||nn==0,a=NaN;return,end
[ss,ord]=sort(score);r=zeros(size(score));i=1;while i<=numel(ss),j=i;while j<numel(ss)&&ss(j+1)==ss(i),j=j+1;end;rr=(i+j)/2;r(ord(i:j))=rr;i=j+1;end;a=(sum(r(pos))-np*(np+1)/2)/(np*nn);
end
function q=qlinear_local(x,p),x=sort(double(x(:)));x=x(isfinite(x));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x(1);return,end;t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end,end

function write_report_local(file,A1,A2,A3,A4,Weak,WeakTransition,WeakAUC)
fid=fopen(file,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'TFKM/TRKM FRESH THEORY AUDITS V1\n');fprintf(fid,'=================================\n\n');
fprintf(fid,'These are truth-open audits on the existing 1230 DGPs. They do not define new tuned selectors.\n\n');
fprintf(fid,'AUDIT 1 - STRICT TRUE-STATE PREREQUISITE / ONE-STEP SELF-MAP\n%s\n',evalc('disp(A1)'));
fprintf(fid,'Interpretation: positive sample eigengap + positive inlier assignment margin + positive trim-order margin are finite-sample prerequisites for the local convergence theorem. The direct endpoint is whether one update from the truth-retained state recovers the truth-retained state. Labels on trimmed points are deliberately not part of this quotient-state diagnostic.\n\n');
fprintf(fid,'AUDIT 2 - METHOD-SPECIFIC SPECTRAL ACCESSIBILITY\n%s\n',evalc('disp(A2)'));
fprintf(fid,'Interpretation: Gspec = population generalized eigengap / operator error. The theorem needs error/gap -> 0, not merely a detectable spike. TFKM also reports its population ordering condition.\n\n');
fprintf(fid,'AUDIT 3 - EFFECTIVE PROJECTED EXACT-RECOVERY MARGIN\n%s\n',evalc('disp(A3)'));
fprintf(fid,'Interpretation: GPeff combines intrinsic centroid separation with the weakest retained signal fraction cstar and orthogonal variance rho. The direct response is exact nearest-population-centroid classification after projection onto the true-state REG-HD subspace.\n\n');
fprintf(fid,'AUDIT 4 - RESIDUAL TRIMMING SEPARATION\n%s\n',evalc('disp(A4)'));
fprintf(fid,'Interpretation: oracle gap b = min(outlier residual)-max(inlier residual); Gtrim=b/(2 max residual perturbation) is a finite sufficient order-preservation diagnostic. Diffuse leverage is the theorem-facing target; bridge is a negative/control geometry, clean is vacuous.\n\n');
fprintf(fid,'WEAK DIRECT SUMMARY\n%s\n',evalc('disp(Weak)'));
if ~isempty(WeakTransition),fprintf(fid,'WEAK THEORY TRANSITION + A4\n%s\n',evalc('disp(WeakTransition)'));end
if ~isempty(WeakAUC),fprintf(fid,'WEAK DIRECT-THEORY AUC\n%s\n',evalc('disp(WeakAUC)'));end
end

function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found. Pass Documents explicitly.');S=load(fullfile(outroot,'study.mat'));
end
