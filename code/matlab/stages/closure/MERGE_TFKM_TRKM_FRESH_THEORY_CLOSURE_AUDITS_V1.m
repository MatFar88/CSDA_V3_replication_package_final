function OUT = MERGE_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1(documents_root)
%MERGE_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1 Merge the three closure audits.
if nargin<1,documents_root="";end
[~,outroot]=load_study_local(documents_root);
outdir=fullfile(outroot,'fresh_theory_closure_audits_v1');assert(isfolder(outdir),'Closure-audit output directory missing.');

D=dir(fullfile(outdir,'worker_*_of_*_grank.csv'));assert(~isempty(D),'No g/r worker files found.');GR=table();
for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');GR=[GR;Q];end %#ok<AGROW>
GR=unique(GR,'rows','stable');keys=string(GR.scenario_id)+"|"+string(GR.replicate);[~,ia]=unique(keys,'stable');GR=GR(ia,:);GR=sortrows(GR,{'scenario_id','replicate'});
assert(height(GR)==1230,'Expected 1230 structural rows, found %d.',height(GR));

D=dir(fullfile(outdir,'worker_*_of_*_local_return.csv'));assert(~isempty(D),'No local-return worker files found.');LR=table();
for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');LR=[LR;Q];end %#ok<AGROW>
LR=unique(LR,'rows','stable');keys=string(LR.method)+"|"+string(LR.scenario_id)+"|"+string(LR.replicate)+"|"+string(LR.perturb_type)+"|"+string(LR.perturb_count);[~,ia]=unique(keys,'stable');LR=LR(ia,:);
LR=sortrows(LR,{'scenario_id','replicate','method','perturb_type','perturb_count'});
assert(height(LR)==24600,'Expected 24600 local-return rows, found %d.',height(LR));

writetable(GR,fullfile(outdir,'CLOSURE_GRANK_ALL.csv'));save(fullfile(outdir,'CLOSURE_GRANK_ALL.mat'),'GR','-v7.3');
writetable(LR,fullfile(outdir,'CLOSURE_LOCAL_RETURN_ALL.csv'));save(fullfile(outdir,'CLOSURE_LOCAL_RETURN_ALL.mat'),'LR','-v7.3');

% -------------------------------------------------------------------------
% Baseline contract with the already merged four-audit package.
% -------------------------------------------------------------------------
base=LR(string(LR.perturb_type)=="BASELINE" & double(LR.perturb_count)==0,:);
assert(height(base)==2460,'Expected 2460 baseline method-units, found %d.',height(base));
oldfile=fullfile(outroot,'fresh_theory_audits_v1','THEORY_AUDITS_ALL.csv');
assert(isfile(oldfile),'Previous THEORY_AUDITS_ALL.csv is required for the baseline contract.');
OLD=readtable(oldfile,'TextType','string');
K=OLD(:,{'method','scenario_id','replicate','assignment_inliers_exact','one_step_trim_exact','one_step_truth_state_recovered','sample_generalized_gap'});
B=innerjoin(base,K,'Keys',{'method','scenario_id','replicate'});assert(height(B)==2460,'Baseline contract join mismatch.');
B.assignment_match=logical(B.return_assignment_exact)==logical(B.assignment_inliers_exact);
B.trim_match=logical(B.return_trim_exact)==logical(B.one_step_trim_exact);
B.full_match=logical(B.return_full_state_exact)==logical(B.one_step_truth_state_recovered);
B.gap_match=abs(double(B.sample_generalized_gap_after_perturb)-double(B.sample_generalized_gap)) <= 1e-10*(1+abs(double(B.sample_generalized_gap)));
assert(all(B.assignment_match) && all(B.trim_match) && all(B.full_match) && all(B.gap_match), ...
    'BASELINE CONTRACT FAIL: closure local map does not reproduce Audit 1 exactly.');
writetable(B,fullfile(outdir,'AUDIT7_BASELINE_CONTRACT.csv'));

% Audit A summary: empirical g path.
A5=group_summary_grank_local(GR);
writetable(A5,fullfile(outdir,'AUDIT5_GPATH_SUMMARY.csv'));

% Audit B summary: literal rank replay.
A6=group_summary_rank_local(GR);
writetable(A6,fullfile(outdir,'AUDIT6_RANK_LITERAL_SUMMARY.csv'));

% Audit C summary: perturbation-return curve (applicable rows only).
A7=group_summary_local_return_local(LR);
writetable(A7,fullfile(outdir,'AUDIT7_LOCAL_RETURN_SUMMARY.csv'));

write_report_local(fullfile(outdir,'THREE_CLOSURE_AUDITS_REPORT.txt'),A5,A6,A7,B);

fprintf('\n============================================================\n');
fprintf('THREE THEORY CLOSURE AUDITS MERGED\n');
fprintf('1230/1230 structural units | 24600/24600 local-return rows\n');
fprintf('Baseline contract with Audit 1: PASS 2460/2460 (assignment, trim, full state, eigengap)\n');
fprintf('============================================================\n');
fprintf('\nAUDIT 5 - EMPIRICAL g-PATH\n');disp(A5);
fprintf('\nAUDIT 6 - LITERAL CROSS-FITTED RANK REPLAY\n');disp(A6);
fprintf('\nAUDIT 7 - LOCAL PERTURBATION/RETURN\n');disp(A7);
fprintf('Output: %s\n',outdir);
OUT=struct('grank',GR,'local_return',LR,'gpath_summary',A5,'rank_summary',A6,'local_summary',A7,'baseline_contract',B,'output_root',string(outdir));
end

function S=group_summary_grank_local(T)
[G,panel,cont,g0]=findgroups(string(T.panel),string(T.contamination),double(T.g0));
N=splitapply(@numel,T.scenario_id,G);
sep=splitapply(@mean,double(T.g_change_point_separated),G);
medG=splitapply(@(x)median(double(x),'omitnan'),T.Gg_emp,G);
q10G=splitapply(@(x)qlinear_local(x,.10),T.Gg_emp,G);
medPre=splitapply(@(x)median(double(x),'omitnan'),T.min_pre_gain,G);
medPost=splitapply(@(x)median(double(x),'omitnan'),T.first_overfit_gain,G);
S=table(panel,cont,g0,N,sep,medG,q10G,medPre,medPost,'VariableNames', ...
    {'panel','contamination','g0','N','rate_empirical_change_point_separated','median_Gg_emp','q10_Gg_emp','median_min_pre_gain','median_first_overfit_gain'});
end

function S=group_summary_rank_local(T)
[G,r0,g0]=findgroups(double(T.r0),double(T.g0));N=splitapply(@numel,T.scenario_id,G);
weyl=splitapply(@mean,double(T.weyl_sufficient_pass),G);
exact=splitapply(@mean,double(T.rank_exact_midpoint),G);
exactPass=splitapply(@conditional_mean_local,double(T.rank_exact_midpoint),double(T.weyl_sufficient_pass),G);
exactFail=splitapply(@conditional_mean_not_local,double(T.rank_exact_midpoint),double(T.weyl_sufficient_pass),G);
medG=splitapply(@(x)median(double(x),'omitnan'),T.G_r_cross,G);
medR=splitapply(@(x)median(double(x),'omitnan'),T.rhat_midpoint,G);
S=table(r0,g0,N,medG,weyl,exact,exactPass,exactFail,medR,'VariableNames', ...
    {'r0','g0','N','median_G_r_cross','rate_weyl_sufficient_pass','rate_rank_exact_midpoint','rank_exact_given_weyl_pass','rank_exact_given_weyl_fail','median_rhat_midpoint'});
end

function S=group_summary_local_return_local(T)
A=T(logical(T.applicable),:);
[G,method,panel,cont,ptype,k]=findgroups(string(A.method),string(A.panel),string(A.contamination),string(A.perturb_type),double(A.perturb_count));
N=splitapply(@numel,A.scenario_id,G);
assign=splitapply(@mean,double(A.return_assignment_exact),G);trim=splitapply(@mean,double(A.return_trim_exact),G);full=splitapply(@mean,double(A.return_full_state_exact),G);
ari=splitapply(@(x)mean(double(x),'omitnan'),A.return_inlier_ARI,G);dist=splitapply(@(x)median(double(x),'omitnan'),A.designed_state_distance,G);
gap=splitapply(@(x)median(double(x),'omitnan'),A.sample_generalized_gap_after_perturb,G);
S=table(method,panel,cont,ptype,k,N,dist,assign,trim,full,ari,gap,'VariableNames', ...
    {'method','panel','contamination','perturb_type','perturb_count','N','median_designed_state_distance','rate_assignment_return_exact','rate_trim_return_exact','rate_full_state_return_exact','mean_return_inlier_ARI','median_postperturb_sample_gap'});
end

function v=conditional_mean_local(a,b),q=logical(b);if ~any(q),v=NaN;else,v=mean(double(a(q)));end,end
function v=conditional_mean_not_local(a,b),q=~logical(b);if ~any(q),v=NaN;else,v=mean(double(a(q)));end,end
function q=qlinear_local(x,p),x=sort(double(x(:)));x=x(isfinite(x));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x(1);return,end;t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end,end

function write_report_local(file,A5,A6,A7,B)
fid=fopen(file,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'TFKM/TRKM FRESH THEORY CLOSURE AUDITS V1\n');fprintf(fid,'==========================================\n\n');
fprintf(fid,'Truth-open post-fit diagnostics on the frozen 1230-DGP study. No selector is retuned.\n');
fprintf(fid,'Baseline local-map contract: %d/%d exact matches with previous Audit 1.\n\n',sum(B.full_match),height(B));
fprintf(fid,'AUDIT 5 - EMPIRICAL COMPONENT-RESPECTING g PATH\n%s\n',evalc('disp(A5)'));
fprintf(fid,'Interpretation: Gg_emp>1 means the realized component-respecting marginal gains before g0 dominate the first empirical overfit split gain. It is a diagnostic of the theorem mechanism, not a new truth-blind selector.\n\n');
fprintf(fid,'AUDIT 6 - LITERAL CROSS-FITTED RANK REPLAY\n%s\n',evalc('disp(A6)'));
fprintf(fid,'Interpretation: rhat_midpoint counts eigenvalues of Bhat_x above beta_min/2. This truth-open midpoint directly replays the proposition but is not an operational selector because beta_min is unknown in practice.\n\n');
fprintf(fid,'AUDIT 7 - LOCAL BOUNDARY-DIRECTED PERTURBATION/RETURN\n%s\n',evalc('disp(A7)'));
fprintf(fid,'Interpretation: the curve is a deterministic local stress test. It does not estimate a formal attraction radius or prove uniform return over a metric ball.\n');
end

function [root,outroot]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found. Pass Documents explicitly.');
end
