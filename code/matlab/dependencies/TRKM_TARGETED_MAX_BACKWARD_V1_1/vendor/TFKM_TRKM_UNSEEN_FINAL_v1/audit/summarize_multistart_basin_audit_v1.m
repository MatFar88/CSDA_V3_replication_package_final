function report=summarize_multistart_basin_audit_v1(audit_root)
%SUMMARIZE_MULTISTART_BASIN_AUDIT_V1 Merge pool-level audit outputs.
pool_root=fullfile(audit_root,'pools'); summary_root=fullfile(audit_root,'summary');
if ~isfolder(summary_root),mkdir(summary_root);end
R=read_many(pool_root,'*__rules.csv'); S=read_many(pool_root,'*__starts.csv'); F=read_many(pool_root,'*__families.csv');
assert(~isempty(R),'No completed audit pools found.');
writetable(R,fullfile(summary_root,'all_pool_rules.csv'));
writetable(S,fullfile(summary_root,'all_pool_starts.csv'));
writetable(F,fullfile(summary_root,'all_family_summary.csv'));

% Aggregate rule performance, keeping truth-blind and oracle rows distinct.
keys=unique(R(:,{'method','target_type','rule'}),'rows','stable'); A=table();
for i=1:height(keys)
    q=string(R.method)==string(keys.method(i)) & string(R.target_type)==string(keys.target_type(i)) & string(R.rule)==string(keys.rule(i));
    X=R(q,:);
    row=table(string(keys.method(i)),string(keys.target_type(i)),string(keys.rule(i)),height(X), ...
        mean(double(X.selected_ari),'omitnan'),min(double(X.selected_ari),[],'omitnan'), ...
        mean(double(X.selector_regret),'omitnan'),max(double(X.selector_regret),[],'omitnan'), ...
        mean(double(X.selected_ari)>=.90,'omitnan'),mean(double(X.selected_ari)>=.95,'omitnan'), ...
        'VariableNames',{'method','target_type','rule','n_pools','mean_selected_ari','min_selected_ari', ...
        'mean_oracle_regret','max_oracle_regret','rate_ari_ge_090','rate_ari_ge_095'});
    A=[A;row]; %#ok<AGROW>
end
writetable(A,fullfile(summary_root,'rule_summary.csv'));

% Pair each TFKM challenger against frozen current on the same fixed triplet.
T=R(string(R.method)=="REG_TFKM_HD" & string(R.target_type)=="CURRENT_SELECTED",:);
cur=T(string(T.rule)=="CURRENT_SIGN_ALL_STABLE",:);
challenger_names=["RV_CV_ALL","FB_STABLE_ALL","HIER_STABLE4","HIER_RVCV4","ALL_STABLE_NAIVE_DIAGNOSTIC","MINOBJ_ALL_DIAGNOSTIC"];
D=table();
for j=1:numel(challenger_names)
    ch=T(string(T.rule)==challenger_names(j),:);
    if isempty(ch),continue;end
    [cur2,ch2]=align_pool_rows(cur,ch);
    delta=double(ch2.selected_ari)-double(cur2.selected_ari);
    for cont=[23001 23004]
        q=double(ch2.scenario_id)==cont;
        if ~any(q),continue;end
        cname="clean";if cont==23004,cname="bridge";end
        row=table(challenger_names(j),cname,sum(q),mean(delta(q),'omitnan'),min(delta(q),[],'omitnan'),max(delta(q),[],'omitnan'), ...
            mean(delta(q)>.01),mean(delta(q)<-.01), ...
            'VariableNames',{'rule','contamination','n_pools','mean_delta_ari_vs_current','min_delta_ari_vs_current', ...
            'max_delta_ari_vs_current','rate_improve_gt_001','rate_degrade_gt_001'});
        D=[D;row]; %#ok<AGROW>
    end
end
writetable(D,fullfile(summary_root,'tfkm_challenger_vs_current.csv'));

% Pool-level availability diagnostics from frozen-current rows.
P=R(ismember(string(R.rule),["CURRENT_SIGN_ALL_STABLE","CURRENT_RV_CV"]),:);
P.material_basin_gap=double(P.basin_gap_all_vs_sign)>=.05;
P.material_selector_regret=double(P.selector_regret)>=.05;
writetable(P,fullfile(summary_root,'current_pool_diagnostics.csv'));

% Text report.
report_file=fullfile(summary_root,'BASIN_AUDIT_REPORT.txt'); fid=fopen(report_file,'w');
cl=onCleanup(@() fclose(fid));
fprintf(fid,'TFKM/TRKM MULTISTART BASIN AUDIT v1\n');
fprintf(fid,'=====================================\n\n');
fprintf(fid,'Truth is used only for diagnostics. Challenger scores are truth-blind.\n');
fprintf(fid,'A 0.05 gap is reported only as a forensic flag, never as an estimator threshold.\n\n');

Tc=P(string(P.method)=="REG_TFKM_HD" & string(P.target_type)=="CURRENT_SELECTED",:);
if ~isempty(Tc)
    fprintf(fid,'TFKM CURRENT-SELECTED CLEAN/BRIDGE POOLS\n');
    fprintf(fid,'n pools: %d\n',height(Tc));
    fprintf(fid,'mean current ARI: %.6f\n',mean(double(Tc.selected_ari),'omitnan'));
    fprintf(fid,'mean oracle-best-all ARI: %.6f\n',mean(double(Tc.best_all_ari),'omitnan'));
    fprintf(fid,'mean current selector regret: %.6f\n',mean(double(Tc.selector_regret),'omitnan'));
    fprintf(fid,'max current selector regret: %.6f\n',max(double(Tc.selector_regret),[],'omitnan'));
    fprintf(fid,'oracle best is non-SIGN: %d / %d\n',sum(string(Tc.oracle_best_family)~="SIGN_PCA_TKM"),height(Tc));
    fprintf(fid,'material ALL-vs-SIGN availability gap >= .05: %d / %d\n',sum(double(Tc.basin_gap_all_vs_sign)>=.05),height(Tc));
    fprintf(fid,'material current selector regret >= .05: %d / %d\n\n',sum(double(Tc.selector_regret)>=.05),height(Tc));
end

Tr=P(string(P.method)=="REG_TRKM_HD" & string(P.target_type)=="CURRENT_SELECTED",:);
if ~isempty(Tr)
    fprintf(fid,'TRKM CURRENT RV_CV CLEAN/BRIDGE POOLS\n');
    fprintf(fid,'n pools: %d\n',height(Tr));
    fprintf(fid,'mean current ARI: %.6f\n',mean(double(Tr.selected_ari),'omitnan'));
    fprintf(fid,'mean oracle-best-all ARI: %.6f\n',mean(double(Tr.best_all_ari),'omitnan'));
    fprintf(fid,'mean RV_CV selector regret: %.6f\n',mean(double(Tr.selector_regret),'omitnan'));
    fprintf(fid,'max RV_CV selector regret: %.6f\n',max(double(Tr.selector_regret),[],'omitnan'));
    fprintf(fid,'material RV_CV selector regret >= .05: %d / %d\n\n',sum(double(Tr.selector_regret)>=.05),height(Tr));
end

if ~isempty(D)
    fprintf(fid,'TFKM TRUTH-BLIND CHALLENGERS AT THE SAME FIXED (r,g,alpha)\n');
    for i=1:height(D)
        fprintf(fid,'%-28s %-6s n=%2d mean dARI=%+.5f min=%+.5f max=%+.5f improve>0.01=%.3f degrade>0.01=%.3f\n', ...
            char(D.rule(i)),char(D.contamination(i)),D.n_pools(i),D.mean_delta_ari_vs_current(i), ...
            D.min_delta_ari_vs_current(i),D.max_delta_ari_vs_current(i),D.rate_improve_gt_001(i),D.rate_degrade_gt_001(i));
    end
end
fprintf(fid,'\nInterpretation rule:\n');
fprintf(fid,'1) If best_all >> best_SIGN, the frozen SIGN-only eligibility discarded an available better basin.\n');
fprintf(fid,'2) If best_all is also poor, the fixed triplet itself is structurally incapable on that dataset; basin arbitration is not the main cause.\n');
fprintf(fid,'3) Challenger improvements are fixed-triplet evidence only. They do NOT retroactively revalidate the full r,g,alpha selector.\n');

report=struct('rules',R,'starts',S,'families',F,'rule_summary',A,'challenger_delta',D,'current_diagnostics',P,'report_file',string(report_file));
end

function T=read_many(root,pattern)
d=dir(fullfile(root,pattern));T=table();
for i=1:numel(d)
    q=readtable(fullfile(d(i).folder,d(i).name),'TextType','string');
    if isempty(T),T=q;else,T=[T;q];end %#ok<AGROW>
end
end
function [A,B]=align_pool_rows(A,B)
keyA=make_key(A);keyB=make_key(B);[common,ia,ib]=intersect(keyA,keyB,'stable'); %#ok<ASGLU>
A=A(ia,:);B=B(ib,:);
end
function k=make_key(T)
k=string(double(T.scenario_id))+"|"+string(double(T.replicate))+"|"+string(T.method)+"|"+string(T.target_type)+"|"+ ...
    string(double(T.r))+"|"+string(double(T.g))+"|"+compose('%.6f',double(T.alpha));
end
