function out=merge_final_fresh_forensics_v1(output_root)
%MERGE_FINAL_FRESH_FORENSICS_V1 Summarize why overselected/low-ARI fits fail or remain benign.
pkg=fileparts(mfilename('fullpath'));addpath(genpath(pkg),'-begin');output_root=char(string(output_root));
S=load(fullfile(output_root,'study.mat'));study=S.study;manifest=S.manifest;Sel=readtable(fullfile(output_root,'summary','final_selections.csv'),'TextType','string');
U=table();C=table();D=table();X=table();
for i=1:height(manifest)
    ss=char(manifest.session(i)); rep=double(manifest.replicate(i)); expected=sum(double(Sel.replicate)==rep & (logical(Sel.structural_over)|logical(Sel.low_ari_lt_095)));
    if expected==0, continue; end
    fu=fullfile(output_root,'forensics',sprintf('forensic_units_%s.csv',ss));fc=fullfile(output_root,'forensics',sprintf('forensic_clusters_%s.csv',ss));fd=fullfile(output_root,'forensics',sprintf('forensic_directions_%s.csv',ss));fx=fullfile(output_root,'forensics',sprintf('forensic_contingency_%s.csv',ss));
    assert(isfile(fu),'Missing forensic unit file %s',fu);T=readtable(fu,'TextType','string');assert(height(T)==expected,'Forensic %s incomplete: %d/%d.',ss,height(T),expected);U=[U;T]; %#ok<AGROW>
    if isfile(fc),T=readtable(fc,'TextType','string');C=[C;T];end %#ok<AGROW>
    if isfile(fd),T=readtable(fd,'TextType','string');D=[D;T];end %#ok<AGROW>
    if isfile(fx),T=readtable(fx,'TextType','string');X=[X;T];end %#ok<AGROW>
end
outdir=fullfile(output_root,'forensics');writetable(U,fullfile(outdir,'forensic_units_all.csv'));writetable(C,fullfile(outdir,'forensic_clusters_all.csv'));writetable(D,fullfile(outdir,'forensic_directions_all.csv'));writetable(X,fullfile(outdir,'forensic_contingency_all.csv'));
% Over-selection quality uses ALL selected units as denominator, not only forensic targets.
Over=table();for method=unique(string(Sel.method))'
    Z=Sel(string(Sel.method)==method,:);ov=logical(Z.structural_over);row=table(method,height(Z),mean(double(ov)),sum(ov), ...
        mean_or_nan(double(Z.ari_inlier(ov))),min_or_nan(double(Z.ari_inlier(ov))),mean_or_nan(double(Z.ari_inlier(ov))>=study.forensic_ari_threshold), ...
        mean(double(Z.low_ari_lt_095)),mean(double(Z.positive_regret),'omitnan'), ...
        'VariableNames',{'method','N','structural_over_rate','N_over','mean_ari_over','min_ari_over','frac_over_ari_ge_095','low_ari_rate','mean_positive_regret'});Over=[Over;row]; %#ok<AGROW>
end
writetable(Over,fullfile(outdir,'overselection_quality.csv'));
% Failure signatures among targeted diagnostic units.
F=count_strings(U,'method','failure_signature');writetable(F,fullfile(outdir,'failure_signature_counts.csv'));
if ~isempty(C),EC=C(logical(C.is_extra_cluster),:);ECcount=count_strings(EC,'method','cluster_signature');else,EC=table();ECcount=table();end
writetable(ECcount,fullfile(outdir,'extra_cluster_signature_counts.csv'));
if ~isempty(D),Dcount=count_strings(D,'method','dominant_external_association');else,Dcount=table();end
writetable(Dcount,fullfile(outdir,'extra_direction_association_counts.csv'));
if isempty(U),Low=table();else,Low=U(double(U.recomputed_ari)<study.forensic_ari_threshold,:);Low=sortrows(Low,'recomputed_ari','ascend');end;if ~isempty(Low),writetable(Low,fullfile(outdir,'low_ari_cases.csv'));end
fprintf('\n============================================================\nFRESH FORENSICS MERGED\n============================================================\n');
fprintf('\nOVERSELECTION QUALITY\n');disp(Over);
fprintf('\nLOW-ARI CASES (< %.2f)\n',study.forensic_ari_threshold);if isempty(Low),fprintf('None.\n');else,disp(Low(:,{'scenario_id','replicate','method','selected_r','selected_g','recomputed_ari','ari_after_truth_open_merge','subspace_coverage','failure_signature'}));end
fprintf('\nEXTRA CLUSTER TYPES\n');if isempty(ECcount),fprintf('None.\n');else,disp(ECcount);end
fprintf('\nEXTRA DIRECTION EXTERNAL ASSOCIATION\n');if isempty(Dcount),fprintf('None.\n');else,disp(Dcount);end
out=struct('units',U,'clusters',C,'directions',D,'contingency',X,'overselection_quality',Over,'failure_counts',F,'extra_cluster_counts',ECcount,'extra_direction_counts',Dcount,'low_ari_cases',Low);save(fullfile(outdir,'final_forensics_report.mat'),'out','-v7.3');
end
function T=count_strings(A,g1,g2)
if isempty(A),T=table();return;end;[G,K]=findgroups(A(:,{g1,g2}));T=K;T.count=splitapply(@numel,ones(height(A),1),G);end
function y=mean_or_nan(x),if isempty(x),y=NaN;else,y=mean(x,'omitnan');end,end
function y=min_or_nan(x),if isempty(x),y=NaN;else,y=min(x);end,end
