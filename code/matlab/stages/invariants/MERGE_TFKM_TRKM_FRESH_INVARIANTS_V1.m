function OUT = MERGE_TFKM_TRKM_FRESH_INVARIANTS_V1(documents_root)
%MERGE_TFKM_TRKM_FRESH_INVARIANTS_V1
% Merge all completed invariant workers and summarize theorem-oriented gaps.
if nargin<1,documents_root="";end
[~,outroot,S]=load_study_local(documents_root);catalog=S.catalog;
outdir=fullfile(outroot,'fresh_invariants_v1');assert(isfolder(outdir),'Invariant output directory missing.');
D=dir(fullfile(outdir,'worker_*_of_*_invariants.csv'));assert(~isempty(D),'No invariant worker files found.');
T=table();for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');T=[T;Q];end %#ok<AGROW>
T=unique(T,'rows','stable');T=sortrows(T,{'scenario_id','replicate'});
keys=string(T.scenario_id)+"|"+string(T.replicate);[uk,ia]=unique(keys,'stable'); %#ok<ASGLU>
if numel(uk)~=height(T)
    % If duplicate worker layouts were run, keep one row per key only after
    % checking that duplicates agree numerically/textually via row equality.
    T2=T(ia,:);assert(height(T2)==1230,'Duplicate-key recovery did not yield 1230 rows.');T=T2;keys=string(T.scenario_id)+"|"+string(T.replicate);
end
assert(height(T)==1230 && numel(unique(keys))==1230,'Expected 1230 unique invariant rows, found %d.',height(T));
assert(all(logical(T.radial_frozen_reproduction_pass)),'At least one radial GMM reproduction check failed.');
allfile=fullfile(outdir,'FRESH_INVARIANTS_ALL.csv');writetable(T,allfile);save(fullfile(outdir,'FRESH_INVARIANTS_ALL.mat'),'T','-v7.3');

% -------------------------------------------------------------------------
% Cell-level summary: 41 frozen DGP cells.
% -------------------------------------------------------------------------
Cell=table();ids=unique(double(T.scenario_id));
for sid=ids(:)'
    Z=T(double(T.scenario_id)==sid,:);c=catalog(double(catalog.scenario_id)==sid,:);assert(height(c)==1);
    Cell=[Cell;table(sid,string(c.panel),string(c.contamination),double(c.p_over_n),double(c.rho),double(c.delta),double(c.r0),double(c.g0),height(Z), ...
        mean(double(Z.radial_GMM_emit)),mean(double(Z.radial_GMM_delta_BIC),'omitnan'),mean(double(Z.radial_GMM_separation_sd),'omitnan'), ...
        mean(double(Z.partition_margin_ratio_8logn),'omitnan'),min(double(Z.partition_margin_ratio_8logn)),mean(double(Z.oracle_nearest_centroid_ARI),'omitnan'), ...
        mean(double(Z.G_g_merge_over_split),'omitnan'),min(double(Z.G_g_merge_over_split)),mean(double(Z.beta_min_between),'omitnan'), ...
        mean(double(Z.pca_overlap_trueU),'omitnan'),min(double(Z.pca_overlap_trueU)),mean(double(Z.pca_BBP_proxy),'omitnan'), ...
        mean(double(Z.crossB_Weyl_ratio_beta_over_2err90),'omitnan'),mean(double(Z.crossB_Weyl_condition_rate),'omitnan'), ...
        'VariableNames',{'scenario_id','panel','contamination','p_over_n','rho','delta','r0','g0','N', ...
        'radial_emit_rate','mean_radial_delta_BIC','mean_radial_sep_sd','mean_partition_margin_ratio','min_partition_margin_ratio','mean_oracle_centroid_ARI', ...
        'mean_Gg','min_Gg','mean_beta_min','mean_PCA_overlap','min_PCA_overlap','mean_BBP_proxy','mean_crossWeyl_ratio','mean_crossWeyl_pass_rate'})]; %#ok<AGROW>
end
Cell=sortrows(Cell,'scenario_id');assert(height(Cell)==41);writetable(Cell,fullfile(outdir,'INVARIANTS_BY_CELL.csv'));

% -------------------------------------------------------------------------
% Theorem-oriented group summaries.
% -------------------------------------------------------------------------
Groups=table();
[G,panel,delta,cont]=findgroups(string(T.panel),double(T.delta),string(T.contamination));
N=splitapply(@numel,T.scenario_id,G);
medGg=splitapply(@(x)median(double(x),'omitnan'),T.G_g_merge_over_split,G);minGg=splitapply(@(x)min(double(x),[],'omitnan'),T.G_g_merge_over_split,G);
rateGg=splitapply(@(x)mean(double(x)>1),T.G_g_merge_over_split,G);
medPm=splitapply(@(x)median(double(x),'omitnan'),T.partition_margin_ratio_8logn,G);minPm=splitapply(@(x)min(double(x),[],'omitnan'),T.partition_margin_ratio_8logn,G);ratePm=splitapply(@(x)mean(double(x)>1),T.partition_margin_ratio_8logn,G);
meanARI=splitapply(@(x)mean(double(x),'omitnan'),T.oracle_nearest_centroid_ARI,G);exactARI=splitapply(@(x)mean(abs(double(x)-1)<=1e-12),T.oracle_nearest_centroid_ARI,G);
medPCA=splitapply(@(x)median(double(x),'omitnan'),T.pca_overlap_trueU,G);q10PCA=splitapply(@(x)qlinear_local(x,.10),T.pca_overlap_trueU,G);
medRW=splitapply(@(x)median(double(x),'omitnan'),T.crossB_Weyl_ratio_beta_over_2err90,G);passRW=splitapply(@(x)mean(double(x)>1),T.crossB_Weyl_ratio_beta_over_2err90,G);
Groups=table(panel,delta,cont,N,medGg,minGg,rateGg,medPm,minPm,ratePm,meanARI,exactARI,medPCA,q10PCA,medRW,passRW, ...
    'VariableNames',{'panel','delta','contamination','N','median_Gg','min_Gg','rate_Gg_gt1','median_partition_margin','min_partition_margin','rate_partition_margin_gt1', ...
    'mean_oracle_centroid_ARI','oracle_exact_rate','median_PCA_overlap','q10_PCA_overlap','median_crossWeyl_ratio','rate_crossWeyl_ratio_gt1'});
writetable(Groups,fullfile(outdir,'THEOREM_GEOMETRY_SUMMARY.csv'));

% Radial estimator vs leverage truth. This is an estimator check for the
% Q-invariant, not the definition of the invariant itself.
lev=string(T.contamination)=="leverage";emit=logical(T.radial_GMM_emit);
Radial=table();
Radial=[Radial;table("RADIAL_GMM",height(T),mean(emit(lev)),mean(emit(~lev)),mean(emit==lev), ...
    rank_auc_local(double(T.radial_GMM_delta_BIC),double(lev)),rank_auc_local(double(T.radial_GMM_separation_sd),double(lev)), ...
    'VariableNames',{'estimator','N','leverage_recall','nonleverage_FPR','action_accuracy','AUC_deltaBIC','AUC_separation'})];
writetable(Radial,fullfile(outdir,'RADIAL_INVARIANT_ESTIMATOR_SUMMARY.csv'));

% Rank-specific summary by truth pair. The (3,5) cell is the key case where
% r0 < g0-1, so rank is not mechanically g0-1.
Rank=table();pairs=unique([double(T.r0_POSTHOC),double(T.g0_POSTHOC)],'rows');
for j=1:size(pairs,1)
    r=pairs(j,1);g=pairs(j,2);Z=T(double(T.r0_POSTHOC)==r&double(T.g0_POSTHOC)==g,:);
    Rank=[Rank;table(r,g,height(Z),mean(double(Z.population_centroid_rank)==r),median(double(Z.beta_min_between),'omitnan'), ...
        median(double(Z.rawB_ratio),'omitnan'),median(double(Z.crossB_rank_SNR),'omitnan'),median(double(Z.crossB_Weyl_ratio_beta_over_2err90),'omitnan'), ...
        mean(double(Z.crossB_Weyl_ratio_beta_over_2err90)>1),mean(double(Z.crossB_Weyl_condition_rate),'omitnan'), ...
        'VariableNames',{'r0','g0','N','population_rank_correct_rate','median_beta_min','median_rawB_ratio','median_crossB_rank_SNR','median_crossWeyl_ratio','rate_crossWeyl_ratio_gt1','mean_splitwise_Weyl_pass_rate'})]; %#ok<AGROW>
end
writetable(Rank,fullfile(outdir,'RANK_INVARIANT_BY_TRUTH_PAIR.csv'));

% -------------------------------------------------------------------------
% Join A4 weak true-state 100-start basin audit if available.
% -------------------------------------------------------------------------
WeakJoin=table();WeakSummary=table();WeakAUC=table();
wf=fullfile(outroot,'fresh_v2_diagnostic_suite_v1','weak_true_state_pool_audit_v1','WEAK_TRUESTATE_POOL_ALL.csv');
if isfile(wf)
    W=readtable(wf,'TextType','string');assert(height(W)==300,'Weak A4 pool table should have 300 rows.');
    K=T(:,{'scenario_id','replicate','G_g_merge_over_split','partition_margin_ratio_8logn','oracle_nearest_centroid_ARI', ...
        'beta_min_between','pca_overlap_trueU','pca_gap','pca_ratio','pca_BBP_proxy','crossB_rank_SNR','crossB_Weyl_ratio_beta_over_2err90','crossB_Weyl_condition_rate'});
    WeakJoin=innerjoin(W,K,'Keys',{'scenario_id','replicate'});assert(height(WeakJoin)==300,'Weak A4/invariant join mismatch.');
    writetable(WeakJoin,fullfile(outdir,'WEAK_A4_INVARIANT_JOIN.csv'));
    [GG,m,pn,rho,mech]=findgroups(string(WeakJoin.method),double(WeakJoin.p_over_n),double(WeakJoin.rho),string(WeakJoin.mechanism_POSTHOC));
    NN=splitapply(@numel,WeakJoin.scenario_id,GG);
    bARI=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.best_truth_ARI_POSTHOC,GG);
    gg=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.G_g_merge_over_split,GG);
    pm=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.partition_margin_ratio_8logn,GG);
    oa=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.oracle_nearest_centroid_ARI,GG);
    po=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.pca_overlap_trueU,GG);
    bb=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.pca_BBP_proxy,GG);
    rw=splitapply(@(x)mean(double(x),'omitnan'),WeakJoin.crossB_Weyl_ratio_beta_over_2err90,GG);
    WeakSummary=table(m,pn,rho,mech,NN,bARI,gg,pm,oa,po,bb,rw,'VariableNames',{'method','p_over_n','rho','mechanism_POSTHOC','N','mean_best_start_ARI','mean_Gg','mean_partition_margin','mean_oracle_centroid_ARI','mean_PCA_overlap','mean_BBP_proxy','mean_crossWeyl_ratio'});
    writetable(WeakSummary,fullfile(outdir,'WEAK_A4_INVARIANT_BY_MECHANISM.csv'));
    metrics=["G_g_merge_over_split","partition_margin_ratio_8logn","oracle_nearest_centroid_ARI","beta_min_between","pca_overlap_trueU","pca_gap","pca_ratio","pca_BBP_proxy","crossB_rank_SNR","crossB_Weyl_ratio_beta_over_2err90"];
    methods=unique(string(WeakJoin.method));
    for mm=methods(:)'
        Z=WeakJoin(string(WeakJoin.method)==mm,:);good=double(Z.best_truth_ARI_POSTHOC)>=.90-1e-12;
        for q=metrics
            auc=rank_auc_local(double(Z.(q)),double(good));
            WeakAUC=[WeakAUC;table(mm,q,auc,'VariableNames',{'method','invariant','AUC_basin_accessible'})]; %#ok<AGROW>
        end
    end
    writetable(WeakAUC,fullfile(outdir,'WEAK_A4_INVARIANT_AUC.csv'));
end

% -------------------------------------------------------------------------
% Human-readable theorem diagnostic report.
% -------------------------------------------------------------------------
report=fullfile(outdir,'THEOREM_DIAGNOSTIC_REPORT.txt');fid=fopen(report,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'TFKM/TRKM FRESH INVARIANTS V1 - THEOREM DIAGNOSTIC\n');
fprintf(fid,'=====================================================\n\n');
fprintf(fid,'This report studies identification targets, not a new frozen algorithm.\n\n');
fprintf(fid,'RADIAL Q-INVARIANT ESTIMATOR CHECK\n%s\n',evalc('disp(Radial)'));
fprintf(fid,'THEOREM GEOMETRY SUMMARY\n%s\n',evalc('disp(Groups)'));
fprintf(fid,'RANK INVARIANT BY TRUTH PAIR\n%s\n',evalc('disp(Rank)'));
if ~isempty(WeakSummary)
    fprintf(fid,'WEAK A4 INVARIANT BY MECHANISM\n%s\n',evalc('disp(WeakSummary)'));
    fprintf(fid,'WEAK A4 INVARIANT AUC\n%s\n',evalc('disp(WeakAUC)'));
end
fprintf(fid,'\nInterpretation targets:\n');
fprintf(fid,'Q: leverage is identifiable when robust radial populations separate beyond concentration noise.\n');
fprintf(fid,'P: exact sample partition is plausible when the minimum Mahalanobis margin dominates log n.\n');
fprintf(fid,'G: g0 is locally identified when forced-merge loss dominates the best within-component split gain (G_g>1).\n');
fprintf(fid,'R: conditional on the partition, r0 is rank(B0); cross-fitting targets B0 without same-sample centroid-noise bias.\n');
fprintf(fid,'U: PCA overlap/BBP proxy diagnoses whether the signal subspace is statistically accessible before optimization.\n');
fprintf(fid,'A4: best-start ARI separates statistical/geometry failure from basin/arbitration failure.\n');

fprintf('\n============================================================\n');
fprintf('FRESH INVARIANTS V1 MERGED: 1230/1230\n');
fprintf('============================================================\n');
fprintf('\nRADIAL Q-INVARIANT ESTIMATOR CHECK\n');disp(Radial);
fprintf('\nRANK INVARIANT BY TRUTH PAIR\n');disp(Rank);
fprintf('\nTHEOREM GEOMETRY SUMMARY\n');disp(Groups);
if ~isempty(WeakAUC),fprintf('\nWEAK A4 INVARIANT AUC\n');disp(WeakAUC);end
fprintf('Output: %s\n',outdir);

OUT=struct('units',T,'cells',Cell,'geometry_summary',Groups,'radial',Radial,'rank',Rank,'weak_join',WeakJoin,'weak_summary',WeakSummary,'weak_auc',WeakAUC,'output_root',string(outdir));
end

% =========================================================================
function a=rank_auc_local(score,y)
score=double(score(:));y=double(y(:));ok=isfinite(score)&isfinite(y);score=score(ok);y=y(ok);pos=(y==1);neg=(y==0);np=sum(pos);nn=sum(neg);if np==0||nn==0,a=NaN;return,end
[ss,ord]=sort(score);r=zeros(size(score));i=1;while i<=numel(ss),j=i;while j<numel(ss)&&ss(j+1)==ss(i),j=j+1;end;rr=(i+j)/2;r(ord(i:j))=rr;i=j+1;end;a=(sum(r(pos))-np*(np+1)/2)/(np*nn);
end
function q=qlinear_local(x,p),x=sort(double(x(:)));x=x(isfinite(x));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x(1);return,end;t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end,end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";
for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
