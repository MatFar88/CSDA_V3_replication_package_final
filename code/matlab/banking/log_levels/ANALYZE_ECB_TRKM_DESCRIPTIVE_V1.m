function OUT = ANALYZE_ECB_TRKM_DESCRIPTIVE_V1(top_k)
%ANALYZE_ECB_TRKM_DESCRIPTIVE_V1 Read-only analysis of the ECB TRKM descriptive fit.
%
% This function is explicitly POST-SELECTION and DESCRIPTIVE ONLY.
% It performs NO refit, NO reselection, and reads NO bank names.
%
% It is intended for the ECB application outcome in which the primary TRKM
% selector returns a unique structural state but NO positive trimming action,
% so the pre-specified secondary alpha=0 refit has been run with
% RUN_ECB_TRKM_DESCRIPTIVE_REFIT_V1.
%
% Usage, from ECB_TRKM_APPLICATION_V1:
%   D = ANALYZE_ECB_TRKM_DESCRIPTIVE_V1;
%   D = ANALYZE_ECB_TRKM_DESCRIPTIVE_V1(25);
%
% Anonymous outputs are written to output/analysis_descriptive/.

if nargin<1 || isempty(top_k), top_k=25; end
top_k=double(top_k);
assert(isscalar(top_k) && isfinite(top_k) && top_k>=1 && top_k==floor(top_k), ...
    'top_k must be a positive integer.');

root=fileparts(mfilename('fullpath'));
CFG=ecb_trkm_config_v1(root);
ecb_trkm_path_lock_v1(root);

primary_file=fullfile(CFG.final_root,'ECB_TRKM_PRIMARY_FINAL.mat');
desc_file=fullfile(CFG.descriptive_root,'ECB_TRKM_DESCRIPTIVE_REFIT.mat');
assert(isfile(primary_file),'Primary final output not found. Run FINALIZE first.');
assert(isfile(desc_file),'Descriptive refit not found. Run RUN_ECB_TRKM_DESCRIPTIVE_REFIT_V1 first.');
assert(isfile(CFG.input_file),'Internal ECB363 input not found.');

P=load(primary_file,'Prediction','Audit');
D=load(desc_file,'Descriptive','Labels','Trim','U','fit');
S=load(CFG.input_file,'X','Lab_vars','original_row_index','protocol_version','pair_seed');
Pred=P.Prediction;

assert(~logical(Pred.final_fit_available), ...
    'Primary TRKM already defines a terminal fit; use a primary-fit diagnostic instead.');
assert(~logical(Pred.action_ambiguous),'Primary TRKM output is action-ambiguous.');
assert(logical(Pred.terminal_unique),'A unique retained structural pair is required.');
assert(~logical(Pred.alpha_emitted),'This diagnostic is written for the no-action ECB TRKM outcome.');
assert(string(S.protocol_version)==CFG.version,'Input protocol version drift.');
assert(double(S.pair_seed)==CFG.pair_seed,'Pair-seed drift.');

X=double(S.X);
orig=double(S.original_row_index(:));
labels=double(D.Labels(:));
trim=logical(D.Trim(:));
U=double(D.U);
varcode=labels_to_strings_local(S.Lab_vars,CFG.p);

assert(isequal(size(X),[CFG.n CFG.p]),'Input dimension drift.');
assert(numel(labels)==CFG.n && numel(trim)==CFG.n,'Descriptive fit dimension mismatch.');
r=double(Pred.terminal_r); g=double(Pred.terminal_g);
assert(size(U,1)==CFG.p && size(U,2)==r,'Descriptive subspace dimension mismatch.');
assert(all(~trim),'No-action descriptive alpha=0 refit unexpectedly trims observations.');
assert(abs(double(D.Descriptive.alpha))<=CFG.tol,'Descriptive alpha is not zero.');
assert(double(D.Descriptive.r)==r && double(D.Descriptive.g)==g,'Descriptive structural state drift.');

ret=~trim;
clusters=sort(unique(labels(ret)));
clusters=clusters(isfinite(clusters)&clusters>0);
assert(numel(clusters)==g,'Descriptive fit has %d positive clusters; expected g=%d.',numel(clusters),g);

% -------------------------------------------------------------------------
% 1. Latent geometry and deterministic orientation.
% -------------------------------------------------------------------------
mu=mean(X(ret,:),1);
XR=X(ret,:)-mu;
Uor=U;
sign_vec=ones(1,r);
for h=1:r
    [~,imax]=max(abs(Uor(:,h)));
    if Uor(imax,h)<0
        Uor(:,h)=-Uor(:,h);
        sign_vec(h)=-1;
    end
end
Z=(X-mu)*Uor;
Centroids=zeros(g,r);
ClusterSizes=zeros(g,1);
for k=1:g
    ii=ret & labels==clusters(k);
    ClusterSizes(k)=sum(ii);
    Centroids(k,:)=mean(Z(ii,:),1);
end
ClusterTable=table(clusters,ClusterSizes,ClusterSizes/sum(ret), ...
    'VariableNames',{'cluster','retained_n','retained_fraction'});
for h=1:r
    ClusterTable.(sprintf('latent_centroid_%d',h))=Centroids(:,h);
end

% -------------------------------------------------------------------------
% 2. Retained-data row-space support of the selected subspace.
% -------------------------------------------------------------------------
[~,Ss,V]=svd(XR,'econ');
sval=diag(Ss);
if isempty(sval)
    tol_rank=0; rankX=0;
else
    tol_rank=max(size(XR))*eps(max(sval));
    rankX=sum(sval>tol_rank);
end
Vr=V(:,1:rankX);
row_mass=nan(r,1); null_mass=nan(r,1);
Urow=zeros(size(Uor));
for h=1:r
    uh=Uor(:,h); uh=uh/norm(uh);
    ur=Vr*(Vr'*uh);
    row_mass(h)=min(max(norm(ur)^2,0),1);
    null_mass(h)=max(0,1-row_mass(h));
    if norm(ur)>0, Urow(:,h)=ur/norm(ur); else, Urow(:,h)=nan(size(ur)); end
end
SubspaceSupport=table((1:r)',row_mass,null_mass, ...
    'VariableNames',{'component','rowspace_mass','nullspace_mass'});

% -------------------------------------------------------------------------
% 3. Top raw/data-supported loadings and score correlations.
% -------------------------------------------------------------------------
K=min(top_k,CFG.p);
RawLoadings=table();
DataSupportedLoadings=table();
for h=1:r
    [~,ord]=sort(abs(Uor(:,h)),'descend');ord=ord(1:K);
    T=table(repmat(h,K,1),(1:K)',ord(:),varcode(ord),Uor(ord,h),abs(Uor(ord,h)), ...
        'VariableNames',{'component','rank_abs','variable_index','indicator_code','loading','abs_loading'});
    RawLoadings=[RawLoadings;T]; %#ok<AGROW>
    if all(isfinite(Urow(:,h)))
        [~,od]=sort(abs(Urow(:,h)),'descend');od=od(1:K);
        T=table(repmat(h,K,1),(1:K)',od(:),varcode(od),Urow(od,h),abs(Urow(od,h)), ...
            'VariableNames',{'component','rank_abs','variable_index','indicator_code','loading','abs_loading'});
        DataSupportedLoadings=[DataSupportedLoadings;T]; %#ok<AGROW>
    end
end

% For r=1, realized score correlations are particularly interpretable.
ScoreCorrelations=table();
size_corr=NaN; size_R2=NaN;
if r==1
    score=Z(:,1);
    rho=nan(CFG.p,1);
    for j=1:CFG.p, rho(j)=pair_corr_local(X(:,j),score); end
    [~,ord]=sort(abs(rho),'descend','MissingPlacement','last');
    ord=ord(1:min(K,sum(isfinite(rho))));
    ScoreCorrelations=table((1:numel(ord))',ord(:),varcode(ord),rho(ord),abs(rho(ord)), ...
        'VariableNames',{'rank_abs_corr','variable_index','indicator_code','corr_with_latent_score','abs_corr'});
    size_corr=pair_corr_local(X(:,1),score); size_R2=size_corr^2;
end

% -------------------------------------------------------------------------
% 4. Total-assets association and cluster contrasts.
% -------------------------------------------------------------------------
SizeByCluster=table();
for k=1:g
    ii=labels==clusters(k);
    if r==1, ms=mean(Z(ii,1)); else, ms=NaN; end
    T=table(clusters(k),sum(ii),mean(X(ii,1)),median(X(ii,1)),ms, ...
        'VariableNames',{'cluster','n','mean_log_total_assets','median_log_total_assets','mean_latent_score_1'});
    SizeByCluster=[SizeByCluster;T]; %#ok<AGROW>
end

ContrastDecomposition=table();
if g==2
    i1=labels==clusters(1); i2=labels==clusters(2);
    sd=std(X,0,1);
    d=mean(X(i2,:),1)-mean(X(i1,:),1);
    ds=d./sd; ds(sd<=eps|~isfinite(sd))=NaN;
    [~,ord]=sort(abs(ds),'descend','MissingPlacement','last');
    ord=ord(1:min(K,sum(isfinite(ds))));
    nz1=zeros(numel(ord),1); nz2=nz1; med1=nan(numel(ord),1); med2=med1;
    for h=1:numel(ord)
        j=ord(h);x1=X(i1,j);x2=X(i2,j);
        nz1(h)=mean(x1~=0);nz2(h)=mean(x2~=0);
        if any(x1~=0),med1(h)=median(x1(x1~=0));end
        if any(x2~=0),med2(h)=median(x2(x2~=0));end
    end
    ContrastDecomposition=table((1:numel(ord))',ord(:),varcode(ord),ds(ord)', ...
        nz1,nz2,nz2-nz1,med1,med2, ...
        'VariableNames',{'rank_abs_std_diff','variable_index','indicator_code','standardized_difference_2_minus_1', ...
        'nonzero_rate_cluster_1','nonzero_rate_cluster_2','nonzero_rate_difference_2_minus_1', ...
        'median_nonzero_cluster_1','median_nonzero_cluster_2'});
end

% -------------------------------------------------------------------------
% 5. Projected residual distribution. With alpha=0 this is descriptive only.
% -------------------------------------------------------------------------
R2all=inf(CFG.n,1);
for k=1:g
    Dk=Z-repmat(Centroids(k,:),CFG.n,1);
    R2all=min(R2all,sum(Dk.^2,2));
end
q=[0 .25 .50 .75 .90 .95 1];qr=quantile(R2all,q);
ResidualSummary=array2table([CFG.n mean(R2all) median(R2all) qr(4) qr(5) qr(6) qr(7)], ...
    'VariableNames',{'n','mean_sq_projected_residual','median_sq_projected_residual','q75','q90','q95','max'});

% -------------------------------------------------------------------------
% 6. Frozen TRKM selector paths and candidate diagnostics.
% -------------------------------------------------------------------------
ArchivedPath=readtable(fullfile(CFG.final_root,'ECB_TRKM_ARCHIVED_PER_ALPHA.csv'),'TextType','string');
raw_file=fullfile(CFG.final_root,'ECB_TRKM_RAW_PER_ALPHA.csv');
if isfile(raw_file)
    RawPath=readtable(raw_file,'TextType','string');
else
    RawPath=table();
end
C=load_candidate_rows_local(CFG);
SelectedState=C(double(C.r)==r & double(C.g)==g,:);
SelectedState=sortrows(SelectedState,'alpha');
SelectedStatePath=minimal_candidate_rows_local(SelectedState);
SelectedStatePath.is_descriptive_alpha0=abs(double(SelectedStatePath.alpha))<=CFG.tol;

% Deterministic descriptive refit should reproduce the stored (r,g,alpha=0)
% production candidate. Check labels, trimming and projector, but do not alter
% or select anything based on this check.
row0=C(double(C.r)==r & double(C.g)==g & abs(double(C.alpha))<=CFG.tol,:);
assert(height(row0)==1,'Stored descriptive-state candidate is missing or duplicated.');
k0=double(row0.candidate_index);
Q0=load(fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k0)),'Labels','Trim','U');
keep=~trim & ~logical(Q0.Trim(:));
ari_reproduction=adjusted_rand_index_v1(labels(keep),double(Q0.Labels(keep)));
trim_reproduction=isequal(trim,logical(Q0.Trim(:)));
if size(Q0.U,2)==r
    BU=orth(double(U)); B0=orth(double(Q0.U));
    projector_diff=norm(BU*BU'-B0*B0',2);
else
    projector_diff=NaN;
end
ReproductionCheck=table(k0,ari_reproduction,logical(trim_reproduction),projector_diff, ...
    'VariableNames',{'candidate_index_alpha0','ARI_descriptive_vs_stored','trim_mask_identical','projector_operator_difference'});

% Radial scores were part of the PRIMARY no-action decision; summarize them
% by the descriptive groups without turning them into a new trimming rule.
Rad=readtable(fullfile(CFG.final_root,'ECB_TRKM_RADIAL_SCORES_ANONYMOUS.csv'),'TextType','string');
assert(height(Rad)==CFG.n,'Radial score file has wrong number of rows.');
radq=double(Rad.radial_q);
RadialByCluster=table();
for k=1:g
    ii=labels==clusters(k);x=radq(ii);
    T=table(clusters(k),sum(ii),mean(x),median(x),quantile(x,.90),quantile(x,.95),max(x), ...
        'VariableNames',{'cluster','n','mean_radial_q','median_radial_q','q90_radial_q','q95_radial_q','max_radial_q'});
    RadialByCluster=[RadialByCluster;T]; %#ok<AGROW>
end
PrimaryAction=table(string(Pred.action_status),logical(Pred.alpha_emitted),double(Pred.archived_delta_RV), ...
    logical(Pred.radial_emit),double(Pred.radial_high_weight),double(Pred.radial_delta_BIC),double(Pred.radial_separation_sd), ...
    'VariableNames',{'action_status','alpha_emitted','archived_delta_RV','radial_emit','radial_high_weight','radial_delta_BIC','radial_separation_sd'});

% -------------------------------------------------------------------------
% 7. Optional, explicitly secondary TFKM-vs-TRKM partition comparison.
%    This is run only if the finalized companion TFKM application is found.
% -------------------------------------------------------------------------
CrossMethodSummary=table();CrossTab=table();TrimByTRKM=table();
tfkm_root=fullfile(fileparts(root),'ECB_TFKM_APPLICATION_V1');
tfkm_final=fullfile(tfkm_root,'output','final','ECB_TFKM_PRIMARY_FINAL.mat');
tfkm_input=fullfile(tfkm_root,'output','ECB363_INPUT_INTERNAL.mat');
if isfile(tfkm_final) && isfile(tfkm_input)
    F=load(tfkm_final,'Prediction','FinalLabels','FinalTrim','FinalU');
    I=load(tfkm_input,'original_row_index');
    if logical(F.Prediction.final_fit_available) && ~logical(F.Prediction.final_fit_ambiguous) && ...
            isequal(double(I.original_row_index(:)),orig)
        flab=double(F.FinalLabels(:));ftrim=logical(F.FinalTrim(:));fU=double(F.FinalU);
        common=~ftrim; % TRKM descriptive alpha=0 retains all rows.
        ari_common=adjusted_rand_index_v1(flab(common),labels(common));
        if ~isempty(fU) && size(fU,1)==CFG.p && size(fU,2)==1 && r==1
            cosang=abs((fU(:,1)'*U(:,1))/(norm(fU(:,1))*norm(U(:,1))));
            cosang=min(max(cosang,0),1);angle_deg=acosd(cosang);overlap_sq=cosang^2;
        else
            angle_deg=NaN;overlap_sq=NaN;
        end
        CrossMethodSummary=table(sum(common),sum(ftrim),ari_common,overlap_sq,angle_deg, ...
            'VariableNames',{'TFKM_retained_n','TFKM_trimmed_n','ARI_on_TFKM_retained','rank1_subspace_overlap_squared','rank1_principal_angle_degrees'});

        fcl=sort(unique(flab(common)));tcl=sort(unique(labels(common)));
        CT=zeros(numel(fcl),numel(tcl));
        for i=1:numel(fcl),for j=1:numel(tcl),CT(i,j)=sum(common & flab==fcl(i) & labels==tcl(j));end,end
        [ii,jj]=ndgrid(1:numel(fcl),1:numel(tcl));
        CrossTab=table(fcl(ii(:)),tcl(jj(:)),CT(:), ...
            'VariableNames',{'TFKM_cluster','TRKM_cluster','n_on_TFKM_retained'});

        for j=1:numel(clusters)
            in=(labels==clusters(j));
            T=table(clusters(j),sum(in),sum(in & ftrim),mean(ftrim(in)), ...
                'VariableNames',{'TRKM_cluster','n_total','n_TFKM_trimmed','TFKM_trim_rate'});
            TrimByTRKM=[TrimByTRKM;T]; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
% 8. Compact descriptive summary and exports.
% -------------------------------------------------------------------------
Summary=table("TRKM","ECB_Q4_2014_363x361",CFG.n,CFG.p,r,g,0,CFG.n,0,rankX,CFG.p-rankX, ...
    mean(row_mass),mean(null_mass),double(D.Descriptive.RV_CV),double(D.Descriptive.W_ratio),double(D.Descriptive.R2), ...
    size_corr,size_R2,string(Pred.action_status), ...
    'VariableNames',{'method','dataset','n','p','r_hat','g_hat','descriptive_alpha','retained_n','trimmed_n', ...
    'rank_centered_X','ambient_null_dimension','mean_rowspace_mass','mean_nullspace_mass','RV_CV','W_ratio','R2', ...
    'corr_latent_score_log_total_assets','R2_latent_score_log_total_assets','primary_action_status'});

outroot=fullfile(CFG.output_root,'analysis_descriptive');
if ~isfolder(outroot),mkdir(outroot);end
writetable(Summary,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_SUMMARY.csv'));
writetable(ClusterTable,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_CLUSTER_SIZES.csv'));
writetable(SubspaceSupport,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_ROWSPACE_SUPPORT.csv'));
writetable(RawLoadings,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_RAW_LOADINGS.csv'));
if ~isempty(DataSupportedLoadings),writetable(DataSupportedLoadings,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_DATA_SUPPORTED_LOADINGS.csv'));end
if ~isempty(ScoreCorrelations),writetable(ScoreCorrelations,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_SCORE_CORRELATIONS.csv'));end
writetable(SizeByCluster,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_SIZE_BY_CLUSTER.csv'));
if ~isempty(ContrastDecomposition),writetable(ContrastDecomposition,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_CONTRAST_ZERO_DECOMPOSITION.csv'));end
writetable(ResidualSummary,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_RESIDUAL_SUMMARY.csv'));
writetable(ArchivedPath,fullfile(outroot,'ECB_TRKM_ARCHIVED_PER_ALPHA_COPY.csv'));
if ~isempty(RawPath),writetable(RawPath,fullfile(outroot,'ECB_TRKM_RAW_PER_ALPHA_COPY.csv'));end
writetable(SelectedStatePath,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_STATE_ALPHA_PATH.csv'));
writetable(ReproductionCheck,fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_REPRODUCTION_CHECK.csv'));
writetable(PrimaryAction,fullfile(outroot,'ECB_TRKM_PRIMARY_ACTION_DIAGNOSTIC.csv'));
writetable(RadialByCluster,fullfile(outroot,'ECB_TRKM_RADIAL_BY_DESCRIPTIVE_CLUSTER.csv'));
if ~isempty(CrossMethodSummary)
    writetable(CrossMethodSummary,fullfile(outroot,'ECB_TFKM_TRKM_SECONDARY_COMPARISON.csv'));
    writetable(CrossTab,fullfile(outroot,'ECB_TFKM_TRKM_CROSSTAB_ON_TFKM_RETAINED.csv'));
    writetable(TrimByTRKM,fullfile(outroot,'ECB_TFKM_TRIM_RATE_BY_TRKM_DESCRIPTIVE_CLUSTER.csv'));
end
save(fullfile(outroot,'ECB_TRKM_DESCRIPTIVE_ANALYSIS.mat'), ...
    'Summary','ClusterTable','SubspaceSupport','RawLoadings','DataSupportedLoadings','ScoreCorrelations', ...
    'SizeByCluster','ContrastDecomposition','ResidualSummary','ArchivedPath','RawPath','SelectedStatePath', ...
    'ReproductionCheck','PrimaryAction','RadialByCluster','CrossMethodSummary','CrossTab','TrimByTRKM', ...
    'Centroids','mu','sign_vec','R2all','Uor','Urow','sval','-v7.3');

fprintf('\n============================================================\n');
fprintf('ECB TRKM SECONDARY DESCRIPTIVE ANALYSIS\n');
fprintf('NO REFIT / NO RESELECTION / NO BANK NAMES\n');
fprintf('============================================================\n');
fprintf('Primary selector: %s | alpha emitted=%d\n',string(Pred.action_status),logical(Pred.alpha_emitted));
fprintf('Secondary descriptive state: r=%d | g=%d | alpha=0.000\n',r,g);
fprintf('\nSummary:\n');disp(Summary);
fprintf('\nDescriptive cluster sizes and oriented latent centroids:\n');disp(ClusterTable);
fprintf('\nRow-space support of descriptive subspace:\n');disp(SubspaceSupport);
fprintf('\nArchived per-alpha structural winners used by primary no-action structure:\n');disp(ArchivedPath);
fprintf('\nRAW per-alpha structural winners (diagnostic only under no-action):\n');disp(RawPath);
fprintf('\n(r,g)=(%d,%d) alpha path:\n',r,g);disp(SelectedStatePath);
fprintf('\nDescriptive refit versus stored alpha=0 candidate:\n');disp(ReproductionCheck);
fprintf('\nPrimary radial/action diagnostic:\n');disp(PrimaryAction);
fprintf('\nRadial score distribution by descriptive cluster:\n');disp(RadialByCluster);
fprintf('\nTotal-assets association by descriptive cluster:\n');disp(SizeByCluster);
fprintf('\nTop %d data-supported loadings:\n',K);disp(DataSupportedLoadings);
if ~isempty(ScoreCorrelations),fprintf('\nTop %d correlations with realized descriptive latent score:\n',height(ScoreCorrelations));disp(ScoreCorrelations);end
if ~isempty(ContrastDecomposition),fprintf('\nTop standardized descriptive-cluster contrasts with zero/nonzero decomposition:\n');disp(ContrastDecomposition);end
if ~isempty(CrossMethodSummary)
    fprintf('\nSECONDARY TFKM-vs-TRKM comparison (TRKM descriptive, not primary):\n');disp(CrossMethodSummary);
    fprintf('\nCluster cross-tab on the 273 TFKM-retained observations:\n');disp(CrossTab);
    fprintf('\nTFKM trim rate within each TRKM descriptive cluster:\n');disp(TrimByTRKM);
end
fprintf('\nAnonymous outputs written to:\n  %s\n',outroot);
fprintf('No bank names were read or written. No candidate was refitted or reselected.\n');

OUT=struct('summary',Summary,'cluster_sizes',ClusterTable,'subspace_support',SubspaceSupport, ...
    'raw_loadings',RawLoadings,'data_supported_loadings',DataSupportedLoadings, ...
    'score_correlations',ScoreCorrelations,'size_by_cluster',SizeByCluster, ...
    'contrast_decomposition',ContrastDecomposition,'residual_summary',ResidualSummary, ...
    'archived_path',ArchivedPath,'raw_path',RawPath,'selected_state_path',SelectedStatePath, ...
    'reproduction_check',ReproductionCheck,'primary_action',PrimaryAction,'radial_by_cluster',RadialByCluster, ...
    'cross_method_summary',CrossMethodSummary,'cross_tab',CrossTab,'trim_by_trkm',TrimByTRKM, ...
    'output_root',string(outroot));
end

function C=load_candidate_rows_local(CFG)
C=table();
for k=1:150
    f=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));
    d=fullfile(CFG.candidate_root,sprintf('candidate_%03d.COMPLETE.txt',k));
    assert(isfile(f)&&isfile(d),'Candidate %d is incomplete.',k);
    Q=load(f,'Row');
    if isempty(C),C=Q.Row;else,C=[C;Q.Row];end %#ok<AGROW>
end
C=sortrows(C,'candidate_index');
assert(height(C)==150,'Expected 150 candidate rows.');
end

function T=minimal_candidate_rows_local(Q)
T=table();
for i=1:height(Q)
    vars=Q.Properties.VariableNames;
    row=Q(i,:);
    T=[T;table(get_numeric_local(row,vars,'alpha'),get_numeric_local(row,vars,'r'),get_numeric_local(row,vars,'g'), ...
        get_numeric_local(row,vars,'selector_cv_score'),get_numeric_local(row,vars,'selector_score'), ...
        get_numeric_local(row,vars,'objective'),get_numeric_local(row,vars,'W_ratio'),get_numeric_local(row,vars,'R2'), ...
        get_string_local(row,vars,'selected_init_family'),get_numeric_local(row,vars,'selected_start_id'), ...
        'VariableNames',{'alpha','r','g','RV_CV','selector_score','objective','W_ratio','R2','selected_family','selected_start_id'})]; %#ok<AGROW>
end
end

function x=get_numeric_local(row,vars,nm)
if ismember(nm,vars),x=double(row.(nm));else,x=NaN;end
end
function s=get_string_local(row,vars,nm)
if ismember(nm,vars),s=string(row.(nm));else,s="";end
end
function c=pair_corr_local(x,y)
x=double(x(:));y=double(y(:));ok=isfinite(x)&isfinite(y);x=x(ok);y=y(ok);
if numel(x)<2||std(x)==0||std(y)==0,c=NaN;return,end
x=x-mean(x);y=y-mean(y);c=(x'*y)/sqrt((x'*x)*(y'*y));
end
function s=labels_to_strings_local(L,p)
try,s=string(L(:));catch,if iscell(L),s=string(L(:));else,s=compose('V%d',(1:p)');end,end
if numel(s)~=p
    if ischar(L)&&size(L,1)==p,s=string(cellstr(L));else,s=compose('V%d',(1:p)');end
end
s=s(:);
end
