function OUT = ANALYZE_ANOR_TFKM_TRKM_V1(topK)
%ANALYZE_ANOR_TFKM_TRKM_V1 Read-only joint interpretation of ANOR TFKM/TRKM.
%
% Run from the common parent folder containing:
%   ANOR_TFKM_APPLICATION_V1
%   ANOR_TRKM_APPLICATION_V1
%
% The function performs NO REFIT and NO RESELECTION. It reads the already
% finalized TFKM descriptive fit and TRKM primary terminal fit, verifies that
% both use the same 363x382 cleaned ratio matrix, and constructs anonymous
% interpretation diagnostics. Bank names are never read or written.
%
% Usage:
%   cd('C:\Users\...\Documents\ANOR_TFKM_TRKM_APPLICATION_V1')
%   A = ANALYZE_ANOR_TFKM_TRKM_V1(25);

if nargin < 1 || isempty(topK), topK = 25; end
assert(isscalar(topK) && topK>=1 && isfinite(topK),'topK must be a positive scalar.');
topK = floor(topK);

here = string(pwd);
[root,tfkmRoot,trkmRoot] = resolve_roots_local(here);

% ---------------------------------------------------------------------
% Exact cleaned inputs used by the two applications.
% ---------------------------------------------------------------------
tfInput = find_single_local(tfkmRoot, '*INPUT_INTERNAL.mat', 'TFKM internal input');
trInput = find_single_local(trkmRoot, '*INPUT_INTERNAL.mat', 'TRKM internal input');
ST = load(tfInput,'X','original_row_index');
SR = load(trInput,'X','original_row_index');
assert(isfield(ST,'X') && isfield(SR,'X'),'Internal input files do not contain X.');
XT = double(ST.X); XR = double(SR.X);
assert(isequal(size(XT),[363 382]),'Unexpected TFKM input size: %dx%d.',size(XT,1),size(XT,2));
assert(isequal(size(XR),[363 382]),'Unexpected TRKM input size: %dx%d.',size(XR,1),size(XR,2));
assert(max(abs(XT(:)-XR(:)))<=1e-12,'TFKM and TRKM cleaned ratio matrices differ.');
X = XT;
idxT = double(ST.original_row_index(:)); idxR = double(SR.original_row_index(:));
assert(isequal(idxT,idxR),'TFKM and TRKM original-row maps differ.');
original_row_index = idxT;

% ---------------------------------------------------------------------
% Recover original ANOR file only to obtain log(total assets), which was
% deliberately excluded from the 382 clustering variables.
% ---------------------------------------------------------------------
rawFile = find_source_local(root,'data_ECB_ANOR.mat');
labFile = find_source_local(root,'labels_ECB_ANOR.mat');
R0 = load(rawFile,'DataECB'); L0 = load(labFile,'labelsPlus');
assert(isfield(R0,'DataECB'),'data_ECB_ANOR.mat lacks DataECB.');
Raw = double(R0.DataECB);
assert(isequal(size(Raw),[365 383]),'Unexpected raw ANOR matrix size.');
assert(all(original_row_index>=1 & original_row_index<=365),'Original-row map outside 1:365.');
logAssets = Raw(original_row_index,1);
XrawRatios = Raw(original_row_index,2:end);
ratio_source_max_abs_diff = max(abs(X(:)-XrawRatios(:)));
assert(ratio_source_max_abs_diff<=1e-10, ...
    'Cleaned 382-variable matrix does not match raw columns 2:383 at mapped rows.');

codes = labels_to_strings_local(L0.labelsPlus);
assert(numel(codes)==382,'Expected 382 ratio labels.');

% Duplicate-row provenance in the 365-row source.
[~,~,ic] = unique(Raw,'rows','stable');
cc = accumarray(ic,1);
dupGroups = cell(0,1);
for k=find(cc>1)'
    dupGroups{end+1,1}=find(ic==k)'; %#ok<AGROW>
end

% ---------------------------------------------------------------------
% TFKM secondary descriptive state and primary metadata.
% ---------------------------------------------------------------------
tfDescFile = find_single_local(tfkmRoot,'*TFKM_DESCRIPTIVE_REFIT.mat','TFKM descriptive refit');
TD = load(tfDescFile,'Descriptive','Labels','Trim','U');
assert(all(isfield(TD,{'Labels','Trim','U'})),'TFKM descriptive file contract mismatch.');
yF = double(TD.Labels(:)); trimF = logical(TD.Trim(:)); UF = double(TD.U);
assert(numel(yF)==363 && numel(trimF)==363,'TFKM descriptive vector length mismatch.');
assert(~any(trimF),'Expected alpha=0 TFKM descriptive fit to retain all observations.');
assert(size(UF,1)==382 && size(UF,2)==2,'Expected TFKM descriptive rank 2.');
UF = orient_columns_local(UF);
ZF = X*UF;

tfFinalFile = find_single_local(tfkmRoot,'*TFKM_PRIMARY_FINAL.mat','TFKM primary final');
TP = load(tfFinalFile,'Prediction');
PredF = TP.Prediction;

% ---------------------------------------------------------------------
% TRKM primary terminal state.
% ---------------------------------------------------------------------
trFinalFile = find_single_local(trkmRoot,'*TRKM_PRIMARY_FINAL.mat','TRKM primary final');
RD = load(trFinalFile,'Prediction','FinalLabels','FinalTrim','FinalU','Audit');
assert(all(isfield(RD,{'Prediction','FinalLabels','FinalTrim','FinalU'})), ...
    'TRKM primary final file contract mismatch.');
PredR = RD.Prediction;
yR = double(RD.FinalLabels(:)); trimR = logical(RD.FinalTrim(:)); UR = double(RD.FinalU);
assert(numel(yR)==363 && numel(trimR)==363,'TRKM terminal vector length mismatch.');
assert(size(UR,1)==382 && size(UR,2)==1,'Expected TRKM terminal rank 1.');
UR = orient_columns_local(UR);
ZR = X*UR;

% ---------------------------------------------------------------------
% Basic selected-state summary.
% ---------------------------------------------------------------------
Summary = table( ...
    string(PredF.action_status), logical(PredF.alpha_emitted), double(PredF.alpha_hat), ...
    double(PredF.terminal_r), double(PredF.terminal_g), logical(PredF.final_fit_available), ...
    string(PredR.action_status), logical(PredR.alpha_emitted), double(PredR.alpha_hat), ...
    double(PredR.terminal_r), double(PredR.terminal_g), logical(PredR.final_fit_available), ...
    sum(~trimR), sum(trimR), ratio_source_max_abs_diff, ...
    'VariableNames',{'TFKM_action','TFKM_alpha_emitted','TFKM_alpha_hat','TFKM_r','TFKM_g','TFKM_primary_fit', ...
    'TRKM_action','TRKM_alpha_emitted','TRKM_alpha_hat','TRKM_r','TRKM_g','TRKM_primary_fit', ...
    'TRKM_retained_n','TRKM_trimmed_n','ratio_source_max_abs_diff'});

% ---------------------------------------------------------------------
% Cluster sizes and latent centroids.
% ---------------------------------------------------------------------
TFKMCluster = cluster_summary_local(yF,false(363,1),ZF,logAssets,'TFKM');
TRKMCluster = cluster_summary_local(yR,trimR,ZR,logAssets,'TRKM');

% Cross-tab and ARI on the TRKM-retained subset.
keepR = ~trimR;
ari_retained = adjusted_rand_local(yF(keepR),yR(keepR));
CrossTab = crosstab_table_local(yF(keepR),yR(keepR));
TrimByTFKM = trim_rate_by_group_local(yF,trimR,'TFKM_cluster');

% ---------------------------------------------------------------------
% Subspace relation and row-space support.
% ---------------------------------------------------------------------
QF = orth(UF); QR = orth(UR);
rank1_overlap_sq = norm(QF'*QR,'fro')^2;
rank1_overlap_sq = min(max(rank1_overlap_sq,0),1);
principal_angle_deg = acosd(sqrt(rank1_overlap_sq));
SubspaceComparison = table(ari_retained,rank1_overlap_sq,principal_angle_deg, ...
    'VariableNames',{'ARI_on_TRKM_retained','TRKM_rank1_in_TFKM_rank2_overlap_squared','principal_angle_degrees'});

RowSupportTFKM = rowspace_support_local(X,false(363,1),UF,'TFKM');
RowSupportTRKM = rowspace_support_local(X,trimR,UR,'TRKM');

% ---------------------------------------------------------------------
% Log-assets association. This variable was not used in clustering.
% ---------------------------------------------------------------------
AssetAssocTFKM = asset_assoc_local(ZF,logAssets,'TFKM');
AssetAssocTRKM = asset_assoc_local(ZR,logAssets,'TRKM');
AssetByTRKMStatus = table( ...
    ["retained";"trimmed"], ...
    [sum(~trimR);sum(trimR)], ...
    [mean(logAssets(~trimR));mean(logAssets(trimR))], ...
    [median(logAssets(~trimR));median(logAssets(trimR))], ...
    'VariableNames',{'status','n','mean_log_total_assets','median_log_total_assets'});

% ---------------------------------------------------------------------
% Ratio-variable interpretation.
% ---------------------------------------------------------------------
TFKMScoreCorrelations = score_correlations_local(X,ZF,codes,'TFKM',topK);
TRKMScoreCorrelations = score_correlations_local(X,ZR,codes,'TRKM',topK);
TFKMLoadings = top_loadings_local(UF,codes,'TFKM',topK);
TRKMLoadings = top_loadings_local(UR,codes,'TRKM',topK);

TFKMContrasts = multigroup_contrasts_local(X,yF,false(363,1),codes,topK);
TRKMContrasts = twogroup_contrasts_local(X,yR,trimR,codes,topK,'TRKM_cluster_2_minus_1');
TrimContrasts = binary_status_contrasts_local(X,trimR,codes,topK);

% Anonymous rows in the tiny TRKM retained cluster.
retLabs = yR(~trimR); uLabs = unique(retLabs(:))';
counts = arrayfun(@(g)sum(retLabs==g),uLabs);
[~,imin] = min(counts); smallLabel = uLabs(imin);
smallMask = (~trimR) & (yR==smallLabel);
TRKMSmallClusterAnonymous = table(find(smallMask),original_row_index(smallMask),logAssets(smallMask), ...
    'VariableNames',{'clean_row_index','original_365_row','log_total_assets'});

% Radial scores, if stored by the finalizer.
RadialSummary = table();
if isfield(RD,'Audit') && isfield(RD.Audit,'radial_q') && numel(RD.Audit.radial_q)==363
    q = double(RD.Audit.radial_q(:));
    RadialSummary = radial_summary_local(q,yR,trimR);
end

% ---------------------------------------------------------------------
% Output.
% ---------------------------------------------------------------------
outdir = fullfile(root,'ANOR_POSTSELECTION_ANALYSIS_V1_output');
if ~isfolder(outdir), mkdir(outdir); end
writetable(Summary,fullfile(outdir,'SUMMARY.csv'));
writetable(TFKMCluster,fullfile(outdir,'TFKM_DESCRIPTIVE_CLUSTER_SUMMARY.csv'));
writetable(TRKMCluster,fullfile(outdir,'TRKM_PRIMARY_CLUSTER_SUMMARY.csv'));
writetable(CrossTab,fullfile(outdir,'TFKM_TRKM_CROSSTAB_ON_TRKM_RETAINED.csv'));
writetable(TrimByTFKM,fullfile(outdir,'TRKM_TRIM_RATE_BY_TFKM_CLUSTER.csv'));
writetable(SubspaceComparison,fullfile(outdir,'SUBSPACE_AND_PARTITION_COMPARISON.csv'));
writetable(RowSupportTFKM,fullfile(outdir,'TFKM_ROWSPACE_SUPPORT.csv'));
writetable(RowSupportTRKM,fullfile(outdir,'TRKM_ROWSPACE_SUPPORT.csv'));
writetable(AssetAssocTFKM,fullfile(outdir,'TFKM_SCORE_ASSET_ASSOCIATION.csv'));
writetable(AssetAssocTRKM,fullfile(outdir,'TRKM_SCORE_ASSET_ASSOCIATION.csv'));
writetable(AssetByTRKMStatus,fullfile(outdir,'TRKM_RETAINED_TRIMMED_ASSETS.csv'));
writetable(TFKMScoreCorrelations,fullfile(outdir,'TFKM_TOP_SCORE_CORRELATIONS.csv'));
writetable(TRKMScoreCorrelations,fullfile(outdir,'TRKM_TOP_SCORE_CORRELATIONS.csv'));
writetable(TFKMLoadings,fullfile(outdir,'TFKM_TOP_LOADINGS.csv'));
writetable(TRKMLoadings,fullfile(outdir,'TRKM_TOP_LOADINGS.csv'));
writetable(TFKMContrasts,fullfile(outdir,'TFKM_TOP_MULTIGROUP_CONTRASTS.csv'));
writetable(TRKMContrasts,fullfile(outdir,'TRKM_TOP_CLUSTER_CONTRASTS.csv'));
writetable(TrimContrasts,fullfile(outdir,'TRKM_TOP_TRIMMED_VS_RETAINED_CONTRASTS.csv'));
writetable(TRKMSmallClusterAnonymous,fullfile(outdir,'TRKM_SMALL_CLUSTER_ANONYMOUS.csv'));
if ~isempty(RadialSummary),writetable(RadialSummary,fullfile(outdir,'TRKM_RADIAL_SUMMARY.csv'));end

Dup = duplicate_table_local(dupGroups);
writetable(Dup,fullfile(outdir,'RAW_DUPLICATE_ROW_GROUPS.csv'));

save(fullfile(outdir,'ANOR_POSTSELECTION_ANALYSIS_V1.mat'), ...
    'Summary','TFKMCluster','TRKMCluster','CrossTab','TrimByTFKM','SubspaceComparison', ...
    'RowSupportTFKM','RowSupportTRKM','AssetAssocTFKM','AssetAssocTRKM','AssetByTRKMStatus', ...
    'TFKMScoreCorrelations','TRKMScoreCorrelations','TFKMLoadings','TRKMLoadings', ...
    'TFKMContrasts','TRKMContrasts','TrimContrasts','TRKMSmallClusterAnonymous','RadialSummary','Dup','-v7.3');

fprintf('\n============================================================\n');
fprintf('ANOR TFKM/TRKM READ-ONLY POST-SELECTION ANALYSIS\n');
fprintf('NO REFIT / NO RESELECTION / NO BANK NAMES\n');
fprintf('============================================================\n');
fprintf('Verified common input: 363 x 382 ratios; raw col 1 is held out as log(total assets).\n');
fprintf('Raw exact duplicate groups: '); print_duplicate_groups_local(dupGroups);
fprintf('\nSelected-state summary:\n'); disp(Summary);
fprintf('\nTFKM descriptive cluster summary:\n'); disp(TFKMCluster);
fprintf('\nTRKM primary retained-cluster summary:\n'); disp(TRKMCluster);
fprintf('\nTFKM vs TRKM on the %d TRKM-retained observations:\n',sum(keepR)); disp(SubspaceComparison); disp(CrossTab);
fprintf('\nTRKM trim rate within each TFKM descriptive cluster:\n'); disp(TrimByTFKM);
fprintf('\nTFKM row-space support:\n'); disp(RowSupportTFKM);
fprintf('\nTRKM row-space support:\n'); disp(RowSupportTRKM);
fprintf('\nLatent-score association with held-out log(total assets):\n'); disp(AssetAssocTFKM); disp(AssetAssocTRKM);
fprintf('\nLog(total assets) by TRKM retained/trimmed status:\n'); disp(AssetByTRKMStatus);
fprintf('\nTop %d TFKM score correlations:\n',topK); disp(TFKMScoreCorrelations);
fprintf('\nTop %d TRKM score correlations:\n',topK); disp(TRKMScoreCorrelations);
fprintf('\nTop %d TFKM multigroup contrasts:\n',topK); disp(TFKMContrasts);
fprintf('\nTop %d TRKM two-group contrasts:\n',topK); disp(TRKMContrasts);
fprintf('\nTop %d TRKM trimmed-vs-retained contrasts:\n',topK); disp(TrimContrasts);
fprintf('\nAnonymous rows in the smaller TRKM retained cluster (n=%d):\n',height(TRKMSmallClusterAnonymous)); disp(TRKMSmallClusterAnonymous);
if ~isempty(RadialSummary),fprintf('\nTRKM radial score summary:\n');disp(RadialSummary);end
fprintf('\nOutputs written to:\n  %s\n',outdir);
fprintf('No bank names were read or written.\n');

OUT = struct('summary',Summary,'tfkm_cluster',TFKMCluster,'trkm_cluster',TRKMCluster, ...
    'cross_tab',CrossTab,'trim_by_tfkm',TrimByTFKM,'subspace',SubspaceComparison, ...
    'row_support_tfkm',RowSupportTFKM,'row_support_trkm',RowSupportTRKM, ...
    'asset_tfkm',AssetAssocTFKM,'asset_trkm',AssetAssocTRKM,'trim_assets',AssetByTRKMStatus, ...
    'tfkm_score_corr',TFKMScoreCorrelations,'trkm_score_corr',TRKMScoreCorrelations, ...
    'tfkm_loadings',TFKMLoadings,'trkm_loadings',TRKMLoadings, ...
    'tfkm_contrasts',TFKMContrasts,'trkm_contrasts',TRKMContrasts,'trim_contrasts',TrimContrasts, ...
    'small_trkm_cluster',TRKMSmallClusterAnonymous,'radial',RadialSummary,'output_root',string(outdir));
end

% =====================================================================
% Helpers
% =====================================================================
function [root,tf,tr] = resolve_roots_local(here)
cands = [here; string(fileparts(here))];
root="";tf="";tr="";
for c=cands(:)'
    a=fullfile(c,'ANOR_TFKM_APPLICATION_V1'); b=fullfile(c,'ANOR_TRKM_APPLICATION_V1');
    if isfolder(a)&&isfolder(b),root=c;tf=string(a);tr=string(b);return;end
end
error('Could not find sibling ANOR_TFKM_APPLICATION_V1 and ANOR_TRKM_APPLICATION_V1. Run from their common parent.');
end

function f = find_single_local(root,pattern,label)
D=dir(fullfile(root,'**',pattern));D=D(~[D.isdir]);
assert(~isempty(D),'Cannot locate %s under %s.',label,root);
% Prefer output files over any copied source/reference file.
paths=string(fullfile({D.folder},{D.name}));
score=contains(lower(paths),'output');
if any(score),paths=paths(score);end
assert(numel(paths)==1,'Expected one %s; found %d.',label,numel(paths));
f=paths(1);
end

function f = find_source_local(root,name)
cands=[string(fullfile(root,name)); string(fullfile(fileparts(root),name)); ...
       string(fullfile(getenv('USERPROFILE'),'Documents',name)); ...
       string(fullfile(root,'ANOR_TFKM_APPLICATION_V1',name)); ...
       string(fullfile(root,'ANOR_TRKM_APPLICATION_V1',name))];
f="";
for c=cands(:)',if isfile(c),f=c;return;end,end
D=dir(fullfile(root,'**',name));
if ~isempty(D),f=string(fullfile(D(1).folder,D(1).name));return;end
error('Cannot locate required source file %s. Put it in the common ANOR parent folder.',name);
end

function s=labels_to_strings_local(x)
x=x(:);s=strings(numel(x),1);
for i=1:numel(x)
    z=x{i};
    if iscell(z),z=z{1};end
    s(i)=string(z);
end
end

function U=orient_columns_local(U)
for j=1:size(U,2)
    [~,k]=max(abs(U(:,j)));
    if U(k,j)<0,U(:,j)=-U(:,j);end
end
end

function T=cluster_summary_local(y,trim,Z,assets,method)
keep=~trim; gs=unique(y(keep));T=table();
for g=gs(:)'
    m=keep & y==g; vals={string(method),double(g),sum(m),sum(m)/sum(keep),mean(assets(m)),median(assets(m))};
    row=table(vals{1},vals{2},vals{3},vals{4},vals{5},vals{6}, ...
        'VariableNames',{'method','cluster','n','retained_fraction','mean_log_total_assets','median_log_total_assets'});
    for j=1:size(Z,2),row.(sprintf('latent_centroid_%d',j))=mean(Z(m,j));end
    T=[T;row]; %#ok<AGROW>
end
end

function T=crosstab_table_local(a,b)
ga=unique(a);gb=unique(b);T=table();
for i=1:numel(ga),for j=1:numel(gb)
    T=[T;table(ga(i),gb(j),sum(a==ga(i)&b==gb(j)), ...
        'VariableNames',{'TFKM_cluster','TRKM_cluster','n_on_TRKM_retained'})]; %#ok<AGROW>
end,end
end

function T=trim_rate_by_group_local(y,trim,name)
gs=unique(y);T=table();
for g=gs(:)'
    m=y==g;row=table(double(g),sum(m),sum(trim(m)),mean(trim(m)), ...
        'VariableNames',{name,'n_total','n_TRKM_trimmed','TRKM_trim_rate'});T=[T;row]; %#ok<AGROW>
end
end

function ari=adjusted_rand_local(a,b)
a=a(:);b=b(:);assert(numel(a)==numel(b));n=numel(a);
[~,~,ia]=unique(a);[~,~,ib]=unique(b);C=accumarray([ia ib],1);
comb2=@(x)x.*(x-1)/2;
sumij=sum(comb2(C(:)));sa=sum(comb2(sum(C,2)));sb=sum(comb2(sum(C,1)'));
tot=comb2(n);expected=sa*sb/tot;mx=.5*(sa+sb);den=mx-expected;
if abs(den)<eps,ari=1;else,ari=(sumij-expected)/den;end
end

function T=rowspace_support_local(X,trim,U,method)
R=X(~trim,:);Rc=R-mean(R,1);[~,S,V]=svd(Rc,'econ');s=diag(S);tol=max(size(Rc))*eps(max([s;1]));rk=sum(s>tol);Vr=V(:,1:rk);
T=table();
for j=1:size(U,2)
    u=U(:,j);u=u/norm(u);mass=sum((Vr'*u).^2);mass=min(max(mass,0),1);
    T=[T;table(string(method),j,size(R,1),size(X,2),rk,size(X,2)-rk,mass,1-mass, ...
        'VariableNames',{'method','component','retained_n','p','rank_centered_X','ambient_null_dimension','rowspace_mass','nullspace_mass'})]; %#ok<AGROW>
end
end

function T=asset_assoc_local(Z,a,method)
T=table();
for j=1:size(Z,2)
    r=corr(Z(:,j),a,'Rows','complete');T=[T;table(string(method),j,r,r^2, ...
        'VariableNames',{'method','component','corr_latent_score_log_total_assets','R2_latent_score_log_total_assets'})]; %#ok<AGROW>
end
end

function T=score_correlations_local(X,Z,codes,method,K)
T=table();
for j=1:size(Z,2)
    rr=nan(size(X,2),1);for q=1:size(X,2),rr(q)=corr(X(:,q),Z(:,j),'Rows','complete');end
    [~,o]=sort(abs(rr),'descend','MissingPlacement','last');o=o(1:min(K,numel(o)));
    A=table(repmat(string(method),numel(o),1),repmat(j,numel(o),1),(1:numel(o))',o(:),codes(o),rr(o),abs(rr(o)), ...
        'VariableNames',{'method','component','rank_abs_corr','variable_index','indicator_code','corr_with_latent_score','abs_corr'});T=[T;A]; %#ok<AGROW>
end
end

function T=top_loadings_local(U,codes,method,K)
T=table();
for j=1:size(U,2)
    [~,o]=sort(abs(U(:,j)),'descend');o=o(1:min(K,numel(o)));
    A=table(repmat(string(method),numel(o),1),repmat(j,numel(o),1),(1:numel(o))',o(:),codes(o),U(o,j),abs(U(o,j)), ...
        'VariableNames',{'method','component','rank_abs','variable_index','indicator_code','loading','abs_loading'});T=[T;A]; %#ok<AGROW>
end
end

function T=multigroup_contrasts_local(X,y,trim,codes,K)
keep=~trim;gs=unique(y(keep));p=size(X,2);best=zeros(p,1);g1=zeros(p,1);g2=zeros(p,1);m1=zeros(p,1);m2=zeros(p,1);nz1=zeros(p,1);nz2=zeros(p,1);
for q=1:p
    v=X(:,q);bestq=-Inf;
    for i=1:numel(gs),for j=i+1:numel(gs)
        a=v(keep&y==gs(i));b=v(keep&y==gs(j));sp=sqrt((var(a,1)+var(b,1))/2);d=0;if sp>0,d=(mean(b)-mean(a))/sp;end
        if abs(d)>abs(bestq) || ~isfinite(bestq),bestq=d;g1(q)=gs(i);g2(q)=gs(j);m1(q)=mean(a);m2(q)=mean(b);nz1(q)=mean(a~=0);nz2(q)=mean(b~=0);end
    end,end
    best(q)=bestq;
end
[~,o]=sort(abs(best),'descend');o=o(1:min(K,p));
T=table((1:numel(o))',o(:),codes(o),g1(o),g2(o),best(o),abs(best(o)),m1(o),m2(o),nz1(o),nz2(o),nz2(o)-nz1(o), ...
    'VariableNames',{'rank_abs_std_diff','variable_index','indicator_code','cluster_1','cluster_2','standardized_difference_2_minus_1','abs_standardized_difference','mean_cluster_1','mean_cluster_2','nonzero_rate_cluster_1','nonzero_rate_cluster_2','nonzero_rate_difference_2_minus_1'});
end

function T=twogroup_contrasts_local(X,y,trim,codes,K,label)
keep=~trim;gs=unique(y(keep));assert(numel(gs)==2,'Expected two retained groups.');aidx=keep&y==gs(1);bidx=keep&y==gs(2);T=contrast_core_local(X,aidx,bidx,codes,K,label);T.cluster_1=repmat(gs(1),height(T),1);T.cluster_2=repmat(gs(2),height(T),1);
end

function T=binary_status_contrasts_local(X,trim,codes,K)
T=contrast_core_local(X,~trim,trim,codes,K,'trimmed_minus_retained');
end

function T=contrast_core_local(X,aidx,bidx,codes,K,label)
p=size(X,2);d=zeros(p,1);ma=zeros(p,1);mb=zeros(p,1);nza=zeros(p,1);nzb=zeros(p,1);
for q=1:p
    a=X(aidx,q);b=X(bidx,q);sp=sqrt((var(a,1)+var(b,1))/2);if sp>0,d(q)=(mean(b)-mean(a))/sp;end
    ma(q)=mean(a);mb(q)=mean(b);nza(q)=mean(a~=0);nzb(q)=mean(b~=0);
end
[~,o]=sort(abs(d),'descend');o=o(1:min(K,p));
T=table(repmat(string(label),numel(o),1),(1:numel(o))',o(:),codes(o),d(o),abs(d(o)),ma(o),mb(o),nza(o),nzb(o),nzb(o)-nza(o), ...
    'VariableNames',{'contrast','rank_abs_std_diff','variable_index','indicator_code','standardized_difference_B_minus_A','abs_standardized_difference','mean_A','mean_B','nonzero_rate_A','nonzero_rate_B','nonzero_rate_difference_B_minus_A'});
end

function T=radial_summary_local(q,y,trim)
T=table();
status=["retained";"trimmed"];
for s=1:2
    if s==1,m=~trim;else,m=trim;end
    row=table(status(s),sum(m),mean(q(m)),median(q(m)),quantile(q(m),.9),quantile(q(m),.95),max(q(m)), ...
        'VariableNames',{'status','n','mean_radial_q','median_radial_q','q90_radial_q','q95_radial_q','max_radial_q'});T=[T;row]; %#ok<AGROW>
end
gs=unique(y(~trim));
for g=gs(:)'
    m=(~trim)&y==g;row=table("retained_cluster_"+string(g),sum(m),mean(q(m)),median(q(m)),quantile(q(m),.9),quantile(q(m),.95),max(q(m)), ...
        'VariableNames',{'status','n','mean_radial_q','median_radial_q','q90_radial_q','q95_radial_q','max_radial_q'});T=[T;row]; %#ok<AGROW>
end
end

function T=duplicate_table_local(groups)
T=table();
for k=1:numel(groups)
    v=groups{k};T=[T;table(k,string(sprintf('%d ',v)),numel(v), ...
        'VariableNames',{'duplicate_group','raw_rows','multiplicity'})]; %#ok<AGROW>
end
end

function print_duplicate_groups_local(groups)
for k=1:numel(groups)
    fprintf('{');fprintf('%d',groups{k}(1));for j=2:numel(groups{k}),fprintf(',%d',groups{k}(j));end;fprintf('}');if k<numel(groups),fprintf(' ');end
end
end
