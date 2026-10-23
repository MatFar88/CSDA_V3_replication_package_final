function OUT = COMPARE_ECB_TFKM_TRKM_PRIMARY_V1(tfkm_root)
%COMPARE_ECB_TFKM_TRKM_PRIMARY_V1 Read-only primary-output comparison.
% No refits, no reselection, no bank names.

root=fileparts(mfilename('fullpath'));
ecb_trkm_path_lock_v1(root);
if nargin<1 || strlength(string(tfkm_root))==0
    tfkm_root=fullfile(fileparts(root),'ECB_TFKM_APPLICATION_V1');
end
tfkm_root=char(string(tfkm_root));
trfile=fullfile(root,'output','final','ECB_TRKM_PRIMARY_FINAL.mat');
tffile=fullfile(tfkm_root,'output','final','ECB_TFKM_PRIMARY_FINAL.mat');
assert(isfile(trfile),'TRKM primary final file not found. Run FINALIZE first.');
assert(isfile(tffile),'TFKM primary final file not found: %s',tffile);

TR=load(trfile,'Prediction','FinalLabels','FinalTrim','FinalU');
TF=load(tffile,'Prediction','FinalLabels','FinalTrim','FinalU');
P1=TF.Prediction;P2=TR.Prediction;
Summary=table("TFKM","TRKM",string(P1.action_status),string(P2.action_status), ...
    logical(P1.alpha_emitted),logical(P2.alpha_emitted),double(P1.alpha_hat),double(P2.alpha_hat), ...
    double(P1.terminal_r),double(P2.terminal_r),double(P1.terminal_g),double(P2.terminal_g), ...
    logical(P1.final_fit_available),logical(P2.final_fit_available), ...
    'VariableNames',{'method_1','method_2','TFKM_action','TRKM_action','TFKM_alpha_emitted','TRKM_alpha_emitted', ...
    'TFKM_alpha_hat','TRKM_alpha_hat','TFKM_terminal_r','TRKM_terminal_r','TFKM_terminal_g','TRKM_terminal_g', ...
    'TFKM_primary_fit','TRKM_primary_fit'});

FitComparison=table();SubspaceComparison=table();
if logical(P1.final_fit_available) && logical(P2.final_fit_available) && ...
        isfield(TF,'FinalLabels') && isfield(TR,'FinalLabels') && ~isempty(TF.FinalLabels) && ~isempty(TR.FinalLabels)
    t1=logical(TF.FinalTrim(:));t2=logical(TR.FinalTrim(:));
    l1=double(TF.FinalLabels(:));l2=double(TR.FinalLabels(:));
    assert(numel(t1)==numel(t2),'Primary fits have different sample sizes.');
    keep=~t1&~t2;
    if sum(keep)>1,ari=adjusted_rand_index_v1(l1(keep),l2(keep));else,ari=NaN;end
    un=t1|t2;if any(un),jac=sum(t1&t2)/sum(un);else,jac=1;end
    FitComparison=table(numel(t1),sum(~t1),sum(~t2),sum(keep),sum(t1),sum(t2),sum(t1&t2),sum(t1|t2),jac,ari, ...
        'VariableNames',{'n','TFKM_retained','TRKM_retained','mutually_retained','TFKM_trimmed','TRKM_trimmed', ...
        'trim_intersection','trim_union','trim_Jaccard','ARI_mutually_retained'});
end

if logical(P1.final_fit_available) && logical(P2.final_fit_available) && ...
        isfield(TF,'FinalU') && isfield(TR,'FinalU') && ~isempty(TF.FinalU) && ~isempty(TR.FinalU)
    U1=orth(double(TF.FinalU));U2=orth(double(TR.FinalU));
    sv=svd(U1'*U2);sv=min(max(sv,0),1);ang=acos(sv)*180/pi;
    overlap=sum(sv.^2)/max(1,min(size(U1,2),size(U2,2)));
    Pgap=norm(U1*U1'-U2*U2',2);
    SubspaceComparison=table(size(U1,2),size(U2,2),min(ang),max(ang),mean(ang),overlap,Pgap, ...
        'VariableNames',{'TFKM_rank','TRKM_rank','min_principal_angle_deg','max_principal_angle_deg', ...
        'mean_principal_angle_deg','normalized_projector_overlap','projector_operator_gap'});
end

outdir=fullfile(root,'output','comparison');if ~isfolder(outdir),mkdir(outdir);end
writetable(Summary,fullfile(outdir,'ECB_TFKM_TRKM_PRIMARY_STRUCTURE.csv'));
if ~isempty(FitComparison),writetable(FitComparison,fullfile(outdir,'ECB_TFKM_TRKM_PRIMARY_PARTITION_COMPARISON.csv'));end
if ~isempty(SubspaceComparison),writetable(SubspaceComparison,fullfile(outdir,'ECB_TFKM_TRKM_PRIMARY_SUBSPACE_COMPARISON.csv'));end

fprintf('\n============================================================\n');
fprintf('ECB PRIMARY TFKM/TRKM READ-ONLY COMPARISON\n');
fprintf('NO REFITS / NO RESELECTION / NO BANK NAMES\n');
fprintf('============================================================\n');
disp(Summary);if ~isempty(FitComparison),disp(FitComparison);end;if ~isempty(SubspaceComparison),disp(SubspaceComparison);end
OUT=struct('summary',Summary,'fit_comparison',FitComparison,'subspace_comparison',SubspaceComparison,'output_root',string(outdir));
end
