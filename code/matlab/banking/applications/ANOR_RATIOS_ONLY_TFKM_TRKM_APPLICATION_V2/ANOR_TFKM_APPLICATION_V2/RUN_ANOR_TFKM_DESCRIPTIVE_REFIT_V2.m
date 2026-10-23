function OUT = RUN_ANOR_TFKM_DESCRIPTIVE_REFIT_V2()
%RUN_ANOR_TFKM_DESCRIPTIVE_REFIT_V2 Pre-specified secondary descriptive fit.
%
% This is NOT part of the primary selector output. It is used only when the
% primary procedure retains a unique structural pair but does not define a
% terminal partition. If no positive trimming action was selected, alpha=0.
% If positive trimming was selected, the emitted alpha is used.

root=fileparts(mfilename('fullpath'));
CFG=anor_tfkm_config_v2(root);
anor_tfkm_path_lock_v2(root);
primary_file=fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_FINAL.mat');
assert(isfile(primary_file),'Run MERGE_FINALIZE_ANOR_TFKM_V2 first.');
P=load(primary_file,'Prediction');
Pred=P.Prediction;

if logical(Pred.final_fit_available)
    fprintf('Primary selector already defines a terminal fit; no descriptive refit is needed.\n');
    OUT=struct('run',false,'reason',"PRIMARY_FIT_AVAILABLE");
    return
end
assert(~logical(Pred.action_ambiguous),'Descriptive refit forbidden after ambiguous trimming action.');
assert(logical(Pred.terminal_unique),'Descriptive refit requires a unique retained structural pair.');
r=double(Pred.terminal_r);g=double(Pred.terminal_g);
assert(isfinite(r)&&isfinite(g),'Unique structural pair is not numerically available.');
if logical(Pred.alpha_emitted)
    a=double(Pred.alpha_hat);
    alpha_source="PRIMARY_EMITTED_ALPHA";
else
    a=0;
    alpha_source="NO_ACTION_ALPHA_ZERO";
end

S=load(CFG.input_file,'X','options','pair_seed','original_row_index');X=double(S.X);options=S.options;original_row_index=double(S.original_row_index(:));
assert(double(S.pair_seed)==CFG.pair_seed,'Pair-seed drift.');

fprintf('\n============================================================\n');
fprintf('ANOR TFKM SECONDARY DESCRIPTIVE REFIT\n');
fprintf('NOT a primary selector output\n');
fprintf('state=(r,g,alpha)=(%d,%d,%.3f) | %s\n',r,g,a,alpha_source);
fprintf('============================================================\n');

t0=tic;
fit=run_hd_multistart_top_hier_stable4_unseen_v1(X,r,g,a,options,CFG.pair_seed);
runtime_seconds=toc(t0);
assert(fit.valid,'Descriptive refit returned an invalid fit: %s',string(fit.message));
st=compute_external_alpha_fit_stats_v1(X,fit,"REG_TFKM_HD");

if ~isfolder(CFG.descriptive_root),mkdir(CFG.descriptive_root);end
Descriptive=table("TFKM",CFG.dataset,r,g,a,string(alpha_source), ...
    string(fit.selected_init_family),double(fit.selected_start_id),double(fit.hier_cross_H_score), ...
    double(fit.hier_cross_ari_score),double(fit.hier_cross_j_score),double(fit.objective), ...
    double(st.W_ratio),double(st.R2),double(st.explained),runtime_seconds,true, ...
    'VariableNames',{'method','dataset','r','g','alpha','alpha_source','selected_family','selected_start_id', ...
    'H','A','J','objective','W_ratio','R2','explained','runtime_seconds','descriptive_only'});
Labels=double(fit.labels(:));Trim=logical(fit.trim(:));U=double(fit.U);
save(fullfile(CFG.descriptive_root,'ANOR_TFKM_DESCRIPTIVE_REFIT.mat'), ...
    'Descriptive','Labels','Trim','U','fit','-v7.3');
writetable(Descriptive,fullfile(CFG.descriptive_root,'ANOR_TFKM_DESCRIPTIVE_REFIT_SUMMARY.csv'));
Anon=table((1:CFG.n)',original_row_index,Labels,Trim,'VariableNames',{'bank_index','original_365_row','cluster','trimmed'});
writetable(Anon,fullfile(CFG.descriptive_root,'ANOR_TFKM_DESCRIPTIVE_PARTITION_ANONYMOUS.csv'));

fprintf('Descriptive refit complete. H=%.6f | family=%s | %.1fs\n', ...
    double(fit.hier_cross_H_score),string(fit.selected_init_family),runtime_seconds);
OUT=struct('run',true,'summary',Descriptive,'output_root',CFG.descriptive_root);
end
