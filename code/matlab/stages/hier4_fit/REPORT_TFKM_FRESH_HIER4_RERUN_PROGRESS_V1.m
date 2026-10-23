function OUT = REPORT_TFKM_FRESH_HIER4_RERUN_PROGRESS_V1(documents_root)
% Completion-only HIER4 rerun progress. No truth metrics.
if nargin<1||strlength(string(documents_root))==0,root=string(fullfile(getenv('USERPROFILE'),'Documents'));else,root=string(documents_root);end
src=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');assert(isfile(fullfile(src,'study.mat')),'Original fresh study not found.');S=load(fullfile(src,'study.mat'),'catalog','CFG');
outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');R=table();
for rep=1:S.CFG.n_replicates
    n=0;for sid=double(S.catalog.scenario_id)',n=n+isfile(fullfile(outroot,'fits','TFKM',sprintf('s%d',sid),sprintf('rep%02d',rep),'COMPLETE.txt'));end
    R=[R;table(rep,n,height(S.catalog),n==height(S.catalog),'VariableNames',{'replicate','N_complete_cells','N_total_cells','replicate_complete'})]; %#ok<AGROW>
end
fprintf('\n============================================================\nTFKM HIER4 CORRECTION-RERUN PROGRESS - NO TRUTH\n============================================================\n');
fprintf('TFKM-HIER4: %d/%d unit cubes complete | full reps=%d/%d\n',sum(R.N_complete_cells),height(S.catalog)*S.CFG.n_replicates,sum(R.replicate_complete),S.CFG.n_replicates);
OUT=struct('progress',R,'output_root',string(outroot));
end
