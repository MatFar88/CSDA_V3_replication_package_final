function OUT = REPORT_TFKM_TRKM_FRESH_PROGRESS_V1(documents_root)
%REPORT_TFKM_TRKM_FRESH_PROGRESS_V1 Completion-only report. No truth metrics.
if nargin<1||strlength(string(documents_root))==0,root=string(fullfile(getenv('USERPROFILE'),'Documents'));else,root=string(documents_root);end
outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');assert(isfile(fullfile(outroot,'study.mat')),'Study not built.');S=load(fullfile(outroot,'study.mat'),'catalog','CFG');
methods=["TFKM","TRKM"];R=table();
for m=methods
    for rep=1:S.CFG.n_replicates
        n=0;
        for sid=double(S.catalog.scenario_id)'
            f=fullfile(outroot,'fits',char(m),sprintf('s%d',sid),sprintf('rep%02d',rep),'COMPLETE.txt');n=n+isfile(f);
        end
        R=[R;table(m,rep,n,height(S.catalog),n==height(S.catalog),'VariableNames',{'method','replicate','N_complete_cells','N_total_cells','replicate_complete'})]; %#ok<AGROW>
    end
end
fprintf('\n============================================================\n');fprintf('FRESH GENERALIZATION FIT PROGRESS - NO TRUTH\n');fprintf('============================================================\n');
for m=methods
    Z=R(R.method==m,:);fprintf('%s: %d/%d unit cubes complete | full reps=%d/%d\n',m,sum(Z.N_complete_cells),height(S.catalog)*S.CFG.n_replicates,sum(Z.replicate_complete),S.CFG.n_replicates);
end
OUT=struct('progress',R,'output_root',string(outroot));
end
