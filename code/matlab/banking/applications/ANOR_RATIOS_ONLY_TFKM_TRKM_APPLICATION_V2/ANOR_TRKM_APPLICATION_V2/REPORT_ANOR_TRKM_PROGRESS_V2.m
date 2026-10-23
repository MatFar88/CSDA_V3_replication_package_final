function T = REPORT_ANOR_TRKM_PROGRESS_V2()
%REPORT_ANOR_TRKM_PROGRESS_V2 Read-only progress report for 150 TRKM candidates.
root=fileparts(mfilename('fullpath'));
CFG=anor_trkm_config_v2(root);
Rows=table();
for k=1:150
    f=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));
    d=fullfile(CFG.candidate_root,sprintf('candidate_%03d.COMPLETE.txt',k));
    complete=isfile(f)&&isfile(d);valid=false;r=NaN;g=NaN;a=NaN;runtime=NaN;fam="";rv=NaN;r2=NaN;
    if complete
        Q=load(f,'Row');R=Q.Row;
        valid=logical(R.valid(1));r=double(R.r(1));g=double(R.g(1));a=double(R.alpha(1));runtime=double(R.runtime_seconds(1));
        fam=string(R.selected_init_family(1));rv=double(R.selector_cv_score(1));r2=double(R.R2(1));
    end
    row=table(k,complete,valid,r,g,a,runtime,fam,rv,r2,'VariableNames', ...
        {'candidate_index','complete','valid','r','g','alpha','runtime_seconds','selected_family','RV_CV','R2'});
    Rows=[Rows;row]; %#ok<AGROW>
end
T=Rows;
fprintf('\nANOR TRKM ratios-only progress: %d/150 complete | %d valid | %d invalid.\n', ...
    sum(T.complete),sum(T.complete&T.valid),sum(T.complete&~T.valid));
if any(T.complete),fprintf('Completed runtime sum: %.2f hours.\n',sum(T.runtime_seconds(T.complete),'omitnan')/3600);end
if isfolder(CFG.output_root),writetable(T,fullfile(CFG.output_root,'CANDIDATE_PROGRESS.csv'));end
end
