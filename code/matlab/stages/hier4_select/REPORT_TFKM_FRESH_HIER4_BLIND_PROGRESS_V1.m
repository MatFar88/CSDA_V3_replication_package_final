function R = REPORT_TFKM_FRESH_HIER4_BLIND_PROGRESS_V1(documents_root)
%REPORT_TFKM_FRESH_HIER4_BLIND_PROGRESS_V1 Read-only progress; no truth.
if nargin<1,documents_root="";end
[root,S]=load_study_local(documents_root);N=height(S.manifest);
base=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output','blind_finalization_hier4_v1','TFKM');
D=dir(fullfile(base,'s*','rep*','COMPLETE.txt'));ct=numel(D);
fullreps=0;
ids=double(S.catalog.scenario_id(:));
for rep=1:30
    ok=true;
    for k=1:numel(ids)
        if ~isfile(fullfile(base,sprintf('s%d',ids(k)),sprintf('rep%02d',rep),'COMPLETE.txt')),ok=false;break,end
    end
    fullreps=fullreps+ok;
end
fprintf('\n============================================================\n');
fprintf('TFKM HIER4 BLIND FINALIZATION PROGRESS - NO TRUTH\n');
fprintf('============================================================\n');
fprintf('TFKM-HIER4 BLIND: %d/%d complete | full reps=%d/30\n',ct,N,fullreps);
R=table(ct,N,fullreps,'VariableNames',{'complete','target','full_replicates'});
end
function [root,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";
for c=cands(:)'
    orig=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');hier=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
    if isfile(fullfile(orig,'study.mat'))&&isfolder(hier),root=c;S=load(fullfile(orig,'study.mat'));break,end
end
assert(strlength(root)>0,'Original fresh study + HIER4 correction-rerun output not found.');
end
