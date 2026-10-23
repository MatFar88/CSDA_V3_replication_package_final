function OUT = MERGE_TFKM_HIER4_WITH_TRKM_BLIND_SEAL_V2(documents_root)
%MERGE_TFKM_HIER4_WITH_TRKM_BLIND_SEAL_V2
% Creates a NEW combined frozen seal. It never overwrites the original V1
% seal. Corrected TFKM comes from the HIER4 recovery finalization; TRKM comes
% from the already-completed original blind seal. No truth is evaluated.
if nargin<1,documents_root="";end
[root,S]=load_study_local(documents_root);N=height(S.manifest);assert(N==1230,'Manifest size drift.');
orig=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
hier=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
tfbase=fullfile(hier,'blind_finalization_hier4_v1','TFKM');
trseal=fullfile(orig,'blind_finalization_v1','FROZEN');
trpred=fullfile(trseal,'TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv');assert(isfile(trpred),'Original frozen TRKM prediction file missing.');
% Merge corrected TFKM per-unit predictions/logs.
D=dir(fullfile(tfbase,'s*','rep*','COMPLETE.txt'));assert(numel(D)==N,'Corrected TFKM blind finalization incomplete: %d/%d.',numel(D),N);
PT=table();ET=table();RT=table();
for k=1:numel(D)
    udir=D(k).folder;pf=fullfile(udir,'prediction.csv');ef=fullfile(udir,'edges.csv');rf=fullfile(udir,'reproduction.csv');assert(isfile(pf)&&isfile(ef)&&isfile(rf),'Missing corrected TFKM export in %s',udir);
    p=readtable(pf,'TextType','string');PT=[PT;p]; %#ok<AGROW>
    e=readtable(ef,'TextType','string');if width(e)>3||height(e)>0,ET=[ET;e];end %#ok<AGROW>
    r=readtable(rf,'TextType','string');if width(r)>3||height(r)>0,RT=[RT;r];end %#ok<AGROW>
end
PT=sortrows(PT,{'scenario_id','replicate'});assert(height(PT)==N&&height(unique(PT(:,{'scenario_id','replicate'}),'rows'))==N,'Corrected TFKM prediction key/count mismatch.');assert_no_truth_columns_local(PT);assert_no_truth_columns_local(ET);assert_no_truth_columns_local(RT);
PR=readtable(trpred,'TextType','string');PR=sortrows(PR,{'scenario_id','replicate'});assert(height(PR)==N&&height(unique(PR(:,{'scenario_id','replicate'}),'rows'))==N,'Frozen TRKM prediction key/count mismatch.');assert_no_truth_columns_local(PR);
% Require exact key alignment.
assert(isequal(double(PT.scenario_id),double(PR.scenario_id))&&isequal(double(PT.replicate),double(PR.replicate)),'TFKM/TRKM key alignment failure.');
seal=fullfile(orig,'blind_finalization_v2','FROZEN');if ~isfolder(seal),mkdir(seal);end
writetable(PT,fullfile(seal,'TFKM_HIER4_BLIND_FINAL_PREDICTIONS_FROZEN.csv'));
writetable(PR,fullfile(seal,'TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv'));
if ~isempty(ET),writetable(ET,fullfile(seal,'TFKM_HIER4_BLIND_EDGE_LOG_FROZEN.csv'));end
if ~isempty(RT),writetable(RT,fullfile(seal,'TFKM_HIER4_BLIND_REPRODUCTION_FROZEN.csv'));end
% Carry forward already-frozen TRKM operational logs without recomputation.
for nm=["TRKM_BLIND_EDGE_LOG_FROZEN.csv","TRKM_BLIND_REPRODUCTION_FROZEN.csv"]
    src=fullfile(trseal,char(nm));if isfile(src),copyfile(src,fullfile(seal,char(nm)),'f');end
end
fid=fopen(fullfile(seal,'BLIND_FINALIZATION_SEAL_V2.txt'),'w');assert(fid>0);
fprintf(fid,'TFKM_TRKM_FRESH_BLIND_FINALIZATION_V2\n');
fprintf(fid,'sealed_at=%s\n',char(datetime('now')));
fprintf(fid,'TFKM_source=HIER4_CORRECTION_RERUN + blind_finalization_hier4_v1\n');
fprintf(fid,'TFKM_legacy_SIGN_run_excluded=YES\n');
fprintf(fid,'TRKM_source=original blind_finalization_v1 frozen seal\n');
fprintf(fid,'TFKM_predictions=%d\nTRKM_predictions=%d\n',height(PT),height(PR));
fprintf(fid,'truth_opened=NO\n');
fprintf(fid,'No selector threshold or statistical rule was changed after fresh fitting began; TFKM rerun corrected a pretruth implementation mismatch only.\n');
fclose(fid);
fprintf('\n============================================================\n');
fprintf('COMBINED BLIND FINALIZATION V2 SEALED - TRUTH STILL CLOSED\n');
fprintf('Corrected TFKM predictions: %d\nFrozen TRKM predictions   : %d\n',height(PT),height(PR));
fprintf('Seal: %s\n',seal);fprintf('============================================================\n');
OUT=struct('tfkm',PT,'trkm',PR,'seal_dir',string(seal),'truth_opened',false);
end
function assert_no_truth_columns_local(T)
if isempty(T),return,end;n=lower(string(T.Properties.VariableNames));bad=contains(n,"truth")|contains(n,"posthoc")|contains(n,"outlier")|ismember(n,["r0","g0","c0","alpha0","contamination","rho","delta","p_over_n"]);assert(~any(bad),'TRUTH FIREWALL FAIL: %s',strjoin(n(bad),', '));
end
function [root,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";
for c=cands(:)'
    orig=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');hier=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
    if isfile(fullfile(orig,'study.mat'))&&isfolder(hier),root=c;S=load(fullfile(orig,'study.mat'));break,end
end
assert(strlength(root)>0,'Original fresh study + HIER4 correction-rerun output not found.');
end
