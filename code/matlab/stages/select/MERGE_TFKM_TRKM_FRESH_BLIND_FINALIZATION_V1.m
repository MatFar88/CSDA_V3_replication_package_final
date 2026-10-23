function OUT = MERGE_TFKM_TRKM_FRESH_BLIND_FINALIZATION_V1(documents_root)
%MERGE_TFKM_TRKM_FRESH_BLIND_FINALIZATION_V1
% Merge and seal operational predictions only after all 2460 unit-method
% finalizations are complete. Does NOT evaluate truth.
if nargin<1,documents_root="";end
[root,outroot,S]=load_study_local(documents_root); %#ok<ASGLU>
base=fullfile(outroot,'blind_finalization_v1');N=height(S.manifest);assert(N==1230,'Manifest size drift.');
PredTF=table();PredTR=table();EdgesTF=table();EdgesTR=table();RepTF=table();RepTR=table();
for method=["TFKM","TRKM"]
    D=dir(fullfile(base,char(method),'s*','rep*','COMPLETE.txt'));assert(numel(D)==N,'%s incomplete: %d/%d.',method,numel(D),N);
    P=table();E=table();R=table();
    for k=1:numel(D)
        udir=D(k).folder;pf=fullfile(udir,'prediction.csv');ef=fullfile(udir,'edges.csv');rf=fullfile(udir,'reproduction.csv');assert(isfile(pf)&&isfile(ef)&&isfile(rf),'Missing unit export in %s',udir);
        p=readtable(pf,'TextType','string');P=[P;p]; %#ok<AGROW>
        e=readtable(ef,'TextType','string');if width(e)>3||height(e)>0,E=[E;e];end %#ok<AGROW>
        r=readtable(rf,'TextType','string');if width(r)>3||height(r)>0,R=[R;r];end %#ok<AGROW>
    end
    P=sortrows(P,{'scenario_id','replicate'});assert(height(P)==N,'Merged %s prediction count mismatch.',method);assert(height(unique(P(:,{'scenario_id','replicate'}),'rows'))==N,'Duplicate %s predictions.',method);assert_no_truth_columns_local(P);assert_no_truth_columns_local(E);assert_no_truth_columns_local(R);
    if method=="TFKM",PredTF=P;EdgesTF=E;RepTF=R;else,PredTR=P;EdgesTR=E;RepTR=R;end
end
seal=fullfile(base,'FROZEN');if ~isfolder(seal),mkdir(seal);end
writetable(PredTF,fullfile(seal,'TFKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv'));writetable(PredTR,fullfile(seal,'TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv'));
if ~isempty(EdgesTF),writetable(EdgesTF,fullfile(seal,'TFKM_BLIND_EDGE_LOG_FROZEN.csv'));end;if ~isempty(EdgesTR),writetable(EdgesTR,fullfile(seal,'TRKM_BLIND_EDGE_LOG_FROZEN.csv'));end
if ~isempty(RepTF),writetable(RepTF,fullfile(seal,'TFKM_BLIND_REPRODUCTION_FROZEN.csv'));end;if ~isempty(RepTR),writetable(RepTR,fullfile(seal,'TRKM_BLIND_REPRODUCTION_FROZEN.csv'));end
fid=fopen(fullfile(seal,'BLIND_FINALIZATION_SEAL.txt'),'w');assert(fid>0);fprintf(fid,'TFKM_TRKM_FRESH_BLIND_FINALIZATION_V1\nsealed_at=%s\nTFKM_predictions=%d\nTRKM_predictions=%d\ntruth_opened=NO\nNo selector threshold or rule was changed after fresh fitting began.\n',char(datetime('now')),height(PredTF),height(PredTR));fclose(fid);
fprintf('\n============================================================\n');fprintf('BLIND FINALIZATION SEALED - TRUTH STILL CLOSED\n');fprintf('TFKM predictions: %d\nTRKM predictions: %d\n',height(PredTF),height(PredTR));fprintf('Seal: %s\n',seal);fprintf('============================================================\n');
OUT=struct('tfkm',PredTF,'trkm',PredTR,'seal_dir',string(seal),'truth_opened',false);
end
function assert_no_truth_columns_local(T)
if isempty(T),return,end;n=lower(string(T.Properties.VariableNames));bad=contains(n,"truth")|contains(n,"posthoc")|contains(n,"outlier")|ismember(n,["r0","g0","c0","alpha0","contamination","rho","delta","p_over_n"]);assert(~any(bad),'TRUTH FIREWALL FAIL: %s',strjoin(n(bad),', '));
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end;assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
