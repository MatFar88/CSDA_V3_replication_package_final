function R = REPORT_TFKM_TRKM_FRESH_BLIND_PROGRESS_V1(documents_root)
%REPORT_TFKM_TRKM_FRESH_BLIND_PROGRESS_V1 Counts completed blind finalizations only.
if nargin<1,documents_root="";end
[root,outroot,S]=load_study_local(documents_root); %#ok<ASGLU>
base=fullfile(outroot,'blind_finalization_v1');
N=height(S.manifest);ct=zeros(2,1);methods=["TFKM","TRKM"];
for m=1:2
    D=dir(fullfile(base,char(methods(m)),'s*','rep*','COMPLETE.txt'));ct(m)=numel(D);
end
fprintf('\n============================================================\n');
fprintf('FRESH BLIND FINALIZATION PROGRESS - NO TRUTH\n');
fprintf('============================================================\n');
fprintf('TFKM: %d/%d complete\n',ct(1),N);fprintf('TRKM: %d/%d complete\n',ct(2),N);
R=table(methods',ct,repmat(N,2,1),'VariableNames',{'method','complete','target'});
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end;assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
