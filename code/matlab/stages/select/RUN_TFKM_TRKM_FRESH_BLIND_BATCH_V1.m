function OUT = RUN_TFKM_TRKM_FRESH_BLIND_BATCH_V1(replicates,scenario_ids,method,documents_root)
%RUN_TFKM_TRKM_FRESH_BLIND_BATCH_V1 Restartable truth-blind finalization batch.
% Example:
%   RUN_TFKM_TRKM_FRESH_BLIND_BATCH_V1(1:10,[],"REG_TFKM_HD");
if nargin<1||isempty(replicates),replicates=1:30;end
if nargin<2,scenario_ids=[];end
if nargin<3,error('Provide method.');end
if nargin<4,documents_root="";end
method=string(method);assert(any(method==["REG_TFKM_HD","REG_TRKM_HD"]));
[root,~,S]=load_study_local(documents_root); %#ok<ASGLU>
allids=double(S.catalog.scenario_id(:))';if isempty(scenario_ids),scenario_ids=allids;else,scenario_ids=double(scenario_ids(:))';assert(all(ismember(scenario_ids,allids)),'Unknown scenario id.');end
replicates=double(replicates(:))';assert(all(ismember(replicates,1:30)),'Bad replicate.');
total=numel(replicates)*numel(scenario_ids);done=0;
fprintf('BLIND BATCH START: %d unit-method finalizations | %s | truth sealed.\n',total,method);
for rep=replicates
    for sid=scenario_ids
        RUN_TFKM_TRKM_FRESH_BLIND_UNIT_V1(sid,rep,method,root);done=done+1;
        if mod(done,10)==0||done==total,fprintf('BLIND BATCH PROGRESS %d/%d | %s\n',done,total,method);end
    end
end
OUT=struct('complete',true,'N',total,'method',method);
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end;assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
