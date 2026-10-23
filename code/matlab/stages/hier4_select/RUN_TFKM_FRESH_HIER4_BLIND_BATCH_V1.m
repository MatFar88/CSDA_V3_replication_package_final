function OUT = RUN_TFKM_FRESH_HIER4_BLIND_BATCH_V1(replicates,scenario_ids,documents_root)
%RUN_TFKM_FRESH_HIER4_BLIND_BATCH_V1 Restartable TFKM-only blind finalization.
% Uses only corrected HIER4 fresh cubes. Truth remains sealed.
if nargin<1||isempty(replicates),replicates=1:30;end
if nargin<2,scenario_ids=[];end
if nargin<3,documents_root="";end
[root,S]=load_study_local(documents_root);
allids=double(S.catalog.scenario_id(:))';
if isempty(scenario_ids),scenario_ids=allids;else,scenario_ids=double(scenario_ids(:))';assert(all(ismember(scenario_ids,allids)),'Unknown scenario id.');end
replicates=double(replicates(:))';assert(all(ismember(replicates,1:30)),'Bad replicate.');
total=numel(replicates)*numel(scenario_ids);done=0;
fprintf('TFKM-HIER4 BLIND BATCH START: %d unit finalizations | truth sealed.\n',total);
for rep=replicates
    for sid=scenario_ids
        RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1(sid,rep,root);done=done+1;
        if mod(done,10)==0||done==total,fprintf('TFKM-HIER4 BLIND BATCH PROGRESS %d/%d\n',done,total);end
    end
end
OUT=struct('complete',true,'N',total,'method',"REG_TFKM_HD",'truth_opened',false);
end
function [root,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";
for c=cands(:)'
    orig=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');hier=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
    if isfile(fullfile(orig,'study.mat'))&&isfolder(hier),root=c;S=load(fullfile(orig,'study.mat'));break,end
end
assert(strlength(root)>0,'Original fresh study + HIER4 correction-rerun output not found.');
end
