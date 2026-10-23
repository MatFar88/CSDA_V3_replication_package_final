function OUT = RUN_TFKM_FRESH_HIER4_RERUN_BATCH_V1(reps,scenario_ids,documents_root)
% Restartable corrected TFKM HIER4-only fresh rerun. Truth sealed.
if nargin==1 && strcmpi(string(reps),"SELFTEST")
    RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1("SELFTEST",[]);fprintf('SELFTEST PASS: HIER4 rerun batch wrapper.\n');OUT=struct('selftest_pass',true);return
end
if nargin<1||isempty(reps),reps=1:30;end
if nargin<3,documents_root="";end
if strlength(string(documents_root))>0,root=string(documents_root);else,root=string(fullfile(getenv('USERPROFILE'),'Documents'));end
src=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');assert(isfile(fullfile(src,'study.mat')),'Original fresh study not found.');
S=load(fullfile(src,'study.mat'),'catalog','CFG');
if nargin<2||isempty(scenario_ids),scenario_ids=double(S.catalog.scenario_id)';end
assert(all(ismember(double(scenario_ids),double(S.catalog.scenario_id))),'Unknown scenario IDs.');assert(all(ismember(double(reps),1:S.CFG.n_replicates)),'Invalid reps.');
N=numel(reps)*numel(scenario_ids);done=0;
fprintf('\nHIER4 RERUN BATCH START: %d TFKM unit cubes | truth sealed.\n',N);
for rep=double(reps(:)')
    for sid=double(scenario_ids(:)')
        RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1(sid,rep,root);done=done+1;
        if mod(done,5)==0||done==N,fprintf('HIER4 RERUN BATCH PROGRESS %d/%d | sid=%d rep=%d\n',done,N,sid,rep);end
    end
end
OUT=REPORT_TFKM_FRESH_HIER4_RERUN_PROGRESS_V1(root);
end
