function OUT = RUN_TFKM_TRKM_FRESH_BATCH_V1(reps,scenario_ids,methods,documents_root)
%RUN_TFKM_TRKM_FRESH_BATCH_V1 Restartable batch wrapper for fresh fits.
%
% Recommended first workload batch (truth remains closed):
%   RUN_TFKM_TRKM_FRESH_BATCH_V1(1:10,[],"REG_TFKM_HD");
%   RUN_TFKM_TRKM_FRESH_BATCH_V1(1:10,[],"REG_TRKM_HD");
%
% Then extension:
%   RUN_TFKM_TRKM_FRESH_BATCH_V1(11:30,[],"REG_TFKM_HD");
%   RUN_TFKM_TRKM_FRESH_BATCH_V1(11:30,[],"REG_TRKM_HD");
%
% [] scenario_ids means all 41 cells. methods may be one or both strings.
if nargin==1&&strcmpi(string(reps),"SELFTEST")
    RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1("SELFTEST",[],[]);fprintf('SELFTEST PASS: batch wrapper.\n');OUT=struct('selftest_pass',true);return
end
if nargin<1||isempty(reps),reps=1:10;end
if nargin<4,documents_root="";end
if nargin<3||isempty(methods),methods=["REG_TFKM_HD","REG_TRKM_HD"];end
methods=string(methods(:)');
if strlength(string(documents_root))>0,root=string(documents_root);else,root=string(fullfile(getenv('USERPROFILE'),'Documents'));end
outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
if ~isfile(fullfile(outroot,'study.mat')),BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1(root);end
S=load(fullfile(outroot,'study.mat'),'catalog','CFG');
if nargin<2||isempty(scenario_ids),scenario_ids=double(S.catalog.scenario_id)';end
assert(all(ismember(double(scenario_ids),double(S.catalog.scenario_id))),'Unknown scenario IDs.');
assert(all(ismember(double(reps),1:S.CFG.n_replicates)),'Invalid reps.');
N=numel(reps)*numel(scenario_ids)*numel(methods);done=0;
fprintf('\nBATCH START: %d unit-method cubes. Restartable.\n',N);
for method=methods
    for rep=double(reps(:)')
        for sid=double(scenario_ids(:)')
            RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1(sid,rep,method,root);done=done+1;
            fprintf('BATCH PROGRESS %d/%d | %s sid=%d rep=%d\n',done,N,method,sid,rep);
        end
    end
end
OUT=REPORT_TFKM_TRKM_FRESH_PROGRESS_V1(root);
end
