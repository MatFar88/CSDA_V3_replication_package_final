function OUT = RUN_ECB_TRKM_APPLICATION_V1(mode,worker_id,n_workers,ecb_mat_file,names364_file)
%RUN_ECB_TRKM_APPLICATION_V1 Top-level launcher for the ECB TRKM application.
%
% Recommended first use:
%   RUN_ECB_TRKM_APPLICATION_V1("SOURCECHECK");
%   RUN_ECB_TRKM_APPLICATION_V1("PILOT");
%
% Four independent MATLAB sessions:
%   RUN_ECB_TRKM_APPLICATION_V1("RUN",1,4);
%   RUN_ECB_TRKM_APPLICATION_V1("RUN",2,4);
%   RUN_ECB_TRKM_APPLICATION_V1("RUN",3,4);
%   RUN_ECB_TRKM_APPLICATION_V1("RUN",4,4);
%   RUN_ECB_TRKM_APPLICATION_V1("FINALIZE");
%
% Modes: SOURCECHECK, SELFTEST, FITTEST, BUILD, PILOT, STATUS, RUN,
% FINALIZE, DESCRIPTIVE, ALL.
%
% Locked application contract:
% n=363, p=361, raw row 93 removed as duplicate of raw row 99, original
% log-level scales, no external standardization, paper 150-state grid,
% 100 TRKM starts allocated 40/30/15/15, RV-CV start arbitration, frozen
% TRKM action/finalizer, and application pair seed 20261013.

if nargin<1 || strlength(string(mode))==0,mode="RUN";end
if nargin<2 || isempty(worker_id),worker_id=1;end
if nargin<3 || isempty(n_workers),n_workers=1;end
if nargin<4,ecb_mat_file="";end
if nargin<5,names364_file="";end
mode=upper(string(mode));worker_id=double(worker_id);n_workers=double(n_workers);
assert(worker_id>=1&&n_workers>=1&&worker_id<=n_workers&&worker_id==floor(worker_id)&&n_workers==floor(n_workers),'Require integer 1 <= worker_id <= n_workers.');
root=fileparts(mfilename('fullpath'));CFG=ecb_trkm_config_v1(root);

switch mode
    case "SOURCECHECK"
        paths=ecb_trkm_path_lock_v1(root);T=SELFTEST_ECB_TRKM_PACKAGE_V1(false);[ef,nf]=locate_private_inputs_local(ecb_mat_file,names364_file,false);
        fprintf('\nPrivate-input discovery (no file content printed):\n');
        if strlength(ef)>0,fprintf('  ECB_cov_input.mat: FOUND\n');else,fprintf('  ECB_cov_input.mat: not located automatically\n');end
        if strlength(nf)>0,fprintf('  names_bank_364.mat: FOUND (optional)\n');else,fprintf('  names_bank_364.mat: not located automatically (optional)\n');end
        OUT=struct('pass',true,'paths',paths,'selftest',T,'ecb_input_found',strlength(ef)>0,'names_input_found',strlength(nf)>0);
    case "SELFTEST"
        OUT=SELFTEST_ECB_TRKM_PACKAGE_V1(false);
    case "FITTEST"
        OUT=SELFTEST_ECB_TRKM_PACKAGE_V1(true);
    case "BUILD"
        [ef,nf]=locate_private_inputs_local(ecb_mat_file,names364_file,true);OUT=BUILD_ECB363_INPUT_TRKM_V1(ef,nf);
    case "PILOT"
        ensure_private_input_local(CFG,ecb_mat_file,names364_file);RUN_ECB_TRKM_CANDIDATES_V1(1,150);T=REPORT_ECB_TRKM_PROGRESS_V1;OUT=struct('mode',mode,'progress',T);
    case "STATUS"
        T=REPORT_ECB_TRKM_PROGRESS_V1;OUT=struct('mode',mode,'progress',T);
    case "RUN"
        ensure_private_input_local(CFG,ecb_mat_file,names364_file);RUN_ECB_TRKM_CANDIDATES_V1(worker_id,n_workers);T=REPORT_ECB_TRKM_PROGRESS_V1;ncomplete=sum(logical(T.complete));
        OUT=struct('mode',mode,'worker_id',worker_id,'n_workers',n_workers,'n_complete',ncomplete,'progress',T);
        if ncomplete==150 && n_workers==1
            OUT.finalization=finalize_and_optional_descriptive_local();
        elseif ncomplete==150
            fprintf('\nAll 150 candidates are complete. Run FINALIZE once from one MATLAB session.\n');
        end
    case "FINALIZE"
        T=REPORT_ECB_TRKM_PROGRESS_V1;assert(sum(logical(T.complete))==150,'FINALIZE requires 150/150 complete candidates.');assert(all(logical(T.valid)),'At least one completed candidate is invalid; investigate before finalization.');OUT=finalize_and_optional_descriptive_local();
    case "DESCRIPTIVE"
        OUT=RUN_ECB_TRKM_DESCRIPTIVE_REFIT_V1;
    case "ALL"
        ensure_private_input_local(CFG,ecb_mat_file,names364_file);RUN_ECB_TRKM_CANDIDATES_V1(1,1);T=REPORT_ECB_TRKM_PROGRESS_V1;assert(sum(logical(T.complete))==150&&all(logical(T.valid)),'Serial candidate run did not finish with 150 valid candidates.');OUT=finalize_and_optional_descriptive_local();
    otherwise
        error('Unknown mode %s. Use SOURCECHECK, SELFTEST, FITTEST, BUILD, PILOT, STATUS, RUN, FINALIZE, DESCRIPTIVE, or ALL.',mode);
end
end

function OUT=finalize_and_optional_descriptive_local()
P=MERGE_FINALIZE_ECB_TRKM_V1;OUT=struct('primary',P,'descriptive',struct('run',false,'reason',"NOT_REQUIRED"));Pred=P.prediction;
if ~logical(Pred.final_fit_available)&&logical(Pred.terminal_unique)&&~logical(Pred.action_ambiguous)
    OUT.descriptive=RUN_ECB_TRKM_DESCRIPTIVE_REFIT_V1;
elseif logical(Pred.final_fit_available)
    fprintf('Primary selector defines a terminal fit; no secondary descriptive refit is required.\n');
else
    fprintf('Secondary descriptive refit not run because the primary structural output is ambiguous.\n');OUT.descriptive=struct('run',false,'reason',"PRIMARY_AMBIGUITY");
end
end

function ensure_private_input_local(CFG,ecb_mat_file,names364_file)
if isfile(CFG.input_file)
    S=load(CFG.input_file,'protocol_version','pair_seed','external_standardization');
    assert(isfield(S,'protocol_version')&&string(S.protocol_version)==CFG.version,'Existing internal input belongs to a different protocol version.');
    assert(isfield(S,'pair_seed')&&double(S.pair_seed)==CFG.pair_seed,'Existing internal input has a different pair seed.');
    assert(isfield(S,'external_standardization')&&string(S.external_standardization)=="NONE_LOG_LEVELS_RETAINED",'Existing internal input violates scale contract.');return
end
[ef,nf]=locate_private_inputs_local(ecb_mat_file,names364_file,true);BUILD_ECB363_INPUT_TRKM_V1(ef,nf);
end

function [ecb_file,names_file]=locate_private_inputs_local(ecb_mat_file,names364_file,require_ecb)
if nargin<3,require_ecb=true;end
docs=string(fullfile(getenv('USERPROFILE'),'Documents'));
if strlength(string(ecb_mat_file))>0,ecb_file=string(ecb_mat_file);else
    direct=[string(fullfile(pwd,'ECB_cov_input.mat'));fullfile(docs,'R2_ECB_APPLICATION_V12_CEIL','ECB_cov_input.mat');fullfile(docs,'R2_ECB_APPLICATION_V12_CEIL_package','ECB_cov_input.mat');fullfile(docs,'R2_ECB_APPLICATION_V12_CEIL_package(1)','ECB_cov_input.mat')];
    ecb_file=first_existing_local(direct);if strlength(ecb_file)==0&&isfolder(docs),D=dir(fullfile(docs,'**','ECB_cov_input.mat'));if ~isempty(D),ecb_file=string(fullfile(D(1).folder,D(1).name));end,end
end
if require_ecb,assert(strlength(ecb_file)>0&&isfile(ecb_file),'ECB_cov_input.mat could not be located; pass its full path as argument 4.');end
if strlength(string(names364_file))>0,names_file=string(names364_file);else
    direct=[string(fullfile(pwd,'names_bank_364.mat'));fullfile(docs,'names_bank_364.mat')];names_file=first_existing_local(direct);if strlength(names_file)==0&&isfolder(docs),D=dir(fullfile(docs,'**','names_bank_364.mat'));if ~isempty(D),names_file=string(fullfile(D(1).folder,D(1).name));end,end
end
if strlength(names_file)>0&&~isfile(names_file),names_file="";end
end
function f=first_existing_local(cands),f="";for c=string(cands(:))',if isfile(c),f=c;return,end,end,end
