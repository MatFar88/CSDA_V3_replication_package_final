function OUT = RUN_ANOR_TRKM_APPLICATION_V2(mode,worker_id,n_workers,data_file,labels_file,data_banks_file,names_file)
%RUN_ANOR_TRKM_APPLICATION_V2 Top-level ANOR ratios-only TRKM launcher.
% Fitting uses DataECB(:,2:383) only; source column 1 (log assets) excluded.
if nargin<1||strlength(string(mode))==0,mode="RUN";end
if nargin<2||isempty(worker_id),worker_id=1;end
if nargin<3||isempty(n_workers),n_workers=1;end
if nargin<4,data_file="";end;if nargin<5,labels_file="";end;if nargin<6,data_banks_file="";end;if nargin<7,names_file="";end
mode=upper(string(mode));worker_id=double(worker_id);n_workers=double(n_workers);assert(worker_id>=1&&n_workers>=1&&worker_id<=n_workers&&worker_id==floor(worker_id)&&n_workers==floor(n_workers),'Require integer 1 <= worker_id <= n_workers.');
root=fileparts(mfilename('fullpath'));CFG=anor_trkm_config_v2(root);
switch mode
 case "SOURCECHECK"
  paths=anor_trkm_path_lock_v2(root);T=SELFTEST_ANOR_TRKM_PACKAGE_V2(false);[df,lf,bf,nf]=locate_inputs_local(data_file,labels_file,data_banks_file,names_file,false);
  fprintf('\nPrivate-input discovery (no file content printed):\n');fprintf('  data_ECB_ANOR.mat: %s\n',found_local(df));fprintf('  labels_ECB_ANOR.mat: %s\n',found_local(lf));fprintf('  data_banks.mat: %s (optional alignment check)\n',found_local(bf));fprintf('  names_banks.mat: %s (optional length check)\n',found_local(nf));OUT=struct('pass',true,'paths',paths,'selftest',T);
 case "SELFTEST",OUT=SELFTEST_ANOR_TRKM_PACKAGE_V2(false);
 case "FITTEST",OUT=SELFTEST_ANOR_TRKM_PACKAGE_V2(true);
 case "BUILD",[df,lf,bf,nf]=locate_inputs_local(data_file,labels_file,data_banks_file,names_file,true);OUT=BUILD_ANOR363_RATIOS_INPUT_TRKM_V2(df,lf,bf,nf);
 case "PILOT",ensure_input_local(CFG,data_file,labels_file,data_banks_file,names_file);RUN_ANOR_TRKM_CANDIDATES_V2(1,150);T=REPORT_ANOR_TRKM_PROGRESS_V2;OUT=struct('mode',mode,'progress',T);
 case "STATUS",T=REPORT_ANOR_TRKM_PROGRESS_V2;OUT=struct('mode',mode,'progress',T);
 case "RUN"
  ensure_input_local(CFG,data_file,labels_file,data_banks_file,names_file);RUN_ANOR_TRKM_CANDIDATES_V2(worker_id,n_workers);T=REPORT_ANOR_TRKM_PROGRESS_V2;ncomplete=sum(logical(T.complete));OUT=struct('mode',mode,'n_complete',ncomplete,'progress',T);if ncomplete==150&&n_workers==1,OUT.finalization=finalize_local();elseif ncomplete==150,fprintf('\nAll 150 candidates complete. Run FINALIZE once.\n');end
 case "FINALIZE",T=REPORT_ANOR_TRKM_PROGRESS_V2;assert(sum(logical(T.complete))==150&&all(logical(T.valid)),'FINALIZE requires 150/150 valid candidates.');OUT=finalize_local();
 case "DESCRIPTIVE",OUT=RUN_ANOR_TRKM_DESCRIPTIVE_REFIT_V2;
 case "ALL",ensure_input_local(CFG,data_file,labels_file,data_banks_file,names_file);RUN_ANOR_TRKM_CANDIDATES_V2(1,1);T=REPORT_ANOR_TRKM_PROGRESS_V2;assert(sum(logical(T.complete))==150&&all(logical(T.valid)),'Candidate run incomplete.');OUT=finalize_local();
 otherwise,error('Unknown mode %s.',mode);
end
end
function OUT=finalize_local()
P=MERGE_FINALIZE_ANOR_TRKM_V2;OUT=struct('primary',P,'descriptive',struct('run',false,'reason',"NOT_REQUIRED"));Pred=P.prediction;
if ~logical(Pred.final_fit_available)&&logical(Pred.terminal_unique)&&~logical(Pred.action_ambiguous),OUT.descriptive=RUN_ANOR_TRKM_DESCRIPTIVE_REFIT_V2;elseif logical(Pred.final_fit_available),fprintf('Primary selector defines a terminal fit; no secondary descriptive refit required.\n');end
end
function ensure_input_local(CFG,data_file,labels_file,data_banks_file,names_file)
if isfile(CFG.input_file),S=load(CFG.input_file,'protocol_version','pair_seed','external_standardization','source_column_1_status');assert(string(S.protocol_version)==CFG.version&&double(S.pair_seed)==CFG.pair_seed,'Existing input contract mismatch.');assert(string(S.external_standardization)=="NONE_RATIOS_AS_SUPPLIED"&&string(S.source_column_1_status)=="EXCLUDED_LOG_ASSETS",'Existing input violates ratios-only contract.');return,end
[df,lf,bf,nf]=locate_inputs_local(data_file,labels_file,data_banks_file,names_file,true);BUILD_ANOR363_RATIOS_INPUT_TRKM_V2(df,lf,bf,nf);
end
function [df,lf,bf,nf]=locate_inputs_local(data_file,labels_file,data_banks_file,names_file,required)
if nargin<5,required=true;end;docs=string(fullfile(getenv('USERPROFILE'),'Documents'));df=locate_one_local(data_file,"data_ECB_ANOR.mat",docs);lf=locate_one_local(labels_file,"labels_ECB_ANOR.mat",docs);bf=locate_one_local(data_banks_file,"data_banks.mat",docs);nf=locate_one_local(names_file,"names_banks.mat",docs);if required,assert(strlength(df)>0&&isfile(df),'data_ECB_ANOR.mat not found.');assert(strlength(lf)>0&&isfile(lf),'labels_ECB_ANOR.mat not found.');end
end
function f=locate_one_local(arg,name,docs)
if strlength(string(arg))>0,f=string(arg);return,end;c=[string(fullfile(pwd,name));fullfile(fileparts(pwd),name);fullfile(docs,name)];f="";for x=c(:)',if isfile(x),f=x;return,end,end;if isfolder(docs),D=dir(fullfile(docs,'**',name));if ~isempty(D),f=string(fullfile(D(1).folder,D(1).name));end,end
end
function s=found_local(f),if strlength(f)>0,s='FOUND';else,s='not located automatically';end,end
