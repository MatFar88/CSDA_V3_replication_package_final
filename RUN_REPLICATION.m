function OUT = RUN_REPLICATION(mode,varargin)
%RUN_REPLICATION Entry point for V15 saved-result reporting and explicit reruns.
% CHECK/REPORT never fit an estimator. FULL requires Confirm=true.
%
% RUN_REPLICATION("SETUP", "Confirm", true)     % Windows x64; permits setup downloads
% RUN_REPLICATION("CHECK")
% RUN_REPLICATION("REPORT", "Python", "C:\path\to\python.exe")
% RUN_REPLICATION("SELFTEST")                    % small source-contract tests
% RUN_REPLICATION("FULL", "Workspace", "D:\CSDA_REPLAY", ...
%                 "Stage", "PREPARE", "Confirm", true)
% RUN_REPLICATION("FULL", "Workspace", "D:\CSDA_REPLAY", ...
%                 "Stage", "FIT", "Worker", 1, "Workers", 6, "Confirm", true)
%
% See documentation/FULL_MATLAB.md before starting expensive calculations.
% This wrapper was statically inspected but has not been run in MATLAB in the
% package-building environment. Original numerical routines are unmodified.
if nargin<1,mode="CHECK";end
ip=inputParser;
addParameter(ip,'Python','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Output','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Workspace','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Experiment','FRESH',@(x)ischar(x)||isstring(x));
addParameter(ip,'Session','ALL',@(x)ischar(x)||isstring(x));
addParameter(ip,'Stage','ALL',@(x)ischar(x)||isstring(x));
addParameter(ip,'Worker',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==floor(x));
addParameter(ip,'Workers',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==floor(x));
addParameter(ip,'MaxUnits',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&(isinf(x)||x==floor(x)));
addParameter(ip,'Confirm',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
mode=upper(string(mode));stage=upper(string(o.Stage));
assert(o.Worker<=o.Workers,'Worker must not exceed Workers.');
root=fileparts(mfilename('fullpath'));
if mode=="SETUP"
    assert(ispc,'Automatic setup is provided for Windows x64. See setup_python.sh for other platforms.');
    assert(o.Confirm,'SETUP downloads reporting dependencies if needed. Use Confirm=true to permit these downloads.');
    oldpath=path;restore=onCleanup(@()path(oldpath)); %#ok<NASGU>
    addpath(fullfile(root,'code','matlab','replication'),'-begin');
    verify_package_native(root);
    command=sprintf('cmd /d /s /c ""%s" --yes"',fullfile(root,'SETUP_PYTHON_WINDOWS.cmd'));
    [status,txt]=system(command);fprintf('%s',txt);
    assert(status==0,'SETUP did not finish. See output/setup_logs; no study was fitted.');
    OUT=struct('mode',mode,'setup_exit_code',status);return
end
if strlength(string(o.Python))==0
    if ispc,candidate=fullfile(root,'.venv','Scripts','python.exe');else,candidate=fullfile(root,'.venv','bin','python');end
    if isfile(candidate)
        o.Python=candidate;
    elseif mode=="CHECK"||mode=="REPORT"
        error('CSDA:ReportingSetupMissing', ...
            ['The reporting environment is not configured. On Windows run ' ...
             'RUN_REPLICATION("SETUP","Confirm",true) once; Python need not be installed beforehand. ' ...
             'Alternatively specify a configured interpreter with "Python".']);
    else
        o.Python='python'; % SELFTEST and native FULL do not use this value.
    end
end
oldpath=path;restore=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(fullfile(root,'code','matlab','replication'),'-begin');
if mode=="CHECK"||mode=="REPORT"
    preflight=sprintf('"%s" -I "%s"',char(o.Python),fullfile(root,'tools','check_python_environment.py'));
    [ready,details]=system(preflight);
    if ready~=0
        fprintf('%s',details);
        error('CSDA:ReportingSetupIncomplete', ...
            ['Reporting dependencies are missing or do not match requirements.txt. ' ...
             'On Windows run RUN_REPLICATION("SETUP","Confirm",true). No fit has started.']);
    end
    command=sprintf('"%s" "%s" %s',char(o.Python),fullfile(root,'run_replication.py'),char(mode));
    if strlength(string(o.Output))>0
        command=[command sprintf(' --output "%s"',char(o.Output))]; %#ok<AGROW>
    end
    [status,txt]=system(command);fprintf('%s',txt);
    assert(status==0,'Python %s failed (exit %d). Check the Python path and requirements.txt.',mode,status);
    OUT=struct('mode',mode,'python_exit_code',status,'matlab_version',version,'matlab_release',version('-release'));
    return
end
verify_package_native(root);
stages=dir(fullfile(root,'code','matlab','stages'));
for k=1:numel(stages)
    if stages(k).isdir&&~startsWith(stages(k).name,'.')
        addpath(fullfile(stages(k).folder,stages(k).name),'-begin');
    end
end
if mode=="SELFTEST"
    BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1("SELFTEST");
    RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1("SELFTEST");
    RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1("SELFTEST");
    RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1("SELFTEST");
    RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1("SELFTEST");
    RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1("SELFTEST",[],[]);
    RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1("SELFTEST",[],[]);
    RUN_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1("SELFTEST",[],[]);
    fprintf('Original MATLAB selftests completed; no full study was fitted.\n');
    OUT=struct('mode',mode,'matlab_release',version('-release'),'status','PASS');
    return
end
assert(mode=="FULL",'Mode must be SETUP, CHECK, REPORT, SELFTEST or FULL.');
assert(upper(string(o.Experiment))=="FRESH",'Experiment must be FRESH. FINAL32 reporting uses saved results; see documentation/FIXED32.md for its separate source archive.');
assert(o.Confirm,'FULL requires Confirm=true. REPORT does not refit; FULL may require substantial computation.');
assert(strlength(string(o.Workspace))>0,'FULL requires an explicit new Workspace.');
if stage~="PREPARE" && stage~="SEAL"
    assert(exist('kmeans','file')==2 && exist('pdist','file')==2,'The original fitting/generator code requires kmeans and pdist (Statistics and Machine Learning Toolbox).');
end
allowed=["ALL","PREPARE","DGP_SMOKE","FIT_SMOKE","FIT","SELECT","SEAL","EVALUATE","DIAGNOSTICS","MERGE_DIAGNOSTICS"];
assert(any(stage==allowed),'Unknown Stage. See documentation/FULL_MATLAB.md.');
if stage=="ALL",assert(o.Workers==1,'ALL is serial. Use separate stage invocations when Workers>1.');end
if any(stage==["SEAL","EVALUATE","MERGE_DIAGNOSTICS"])
    assert(o.Worker==1,'Run this merge/evaluation stage once, with Worker=1.');
end
if isfinite(o.MaxUnits)
    assert(any(stage==["FIT","SELECT"]),'MaxUnits is supported only for FIT or SELECT pilot subsets.');
end
workspace=prepare_workspace(root,char(o.Workspace));
lockfile=fullfile(workspace,sprintf('.REPLICATION_%s_%d_OF_%d.lock',char(stage),o.Worker,o.Workers));
lf=java.io.File(lockfile);
assert(lf.createNewFile(),'A matching stage/worker lock exists. Verify that no session is active before removing %s.',lockfile);
release=onCleanup(@()delete_if_exists(lockfile)); %#ok<NASGU>
logdir=fullfile(workspace,'replication_logs');if ~isfolder(logdir),mkdir(logdir);end
stamp=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
logfile=fullfile(logdir,sprintf('%s_%s_worker%d.txt',stamp,char(stage),o.Worker));
diary(logfile);logcleanup=onCleanup(@()diary('off')); %#ok<NASGU>
E=struct('package','CSDA_REPLICATION_PACKAGE_V1','started_at',char(datetime('now')), ...
    'mode','FULL','stage',char(stage),'workspace',workspace,'worker',o.Worker,'workers',o.Workers, ...
    'matlab_version',version,'matlab_release',version('-release'),'computer',computer, ...
    'toolboxes',{ver},'historical_holdout',false,'meaning','Replication of already published/opened design; not new holdout evidence');
write_json(fullfile(logdir,[stamp '_ENVIRONMENT.json']),E);
fprintf('\nREPLICATION workspace: %s\nStage: %s | worker %d/%d\n',workspace,stage,o.Worker,o.Workers);
fprintf('Reference results will not be overwritten. Original numerical rules are unchanged.\n');
sequence=stage;
if stage=="ALL",sequence=["DGP_SMOKE","FIT","SELECT","SEAL","EVALUATE","DIAGNOSTICS","MERGE_DIAGNOSTICS"];end
for step=sequence
    switch step
        case "PREPARE"
            fprintf('Prepared original study configuration and isolated dependencies. No fit started.\n');
        case "DGP_SMOKE"
            RUN_TFKM_TRKM_FRESH_DGP_SMOKE_V1(workspace);
        case "FIT_SMOKE"
            smoke_fixed_candidate(workspace);
        case {"FIT","SELECT"}
            M=readtable(fullfile(workspace,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output','UNIT_MANIFEST.csv'));
            idx=find(mod((1:height(M))'-1,o.Workers)+1==o.Worker);
            if isfinite(o.MaxUnits),idx=idx(1:min(numel(idx),o.MaxUnits));end
            register_worker_layout(workspace,step,o.Workers);
            for j=idx(:)'
                sid=double(M.scenario_id(j));rep=double(M.replicate(j));
                fprintf('\n%s dataset %d/1230: sid=%d rep=%d\n',step,j,sid,rep);
                if step=="FIT"
                    % Deliberately NOT the superseded SIGN-only TFKM branch.
                    RUN_TFKM_FRESH_HIER4_RERUN_UNIT_V1(sid,rep,workspace);
                    RUN_TFKM_TRKM_FRESH_UNIT_FIT_V1(sid,rep,"REG_TRKM_HD",workspace);
                else
                    RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1(sid,rep,workspace);
                    RUN_TFKM_TRKM_FRESH_BLIND_UNIT_V1(sid,rep,"REG_TRKM_HD",workspace);
                end
            end
        case "SEAL"
            collect_current_seal(workspace);
        case "EVALUATE"
            assert(isfile(fullfile(workspace,'REPLICATION_PREDICTIONS_SEALED.json')),'Run SEAL after every dataset has completed.');
            RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1(workspace);
            RUN_TFKM_TRKM_FRESH_POSTHOC_DIAGNOSTIC_V1(workspace);
            RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1(workspace);
        case "DIAGNOSTICS"
            register_worker_layout(workspace,step,o.Workers);
            RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1(o.Worker,o.Workers,workspace);
            RUN_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1(o.Worker,o.Workers,workspace);
            RUN_TFKM_TRKM_FRESH_INVARIANTS_V1(o.Worker,o.Workers,workspace);
            RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1(o.Worker,o.Workers,workspace);
        case "MERGE_DIAGNOSTICS"
            require_all_diagnostic_workers(workspace);
            MERGE_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1(workspace);
            MERGE_TFKM_TRKM_FRESH_THEORY_AUDITS_V1(workspace);
            MERGE_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1(workspace);
            MERGE_TFKM_TRKM_FRESH_INVARIANTS_V1(workspace);
    end
end
OUT=struct('mode',mode,'stage',stage,'workspace',string(workspace),'logfile',string(logfile),'completed_at',char(datetime('now')));
write_json(fullfile(logdir,[stamp '_COMPLETED.json']),OUT);
fprintf('\nRequested stage finished. For numerical comparison use REPORT with --fresh-root pointing to the regenerated original-output tree.\n');
end

function workspace=prepare_workspace(root,workspace)
jroot=java.io.File(root);root=char(jroot.getCanonicalPath());
jworkspace=java.io.File(workspace);workspace=char(jworkspace.getCanonicalPath());
assert(~strcmpi(root,workspace)&&~startsWith(root,[workspace filesep],'IgnoreCase',true),'Workspace must not be the package or its ancestor.');
if startsWith(workspace,[root filesep],'IgnoreCase',true)
    assert(startsWith(workspace,[fullfile(root,'output') filesep],'IgnoreCase',true),'Inside the package, use only output/<new_workspace>.');
end
marker=fullfile(workspace,'REPLICATION_WORKSPACE.json');
if isfolder(workspace)&&~isfile(marker)
    q=dir(workspace);q=q(~ismember({q.name},{'.','..'}));
    assert(isempty(q),'Workspace is not empty and is not a replication workspace. Original result folders must never be used as output.');
end
if ~isfolder(workspace),mkdir(workspace);end
if ~isfile(marker)
    write_json(marker,struct('package','CSDA_REPLICATION_PACKAGE_V1','created_at',char(datetime('now')),'synthetic_study','V15 41-cell replication'));
else
    x=jsondecode(fileread(marker));assert(strcmp(x.package,'CSDA_REPLICATION_PACKAGE_V1'),'Workspace belongs to a different package.');
end
dest=fullfile(workspace,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if ~isfolder(dest),mkdir(dest);end
for name=["study.mat","STUDY_CATALOG.csv","UNIT_MANIFEST.csv"]
    copy_immutable(fullfile(root,'results','reference','fresh',name),fullfile(dest,name));
end
% Stage native path layouts expected by original source; never modify source.
base=fullfile(root,'code','matlab','dependencies');allfiles=dir(fullfile(base,'**','*'));
for i=1:numel(allfiles)
    if allfiles(i).isdir,continue,end
    src=fullfile(allfiles(i).folder,allfiles(i).name);rel=src(numel(base)+2:end);
    copy_immutable(src,fullfile(workspace,rel));
end
end
function copy_immutable(src,dst)
if isfile(dst)
    assert(strcmp(file_sha256(src),file_sha256(dst)),'Existing staged input differs from the package: %s',dst);
else
    parent=fileparts(dst);if ~isfolder(parent),mkdir(parent);end
    copyfile(src,dst);
end
end
function register_worker_layout(root,stage,n)
f=fullfile(root,sprintf('REPLICATION_%s_WORKERS.json',char(stage)));
if isfile(f)
    q=jsondecode(fileread(f));assert(q.workers==n,'Keep the same Workers value for resumed %s stages.',stage);
else
    write_json(f,struct('stage',char(stage),'workers',n));
end
end
function require_all_diagnostic_workers(root)
f=fullfile(root,'REPLICATION_DIAGNOSTICS_WORKERS.json');assert(isfile(f),'No diagnostic worker configuration.');
q=jsondecode(fileread(f));n=q.workers;
out=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
% Detailed mergers additionally check row counts/keys. Here forbid unnoticed
% mixing of old worker layouts before invoking them.
folders={fullfile(out,'fresh_theory_audits_v1'),fullfile(out,'fresh_theory_closure_audits_v1'),fullfile(out,'fresh_invariants_v1'),fullfile(out,'fresh_v2_diagnostic_suite_v1','weak_true_state_pool_audit_v1')};
for j=1:numel(folders)
    for k=1:n
        a=dir(fullfile(folders{j},sprintf('worker_%02d_of_%02d_*.csv',k,n)));
        assert(~isempty(a),'Missing diagnostic worker %d/%d in %s.',k,n,folders{j});
    end
end
end
function delete_if_exists(f)
if isfile(f),delete(f);end
end
function write_json(f,s)
tmp=[tempname(fileparts(f)) '.json'];
fid=fopen(tmp,'w');assert(fid>0,'Cannot open %s',tmp);
fprintf(fid,'%s\n',jsonencode(s));fclose(fid);
[ok,msg]=movefile(tmp,f,'f');assert(ok,'Cannot install JSON metadata: %s',msg);
end
