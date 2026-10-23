function O = AUDIT_TFKM_FRESH_HIER4_BLIND_FINALIZATION_V1(documents_root)
%AUDIT_TFKM_FRESH_HIER4_BLIND_FINALIZATION_V1 Provenance/completeness only; no truth.
if nargin<1,documents_root="";end
[root,S]=load_study_local(documents_root);N=height(S.manifest);assert(N==1230,'Manifest size drift.');
base=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output','blind_finalization_hier4_v1','TFKM');
D=dir(fullfile(base,'s*','rep*','COMPLETE.txt'));
nComplete=numel(D);nPred=0;nCube=0;nRows=0;nRepro=0;nReproFail=0;maxIdErr=0;badTruth=0;
keys=zeros(0,2);
for k=1:numel(D)
    udir=D(k).folder;pf=fullfile(udir,'prediction.csv');rf=fullfile(udir,'reproduction.csv');assert(isfile(pf)&&isfile(rf),'Missing unit export: %s',udir);
    P=readtable(pf,'TextType','string');assert(height(P)==1,'Prediction row count !=1: %s',pf);nPred=nPred+1;
    keys(end+1,:)=[double(P.scenario_id(1)),double(P.replicate(1))]; %#ok<AGROW>
    badTruth=badTruth+has_truth_columns_local(P);
    R=readtable(rf,'TextType','string');
    if width(R)>3
        badTruth=badTruth+has_truth_columns_local(R);nRepro=nRepro+height(R);
        if ismember('reproduction_pass',R.Properties.VariableNames)
            pass=parse_logical_local(R.reproduction_pass);nReproFail=nReproFail+sum(~pass);
        end
    end
end
% Re-audit all corrected HIER4 input cubes so the finalization seal cannot
% accidentally be attached to legacy SIGN_ALL_STABLE inputs.
ids=double(S.catalog.scenario_id(:));
for rep=1:30
    for j=1:numel(ids)
        f=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output','fits','TFKM',sprintf('s%d',ids(j)),sprintf('rep%02d',rep),'cube.mat');
        c=fullfile(fileparts(f),'COMPLETE.txt');assert(isfile(f)&&isfile(c),'Corrected cube missing.');
        Q=load(f,'Candidates','unit_info');C=Q.Candidates;nCube=nCube+1;nRows=nRows+height(C);
        req=["hier_cross_ari_score","hier_cross_j_score","hier_cross_H_score","selector_score"];assert(all(ismember(req,string(C.Properties.VariableNames))),'HIER columns missing.');
        A=double(C.hier_cross_ari_score);J=double(C.hier_cross_j_score);H=double(C.hier_cross_H_score);Sc=double(C.selector_score);
        assert(all(isfinite(A))&&all(isfinite(J))&&all(isfinite(H))&&all(isfinite(Sc)),'Nonfinite HIER evidence.');
        e=max(abs(Sc-(A+.2*J)));maxIdErr=max(maxIdErr,e);
        if isstruct(Q.unit_info)&&isfield(Q.unit_info,'hier_wrapper_path'),assert(contains(string(Q.unit_info.hier_wrapper_path),[filesep 'vendor' filesep 'hier4' filesep],'IgnoreCase',true),'Bad HIER wrapper provenance.');end
    end
end
uniqueKeys=size(unique(keys,'rows'),1);
pass=(nComplete==N)&&(nPred==N)&&(uniqueKeys==N)&&(nCube==N)&&(nRows==184500)&&(nReproFail==0)&&(maxIdErr<=1e-10)&&(badTruth==0);
fprintf('\n============================================================\n');
fprintf('TFKM HIER4 BLIND FINALIZATION AUDIT - TRUTH SEALED\n');
fprintf('============================================================\n');
fprintf('blind COMPLETE markers      : %d / %d\n',nComplete,N);
fprintf('blind prediction rows       : %d / %d\n',nPred,N);
fprintf('unique prediction keys      : %d / %d\n',uniqueKeys,N);
fprintf('corrected HIER4 cubes       : %d / %d\n',nCube,N);
fprintf('corrected HIER4 rows        : %d / 184500\n',nRows);
fprintf('max |selector-(A+.2J)|      : %.3g\n',maxIdErr);
fprintf('reproduction rows checked   : %d\n',nRepro);
fprintf('reproduction failures       : %d\n',nReproFail);
fprintf('truth-column violations     : %d\n',badTruth);
if pass
    fprintf('STATUS: HIER4_BLIND_FINALIZATION_CONTRACT_CONFIRMED\n');
else
    fprintf('STATUS: STOP_BEFORE_COMBINED_SEAL\n');
end
fprintf('============================================================\n');
O=struct('pass',pass,'complete',nComplete,'predictions',nPred,'unique_keys',uniqueKeys,'hier4_cubes',nCube,'hier4_rows',nRows,'max_identity_error',maxIdErr,'reproduction_rows',nRepro,'reproduction_failures',nReproFail,'truth_column_violations',badTruth,'truth_opened',false);
end
function b=has_truth_columns_local(T)
if isempty(T),b=0;return,end;n=lower(string(T.Properties.VariableNames));bad=contains(n,"truth")|contains(n,"posthoc")|contains(n,"outlier")|ismember(n,["r0","g0","c0","alpha0","contamination","rho","delta","p_over_n"]);b=sum(bad);
end
function x=parse_logical_local(v)
if islogical(v),x=v;elseif isnumeric(v),x=logical(v);else,s=lower(strtrim(string(v)));x=ismember(s,["1","true","yes"]);end
end
function [root,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";
for c=cands(:)'
    orig=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');hier=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
    if isfile(fullfile(orig,'study.mat'))&&isfolder(hier),root=c;S=load(fullfile(orig,'study.mat'));break,end
end
assert(strlength(root)>0,'Original fresh study + HIER4 correction-rerun output not found.');
end
