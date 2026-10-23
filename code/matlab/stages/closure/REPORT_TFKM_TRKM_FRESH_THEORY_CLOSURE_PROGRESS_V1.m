function REPORT_TFKM_TRKM_FRESH_THEORY_CLOSURE_PROGRESS_V1(documents_root)
%REPORT_TFKM_TRKM_FRESH_THEORY_CLOSURE_PROGRESS_V1 Compatibility-safe progress report.
if nargin<1,documents_root="";end
[~,outroot]=load_study_local(documents_root);
outdir=fullfile(outroot,'fresh_theory_closure_audits_v1');
if ~isfolder(outdir),fprintf('No closure-audit output directory yet: %s\n',outdir);return,end
D=dir(fullfile(outdir,'worker_*_of_*_grank.csv'));GR=table();
for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');GR=[GR;Q];end %#ok<AGROW>
if ~isempty(GR)
    keys=string(GR.scenario_id)+"|"+string(GR.replicate);nGR=numel(unique(keys));
else,nGR=0;end
D=dir(fullfile(outdir,'worker_*_of_*_local_return.csv'));LR=table();
for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');LR=[LR;Q];end %#ok<AGROW>
if ~isempty(LR)
    keys=string(LR.method)+"|"+string(LR.scenario_id)+"|"+string(LR.replicate)+"|"+string(LR.perturb_type)+"|"+string(LR.perturb_count);
    [~,ia]=unique(keys,'stable');LR=LR(ia,:);keys=keys(ia);nLR=numel(keys);
    ukeys=string(LR.scenario_id)+"|"+string(LR.replicate);u=unique(ukeys);complete=0;
    for j=1:numel(u),if sum(ukeys==u(j))==20,complete=complete+1;end,end
else,nLR=0;complete=0;end
fprintf('Structural g/r units: %d/1230 unique datasets.\n',nGR);
fprintf('Local-return rows:    %d/24600 unique rows.\n',nLR);
fprintf('Fully complete units: %d/1230 datasets.\n',complete);
end

function [root,outroot]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found. Pass Documents explicitly.');
end
