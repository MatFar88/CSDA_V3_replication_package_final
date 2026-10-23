function REPORT_TFKM_TRKM_FRESH_INVARIANTS_PROGRESS_V1(documents_root)
if nargin<1,documents_root="";end
[~,outroot]=load_study_local(documents_root);outdir=fullfile(outroot,'fresh_invariants_v1');
if ~isfolder(outdir),fprintf('No invariant output directory yet.\n');return,end
D=dir(fullfile(outdir,'worker_*_of_*_invariants.csv'));T=table();
for k=1:numel(D)
    Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');T=[T;Q]; %#ok<AGROW>
end
if isempty(T),fprintf('INVARIANTS: 0/1230\n');return,end
keys=unique(string(T.scenario_id)+"|"+string(T.replicate));fprintf('INVARIANTS: %d/1230 unique units complete\n',numel(keys));
end
function [root,outroot]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";
for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found.');
end
