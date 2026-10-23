function REPORT_TFKM_TRKM_FRESH_THEORY_AUDITS_PROGRESS_V1(documents_root)
if nargin<1,documents_root="";end
if strlength(string(documents_root))>0,cands=string(documents_root);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output','fresh_theory_audits_v1');if isfolder(o),outroot=string(o);break,end,end
assert(strlength(outroot)>0,'fresh_theory_audits_v1 not found.');D=dir(fullfile(outroot,'worker_*_of_*_theory_audits.csv'));T=table();for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name) ,'TextType','string');T=[T;Q];end %#ok<AGROW>
if isempty(T),fprintf('Progress: 0/2460 method-units.\n');return,end
keys=string(T.method)+"|"+string(T.scenario_id)+"|"+string(T.replicate);fprintf('Progress: %d/2460 unique method-units.\n',numel(unique(keys)));
if ismember('method',T.Properties.VariableNames),disp(groupsummary(T,'method','numel','scenario_id'));end
end
