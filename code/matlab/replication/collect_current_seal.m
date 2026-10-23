function OUT=collect_current_seal(root)
%COLLECT_CURRENT_SEAL Collection-only adapter for the corrected production path.
% Does not run any estimator or selection rule. It avoids requiring the obsolete
% SIGN-only TFKM output just to merge the current TRKM per-unit predictions.
orig=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
M=readtable(fullfile(orig,'UNIT_MANIFEST.csv'));assert(height(M)==1230,'Manifest drift.');
base=fullfile(orig,'blind_finalization_v1','TRKM');
P=table();E=table();R=table();
for k=1:height(M)
    d=fullfile(base,sprintf('s%d',M.scenario_id(k)),sprintf('rep%02d',M.replicate(k)));
    assert(isfile(fullfile(d,'COMPLETE.txt')),'Incomplete TRKM blind result: %s',d);
    p=readtable(fullfile(d,'prediction.csv'),'TextType','string');
    assert(height(p)==1&&p.scenario_id==M.scenario_id(k)&&p.replicate==M.replicate(k),'Wrong TRKM prediction key.');
    P=[P;p]; %#ok<AGROW>
    e=readtable(fullfile(d,'edges.csv'),'TextType','string');if width(e)>3||height(e)>0,E=[E;e];end %#ok<AGROW>
    r=readtable(fullfile(d,'reproduction.csv'),'TextType','string');if width(r)>3||height(r)>0,R=[R;r];end %#ok<AGROW>
end
P=sortrows(P,{'scenario_id','replicate'});check_no_truth(P);check_no_truth(E);check_no_truth(R);
dest=fullfile(orig,'blind_finalization_v1','FROZEN');if ~isfolder(dest),mkdir(dest);end
write_once(P,fullfile(dest,'TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv'));
if ~isempty(E),write_once(E,fullfile(dest,'TRKM_BLIND_EDGE_LOG_FROZEN.csv'));end
if ~isempty(R),write_once(R,fullfile(dest,'TRKM_BLIND_REPRODUCTION_FROZEN.csv'));end
OUT=MERGE_TFKM_HIER4_WITH_TRKM_BLIND_SEAL_V2(root);
% Legacy folder name is retained solely because original evaluators expect it.
% It is a snapshot of THIS replication run, not a new historical holdout seal.
src=fullfile(orig,'blind_finalization_v2','FROZEN');snapshot=fullfile(orig,'FROZEN_PRETRUTH_20261010');
if ~isfolder(snapshot),mkdir(snapshot);end
F=dir(fullfile(src,'*.csv'));
for k=1:numel(F)
    s=fullfile(F(k).folder,F(k).name);d=fullfile(snapshot,F(k).name);
    if isfile(d),assert(strcmp(file_sha256(s),file_sha256(d)),'Replication prediction snapshot drift: %s',F(k).name);else,copyfile(s,d);end
end
meta=struct('sealed_at',char(datetime('now')),'datasets_per_method',1230,'source','regenerated corrected-TFKM and original-TRKM results','legacy_folder_name','FROZEN_PRETRUTH_20261010','historical_holdout',false);
fid=fopen(fullfile(root,'REPLICATION_PREDICTIONS_SEALED.json'),'w');assert(fid>0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(meta));
end
function write_once(T,f)
if isfile(f)
    Q=readtable(f,'TextType','string');
    assert(isequaln(Q,T),'Previously collected table differs: %s',f);
else
    writetable(T,f);
end
end
function check_no_truth(T)
if isempty(T),return,end
names=lower(string(T.Properties.VariableNames));
bad=contains(names,'truth')|contains(names,'posthoc')|contains(names,'outlier')|ismember(names,["r0","g0","c0","alpha0","contamination","rho","delta","p_over_n"]);
assert(~any(bad),'Truth fields found in predictions/operational log.');
end
