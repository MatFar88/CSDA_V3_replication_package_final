function OUT = AUDIT_TFKM_FRESH_HIER4_RERUN_V1(documents_root)
% Pre-truth provenance/completeness audit for corrected HIER4 rerun.
if nargin<1||strlength(string(documents_root))==0,root=string(fullfile(getenv('USERPROFILE'),'Documents'));else,root=string(documents_root);end
src=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');assert(isfile(fullfile(src,'study.mat')),'Original fresh study missing.');S=load(fullfile(src,'study.mat'),'catalog','CFG','options');
outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
Ncube=0;Nrows=0;NfiniteA=0;NfiniteJ=0;NfiniteBoth=0;Nidentity=0;maxerr=0;Ncomplete=0;bad=strings(0,1);labels=strings(0,1);
for rep=1:S.CFG.n_replicates
 for sid=double(S.catalog.scenario_id)'
  d=fullfile(outroot,'fits','TFKM',sprintf('s%d',sid),sprintf('rep%02d',rep));f=fullfile(d,'cube.mat');
  if ~isfile(f),continue,end
  Ncube=Ncube+1;Q=load(f,'Candidates');C=Q.Candidates;Nrows=Nrows+height(C);
  req={'hier_cross_ari_score','hier_cross_j_score','hier_cross_H_score','selector_score'};
  if ~all(ismember(req,C.Properties.VariableNames)),bad(end+1)=string(f);continue,end
  A=double(C.hier_cross_ari_score);J=double(C.hier_cross_j_score);H=double(C.hier_cross_H_score);Z=double(C.selector_score);
  fa=isfinite(A);fj=isfinite(J);fb=fa&fj&isfinite(H)&isfinite(Z);NfiniteA=NfiniteA+sum(fa);NfiniteJ=NfiniteJ+sum(fj);NfiniteBoth=NfiniteBoth+sum(fb);
  if any(fb),e=abs(Z(fb)-(A(fb)+double(S.options.trim_stability_weight)*J(fb)));Nidentity=Nidentity+sum(fb);maxerr=max(maxerr,max(e));end
  if ismember('internal_selector',C.Properties.VariableNames),labels=[labels;string(C.internal_selector)];end %#ok<AGROW>
  if isfile(fullfile(d,'COMPLETE.txt')),Ncomplete=Ncomplete+1;end
 end
end
fprintf('\n============================================================\nTFKM HIER4 RERUN PROVENANCE AUDIT - TRUTH SEALED\n============================================================\n');
fprintf('cubes found                 : %d / 1230\n',Ncube);fprintf('complete markers            : %d / 1230\n',Ncomplete);fprintf('candidate rows              : %d / 184500\n',Nrows);
fprintf('finite HIER A rows          : %d\n',NfiniteA);fprintf('finite HIER J rows          : %d\n',NfiniteJ);fprintf('finite HIER A+J+H+S rows    : %d\n',NfiniteBoth);fprintf('identity rows checked       : %d\n',Nidentity);fprintf('max |selector-(A+.2J)|      : %.3g\n',maxerr);
if ~isempty(labels),u=unique(labels);fprintf('internal selector labels:\n');for k=1:numel(u),fprintf('  %s : %d\n',u(k),sum(labels==u(k)));end,end
observed_ok=(Nrows>0)&&(NfiniteBoth==Nrows)&&(maxerr<=1e-10)&&isempty(bad);
pass=(Ncube==1230)&&(Ncomplete==1230)&&(Nrows==184500)&&(NfiniteBoth==184500)&&observed_ok;
if pass
    fprintf('STATUS: HIER4_RERUN_CONTRACT_CONFIRMED\n');
elseif observed_ok
    fprintf('STATUS: PARTIAL_HIER4_CONTRACT_OK - RERUN INCOMPLETE, TRUTH REMAINS SEALED\n');
else
    fprintf('STATUS: HIER4_CONTRACT_FAILURE - STOP, DO NOT OPEN TRUTH\n');
end
fprintf('============================================================\n');
OUT=struct('pass',pass,'observed_ok',observed_ok,'N_cubes',Ncube,'N_complete',Ncomplete,'N_rows',Nrows,'N_finite_both',NfiniteBoth,'max_identity_error',maxerr,'bad_files',bad);
end
