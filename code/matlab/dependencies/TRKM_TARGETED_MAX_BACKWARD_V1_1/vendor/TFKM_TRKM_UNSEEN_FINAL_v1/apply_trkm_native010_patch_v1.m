function info=apply_trkm_native010_patch_v1(output_root)
%APPLY_TRKM_NATIVE010_PATCH_V1 Register PATCH1 in an EXISTING fresh output.
% Does not alter/recompute candidate fits, DGP certification, seeds, or data.
output_root=char(string(output_root)); f=fullfile(output_root,'study.mat');
assert(isfile(f),'Missing study.mat: %s',f);
S=load(f); assert(isfield(S,'study'),'study.mat has no study struct.');
study=S.study;
assert(startsWith(string(study.version),"final_fresh_validation_v1"),'Unexpected study version: %s',string(study.version));
study.version="final_fresh_validation_v1_patch1";
study.frozen_trkm_native010="NATIVE010_ONE_STEP";
study.native010_R2_tolerance=.01;
study.native010_iterated=false;
study.patch1_basis="Consumed Confirm4 offline audit; at adoption ~15 fresh candidate rows existed, all REG_TFKM_HD; no TRKM fresh candidate or complete 150-row unit/final selection was used.";
study.patch1_candidate_fits_reusable=true;
S.study=study;
save(f,'-struct','S','-v7.3');
marker=fullfile(output_root,'PIPELINE_PATCH1_TRKM_NATIVE010_ONE_STEP.txt');
fid=fopen(marker,'w'); assert(fid>0,'Cannot write patch marker.');
fprintf(fid,'PATCH1: TRKM NATIVE010_ONE_STEP\n');
fprintf(fid,'Rule: after G_STRICT_W and R_STRICT_NR2, attempt at most one g->g-1 step.\n');
fprintf(fid,'Candidate r=min(current r,new g-1). Accept iff R2_candidate >= R2_ref-0.01 and RV_CV_candidate >= RV_CV_ref.\n');
fprintf(fid,'Then re-apply frozen alpha selector on the resulting (r,g).\n');
fprintf(fid,'No candidate fit, DGP, seed, multistart, or stored result is changed.\n');
fprintf(fid,'Decision basis: consumed Confirm4 offline safeguard audit; at adoption ~15 fresh rows existed, all REG_TFKM_HD; no TRKM fresh candidate/final selection was used.\n');
fclose(fid);
info=struct('study_file',string(f),'marker',string(marker),'version',study.version,'candidate_fits_reusable',true);
fprintf('PATCH1 REGISTERED: TRKM NATIVE010_ONE_STEP. Existing candidate fits remain valid and reusable.\n');
end
