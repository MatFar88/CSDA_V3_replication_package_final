function study=setup_final_fresh_validation_v1(output_root)
%SETUP_FINAL_FRESH_VALIDATION_V1 Build a completely new validation experiment.
pkg=fileparts(mfilename('fullpath')); addpath(genpath(pkg),'-begin');
output_root=char(string(output_root)); if ~isfolder(output_root),mkdir(output_root);end
for d={'results','summary','dgp_certification','forensics'}
    p=fullfile(output_root,d{1}); if ~isfolder(p),mkdir(p);end
end
options=make_generic_study_options_v1("joint");
options.N_starts=100; options.competitor_starts=100;
options.base_seed=314159; % validation namespace only
options.replicates=8;
catalog=build_final_fresh_catalog_v1(options);
methods=table(["REG_TFKM_HD";"REG_TRKM_HD"],'VariableNames',{'method'});
pr=[]; pg=[];
for g=2:6
    for r=1:(g-1), pr(end+1,1)=r; pg(end+1,1)=g; end %#ok<AGROW>
end
pairs=table(pr,pg,'VariableNames',{'r','g'});
alpha_grid=table([0;.025;.05;.075;.10;.125;.15;.20;.25;.30],'VariableNames',{'alpha'});
validation_reps=(1:8)'; certification_reps=(1:8)';
sessions=string(char('A'+(0:numel(validation_reps)-1))');
expected_per_session=height(catalog)*height(methods)*height(pairs)*height(alpha_grid);
manifest=table(sessions,validation_reps,repmat(expected_per_session,numel(validation_reps),1), ...
    'VariableNames',{'session','replicate','expected_candidate_fits'});
study=struct();
study.version="final_fresh_validation_v1_patch1"; study.package_root=string(pkg); study.output_root=string(output_root);
study.validation_base_seed=314159; study.certification_base_seed=271828;
study.validation_reps=validation_reps; study.certification_reps=certification_reps; study.sessions=sessions;
study.forensic_ari_threshold=.95; study.random_reference_B=100; study.noise_pc_count=3;
study.frozen_stage1_tfkm="CURRENT_PLUS_GMAX_R_E"; study.frozen_stage1_trkm="CURRENT";
study.frozen_stage2_g="G_STRICT_W"; study.frozen_stage2_r="R_STRICT_NR2";
study.frozen_trkm_native010="NATIVE010_ONE_STEP"; study.native010_R2_tolerance=.01; study.native010_iterated=false;
study.alpha_after_stage2=true; study.N_starts=100;
save(fullfile(output_root,'study.mat'),'study','options','catalog','methods','pairs','alpha_grid','manifest','-v7.3');
writetable(catalog,fullfile(output_root,'scenario_catalog.csv'));
writetable(manifest,fullfile(output_root,'manifest.csv'));
fprintf('\nFINAL FRESH VALIDATION v1 SET UP\n');
fprintf('Scenarios: %s\n',strjoin(string(catalog.scenario_id)',', '));
fprintf('Clean repair: sid 23001 delta=10 rho=2.\n');
fprintf('Certification seed namespace: %d + sid*100000 + rep\n',study.certification_base_seed);
fprintf('Validation seed namespace:    %d + sid*100000 + rep\n',study.validation_base_seed);
fprintf('Fresh validation reps: %s | sessions %s\n',mat2str(validation_reps'),strjoin(sessions',','));
fprintf('Candidate fits/session: %d | total: %d\n',expected_per_session,expected_per_session*numel(validation_reps));
fprintf('NEXT: certify_fresh_dgp_catalog_v1(output_root). Validation sessions refuse to start until certification passes.\n');
end
