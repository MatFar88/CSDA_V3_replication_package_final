function study=setup_unseen_validation_v1(output_root)
%SETUP_UNSEEN_VALIDATION_V1 Create one NEW root for TFKM + TRKM unseen validation.
%
% Frozen confirmation design:
% - same four DGP definitions as development;
% - new scenario IDs 33001/33004/33005/33007;
% - new seed namespace 1618033;
% - 8 independent replicates;
% - TFKM core = HIER_STABLE4;
% - TRKM core = unchanged RV_CV;
% - final metaparameter selectors are NOT the historical package selectors.

pkg=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(pkg),'-begin');
output_root=char(string(output_root));

assert(~contains(lower(output_root),'final_fresh_validation_v1_output'), ...
    'Use a NEW unseen output root, not the development output folder.');

if ~isfolder(output_root), mkdir(output_root); end
for d={ 'candidates','candidates/TFKM','candidates/TRKM', ...
        '_SEALED_TRUTH_DO_NOT_OPEN_BEFORE_SELECTION', ...
        'operational_selection','truth_evaluation'}
    p=fullfile(output_root,d{1}); if ~isfolder(p), mkdir(p); end
end

options=make_generic_study_options_v1("joint");
options.N_starts=100;
options.competitor_starts=100;
options.base_seed=1618033;
options.replicates=8;

catalog=build_unseen_catalog_v1(options);
methods=table(["REG_TFKM_HD";"REG_TRKM_HD"],'VariableNames',{'method'});

pr=[]; pg=[];
for g=2:6
    for r=1:(g-1)
        pr(end+1,1)=r; pg(end+1,1)=g; %#ok<AGROW>
    end
end
pairs=table(pr,pg,'VariableNames',{'r','g'});
alpha_grid=table([0;.025;.05;.075;.10;.125;.15;.20;.25;.30], ...
    'VariableNames',{'alpha'});

validation_reps=(1:8)';
sessions=string(char('A'+(0:7))');
expected_per_method=height(catalog)*height(pairs)*height(alpha_grid); % 600
manifest=table(sessions,validation_reps, ...
    repmat(expected_per_method,numel(validation_reps),1), ...
    repmat(2*expected_per_method,numel(validation_reps),1), ...
    'VariableNames',{'session','replicate','expected_per_method','expected_total'});

study=struct();
study.version="TFKM_TRKM_UNSEEN_FINAL_v1";
study.package_root=string(pkg);
study.output_root=string(output_root);
study.validation_base_seed=1618033;
study.validation_reps=validation_reps;
study.sessions=sessions;
study.N_starts=100;
study.tfkm_internal_selector="HIER_STABLE4";
study.trkm_internal_selector="RV_CV";
study.tfkm_structural_rule="SUPPORTED_WEAK";
study.trkm_structural_rule="RV_R2_LEX_CONSERVATIVE";
study.tfkm_regime_thresholds=[0.01 0.05];
study.trkm_break_threshold=0.01;
study.alpha_grid=double(alpha_grid.alpha)';
study.truth_firewall=true;
study.truth_open_allowed=false;
study.created=string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

save(fullfile(output_root,'study.mat'), ...
    'study','options','catalog','methods','pairs','alpha_grid','manifest','-v7.3');
writetable(manifest,fullfile(output_root,'manifest.csv'));

fprintf('\n============================================================\n');
fprintf('TFKM/TRKM UNSEEN FINAL v1 CREATED\n');
fprintf('Root: %s\n',output_root);
fprintf('New scenarios: %s\n',mat2str(double(catalog.scenario_id)'));
fprintf('New seed namespace: %d + sid*100000 + rep\n',study.validation_base_seed);
fprintf('Sessions: %s | 600 TFKM + 600 TRKM candidates/session\n',strjoin(sessions',','));
fprintf('IMPORTANT: do not inspect _SEALED_TRUTH... before operational selection is frozen.\n');
fprintf('============================================================\n');
end
