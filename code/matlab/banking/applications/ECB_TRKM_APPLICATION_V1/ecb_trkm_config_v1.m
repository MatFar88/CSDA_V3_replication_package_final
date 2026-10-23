function CFG = ecb_trkm_config_v1(package_root)
%ECB_TRKM_CONFIG_V1 Fixed protocol for the ECB 363 x 361 TRKM application.
%
% This configuration is data-independent and must not be changed after
% inspecting TRKM clustering results. It deliberately reuses the application
% seed already fixed for the companion TFKM analysis.

if nargin<1 || strlength(string(package_root))==0
    package_root = fileparts(mfilename('fullpath'));
end
package_root = char(string(package_root));

CFG = struct();
CFG.version = "ECB_TRKM_APPLICATION_V1";
CFG.package_root = string(package_root);
CFG.vendor_root = string(fullfile(package_root,'vendor'));
CFG.output_root = string(fullfile(package_root,'output'));
CFG.input_file = string(fullfile(CFG.output_root,'ECB363_INPUT_INTERNAL.mat'));

% Primary application sample contract.
CFG.n_raw = 364;
CFG.p = 361;
CFG.duplicate_rows_raw = [93 99];
CFG.drop_raw_row = 93;
CFG.n = 363;

% Exact finite grid used in the paper.
CFG.alpha_grid = [0 .025 .05 .075 .10 .125 .15 .20 .25 .30];
P = zeros(15,2); z = 0;
for g=2:6
    for r=1:(g-1)
        z=z+1; P(z,:)=[r g];
    end
end
CFG.pairs = P;

% Reuse the application-level seed already fixed before the TFKM run.
% This introduces no additional seed choice after seeing TFKM results.
CFG.pair_seed = 20261013;

% Fixed multistart / ridge options from the paper contract.
CFG.options = struct();
CFG.options.N_starts = 100;
CFG.options.hd_init_names = ["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"];
CFG.options.hd_tfkm_allocation = [50 35 10 5]; % retained for bundled-core compatibility
CFG.options.hd_trkm_allocation = [40 30 15 15];
CFG.options.max_iterations = 50;
CFG.options.tauW = 1e-5;
CFG.options.tauT = 1e-3;
CFG.options.trim_stability_weight = .2;

% Frozen TRKM finalization constants.
CFG.tol = 1e-12;
CFG.ari_cutoff = .90;
CFG.trkm_onset_tol = .01;
CFG.trkm_delta_cutoff = .025;
CFG.trkm_score_tol = 1e-8;
CFG.obj_rel_tol = 5e-8;
CFG.pos_grid = CFG.alpha_grid(CFG.alpha_grid>0);

% Output locations.
CFG.candidate_root = string(fullfile(CFG.output_root,'candidates'));
CFG.final_root = string(fullfile(CFG.output_root,'final'));
CFG.descriptive_root = string(fullfile(CFG.output_root,'descriptive_refit'));

end
