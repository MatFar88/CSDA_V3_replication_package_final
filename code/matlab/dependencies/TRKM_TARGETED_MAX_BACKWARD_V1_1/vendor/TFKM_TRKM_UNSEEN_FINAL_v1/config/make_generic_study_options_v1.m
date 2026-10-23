function options = make_generic_study_options_v1(mode)
% Central configuration shared by FIXED, SELECT_RG, SELECT_ALPHA and JOINT.
if nargin < 1, mode = "fixed"; end
options = struct();
options.mode = lower(string(mode));
options.n = 300;
options.p_over_n = [0.1 0.3 0.5 1 2];
options.structures = [2 3; 3 4; 4 5];
options.arms = ["TFKM_favorable","TRKM_favorable"];
options.difficulty = ["difficult","moderate","easy","perfect"];
options.contamination = ["clean","leverage","bridge"];
options.alpha_by_contamination = [0 0.10 0.10];
options.replicates = 100;
options.base_seed = 8142026;
options.max_iterations = 50;
options.tauW = 1e-5;
options.tauT = 1e-3;
options.competitor_starts = 100;
options.N_starts = options.competitor_starts;
options.tclust_model = "VII";
options.refsteps = 15;
options.reftol = 1e-8;
options.tclust_batch_size = 10;
options.tclust_equalweights = false;
options.rpca_ridge = 1e-6;
options.rpca_tol = 1e-7;
options.rpca_max_iter = 200;
options.stop_on_error = false;
options.save_every = 1;
options.trim_stability_weight = 0.2;
options.hd_init_names = ["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"];
options.hd_tfkm_allocation = [50 35 10 5];
options.hd_trkm_allocation = [40 30 15 15];
options.registry_version = "calibrated_v1";
options.queue_count = 8;
options.method_names = ["TFKM_ORIGINAL","TRKM_ORIGINAL","REG_TFKM_HD", ...
    "REG_TRKM_HD","TKM","PCA_TKM","RPCA_TKM","TCLUST"];
% Windowing controls. Runtime calibration can overwrite these values.
options.target_window_minutes = 180;
options.window_safety_factor = 1.25;
options.default_replicates_per_window = 2;
end
