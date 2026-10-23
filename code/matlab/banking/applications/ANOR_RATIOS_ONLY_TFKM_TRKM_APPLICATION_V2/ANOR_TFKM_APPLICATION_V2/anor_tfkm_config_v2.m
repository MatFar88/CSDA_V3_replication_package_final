function CFG = anor_tfkm_config_v2(package_root)
%ANOR_TFKM_CONFIG_V2 Fixed protocol for ANOR ratios-only TFKM application.
% Source data_ECB_ANOR.mat is 365x383: column 1 is log(total assets),
% columns 2:383 are the 382 ratios used here. Raw duplicate copies 93 and
% 108 are removed after exact-pair verification (93==99, 108==126), giving
% a 363x382 ratios-only matrix.

if nargin<1 || strlength(string(package_root))==0
    package_root=fileparts(mfilename('fullpath'));
end
package_root=char(string(package_root));
CFG=struct();
CFG.version="ANOR_TFKM_APPLICATION_V2_RATIOS_ONLY";
CFG.package_root=string(package_root);
CFG.vendor_root=string(fullfile(package_root,'vendor'));
CFG.output_root=string(fullfile(package_root,'output'));
CFG.input_file=string(fullfile(CFG.output_root,'ANOR363x382_RATIOS_ONLY_INPUT_INTERNAL.mat'));

CFG.n_raw=365;
CFG.p_raw=383;
CFG.drop_source_column=1;
CFG.p=382;
CFG.duplicate_pairs_raw=[93 99;108 126];
CFG.drop_raw_rows=[93 108];
CFG.n=363;
CFG.dataset="ANOR_Q4_2014_363x382_RATIOS_ONLY";

CFG.alpha_grid=[0 .025 .05 .075 .10 .125 .15 .20 .25 .30];
P=zeros(15,2);z=0;
for g=2:6
    for r=1:(g-1),z=z+1;P(z,:)=[r g];end
end
CFG.pairs=P;
CFG.pair_seed=20261013;

CFG.options=struct();
CFG.options.N_starts=100;
CFG.options.hd_init_names=["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"];
CFG.options.hd_tfkm_allocation=[50 35 10 5];
CFG.options.max_iterations=50;
CFG.options.tauW=1e-5;
CFG.options.tauT=1e-3;
CFG.options.trim_stability_weight=.2;

CFG.tol=1e-12;
CFG.ari_cutoff=.90;
CFG.trim_weight=.2;
CFG.tfkm_Smed_cutoff=.44;
CFG.tfkm_DS_cutoff=.12;
CFG.tfkm_J_cutoff=.50;
CFG.tfkm_score_tol=5e-8;
CFG.obj_rel_tol=5e-8;

CFG.candidate_root=string(fullfile(CFG.output_root,'candidates'));
CFG.final_root=string(fullfile(CFG.output_root,'final'));
CFG.descriptive_root=string(fullfile(CFG.output_root,'descriptive_refit'));
end
