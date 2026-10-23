root=fileparts(mfilename('fullpath'));
addpath(genpath(root),'-begin');
fprintf('TFKM/TRKM FINAL FRESH VALIDATION v1 loaded from:\n%s\n',root);
fprintf('Frozen HD core + Stage1/Stage2 + TRKM NATIVE010_ONE_STEP PATCH1 + alpha. New DGP IDs and seed namespaces.\n');
key={'run_hd_multistart_top_v1','regularized_trimmed_fit_corrected_v1','select_joint_rg_vote_then_alpha_v1','select_alpha_performance_frozen_v1'};
for j=1:numel(key)
    p=which(key{j});
    assert(startsWith(string(p),string(root)),'Path shadowing: %s resolves to %s',key{j},p);
end
