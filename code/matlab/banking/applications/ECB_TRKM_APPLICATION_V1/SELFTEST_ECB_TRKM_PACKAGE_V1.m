function OUT = SELFTEST_ECB_TRKM_PACKAGE_V1(run_fit_test)
%SELFTEST_ECB_TRKM_PACKAGE_V1 Static contract checks + optional small TRKM fit.
if nargin<1,run_fit_test=false;end
root=fileparts(mfilename('fullpath'));
CFG=ecb_trkm_config_v1(root);
paths=ecb_trkm_path_lock_v1(root);
assert(isequal(CFG.pairs,[1 2;1 3;2 3;1 4;2 4;3 4;1 5;2 5;3 5;4 5;1 6;2 6;3 6;4 6;5 6]),'Pair-grid ordering drift.');
assert(isequal(CFG.alpha_grid,[0 .025 .05 .075 .10 .125 .15 .20 .25 .30]),'Alpha-grid drift.');
assert(isequal(CFG.options.hd_trkm_allocation,[40 30 15 15]),'TRKM start-allocation drift.');
assert(CFG.options.max_iterations==50 && CFG.options.tauW==1e-5 && CFG.options.tauT==1e-3,'Core option drift.');
assert(CFG.pair_seed==20261013,'Pair seed drift.');
assert(CFG.trkm_onset_tol==.01 && CFG.trkm_delta_cutoff==.025 && CFG.ari_cutoff==.90,'TRKM finalizer constant drift.');

fit_test_pass=true;
if run_fit_test
    rng(123,'twister');
    X=[randn(15,8);randn(15,8)+2];
    fit=run_hd_multistart_top_v1(X,1,2,0,"REG_TRKM_HD",CFG.options,CFG.pair_seed);
    fit_test_pass=fit.valid && isfinite(fit.selector_cv_score) && string(fit.selector)=="RV_CV";
    assert(fit_test_pass,'Optional small TRKM fit test failed.');
end
fprintf('SELFTEST PASS: ECB TRKM package contracts%s.\n',ternary_local(run_fit_test,' + small fit',''));
OUT=struct('pass',true,'fit_test_pass',fit_test_pass,'paths',paths);
end
function s=ternary_local(tf,a,b),if tf,s=a;else,s=b;end,end
