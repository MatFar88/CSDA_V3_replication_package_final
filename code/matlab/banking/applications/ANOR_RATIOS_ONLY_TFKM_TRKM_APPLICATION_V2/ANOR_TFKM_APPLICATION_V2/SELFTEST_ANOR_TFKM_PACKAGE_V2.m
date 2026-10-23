function OUT = SELFTEST_ANOR_TFKM_PACKAGE_V2(run_fit_test)
%SELFTEST_ANOR_TFKM_PACKAGE_V2 Static contract checks + optional small fit test.
if nargin<1,run_fit_test=false;end
root=fileparts(mfilename('fullpath'));
CFG=anor_tfkm_config_v2(root);
paths=anor_tfkm_path_lock_v2(root);
assert(isequal(CFG.pairs,[1 2;1 3;2 3;1 4;2 4;3 4;1 5;2 5;3 5;4 5;1 6;2 6;3 6;4 6;5 6]),'Pair-grid ordering drift.');
assert(isequal(CFG.alpha_grid,[0 .025 .05 .075 .10 .125 .15 .20 .25 .30]),'Alpha-grid drift.');
assert(isequal(CFG.options.hd_tfkm_allocation,[50 35 10 5]),'Start-allocation drift.');
assert(CFG.options.max_iterations==50 && CFG.options.tauW==1e-5 && CFG.options.tauT==1e-3,'Core option drift.');
assert(CFG.pair_seed==20261013,'Pair seed drift.');
assert(CFG.n_raw==365 && CFG.p_raw==383 && CFG.n==363 && CFG.p==382,'ANOR dimension contract drift.');
assert(CFG.drop_source_column==1 && isequal(CFG.drop_raw_rows,[93 108]),'ANOR source exclusion/deletion contract drift.');
assert(isequal(CFG.duplicate_pairs_raw,[93 99;108 126]),'ANOR duplicate-pair contract drift.');

fit_test_pass=true;
if run_fit_test
    rng(123,'twister');
    X=[randn(15,8);randn(15,8)+2];
    fit=run_hd_multistart_top_hier_stable4_unseen_v1(X,1,2,0,CFG.options,CFG.pair_seed);
    fit_test_pass=fit.valid && isfinite(fit.selector_score) && ...
        abs(fit.hier_cross_H_score-(fit.hier_cross_ari_score+.2*fit.hier_cross_j_score))<1e-10;
    assert(fit_test_pass,'Optional small fit test failed.');
end
fprintf('SELFTEST PASS: ANOR TFKM package contracts%s.\n', ...
    ternary_local(run_fit_test,' + small fit',''));
OUT=struct('pass',true,'fit_test_pass',fit_test_pass,'paths',paths);
end
function s=ternary_local(tf,a,b),if tf,s=a;else,s=b;end,end
