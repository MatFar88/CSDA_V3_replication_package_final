function OUT = BUILD_ECB363_INPUT_TRKM_V1(ecb_mat_file,names364_file)
%BUILD_ECB363_INPUT_TRKM_V1 Build the locked 363 x 361 confidential ECB input.
%
% No external standardization is applied. Raw row 93 is removed only after
% verifying exact equality with raw row 99. Bank-name contents are never
% copied into generated outputs.

root = fileparts(mfilename('fullpath'));
CFG = ecb_trkm_config_v1(root);
if ~isfolder(CFG.output_root), mkdir(CFG.output_root); end

assert(isfile(ecb_mat_file),'ECB MAT file not found: %s',ecb_mat_file);
if nargin<2, names364_file=""; end

D = load(ecb_mat_file,'X','C','Lab_vars');
assert(isfield(D,'X') && isfield(D,'C') && isfield(D,'Lab_vars'), ...
    'ECB input must contain X, C and Lab_vars.');
Xraw = double(D.X);
assert(isequal(size(Xraw),[CFG.n_raw CFG.p]), ...
    'Unexpected X size: expected %d x %d.',CFG.n_raw,CFG.p);
assert(all(isfinite(Xraw(:))),'X contains nonfinite values.');
assert(isequal(size(D.C),[CFG.p CFG.p]),'Unexpected covariance size.');
assert(numel(D.Lab_vars)==CFG.p,'Unexpected number of variable labels.');

rel_cov = norm(cov(Xraw)-double(D.C),'fro')/max(norm(double(D.C),'fro'),eps);
assert(rel_cov < 1e-10,'Archived covariance does not reproduce cov(X).');

a = CFG.duplicate_rows_raw(1); b = CFG.duplicate_rows_raw(2);
assert(isequal(Xraw(a,:),Xraw(b,:)), ...
    'Expected exact duplicate raw rows %d and %d were not found.',a,b);

names_verified = false;
if strlength(string(names364_file))>0
    assert(isfile(names364_file),'Optional 364-name MAT file not found: %s',names364_file);
    N = load(names364_file);
    assert(isfield(N,'names_bank_top'),'names_bank_364.mat must contain names_bank_top.');
    assert(numel(N.names_bank_top)==CFG.n_raw,'Expected 364 bank-name entries, found %d.',numel(N.names_bank_top));
    names_verified = true;
end

keep = true(CFG.n_raw,1); keep(CFG.drop_raw_row)=false;
X = Xraw(keep,:);
original_row_index = find(keep);
Lab_vars = D.Lab_vars;
assert(isequal(size(X),[CFG.n CFG.p]),'Clean input size mismatch.');
assert(size(unique(X,'rows'),1)==CFG.n,'Exact duplicate row remains after prescribed deletion.');

external_standardization = "NONE_LOG_LEVELS_RETAINED";
protocol_version = CFG.version;
pair_seed = CFG.pair_seed;
alpha_grid = CFG.alpha_grid;
pairs = CFG.pairs;
options = CFG.options;

save(CFG.input_file,'X','Lab_vars','original_row_index','names_verified', ...
    'external_standardization','protocol_version','pair_seed','alpha_grid', ...
    'pairs','options','rel_cov','-v7.3');

Audit = table(CFG.n_raw,CFG.n,CFG.p,CFG.drop_raw_row,a,b,rel_cov, ...
    CFG.p/CFG.n,string(external_standardization),double(pair_seed), ...
    'VariableNames',{'n_raw','n_clean','p','dropped_raw_row', ...
    'duplicate_raw_row_1','duplicate_raw_row_2','archived_cov_rel_fro_error', ...
    'p_over_n_clean','external_standardization','pair_seed'});
writetable(Audit,fullfile(CFG.output_root,'INPUT_CONTRACT_ANONYMOUS.csv'));

fprintf('\nECB TRKM input contract built successfully.\n');
fprintf('  X clean: %d x %d\n',size(X,1),size(X,2));
fprintf('  dropped raw row: %d\n',CFG.drop_raw_row);
fprintf('  raw duplicate verified: %d == %d exactly\n',a,b);
fprintf('  external standardization: NONE (log-level scales retained)\n');
fprintf('  pair seed: %d (same application-level seed as TFKM)\n',CFG.pair_seed);
fprintf('  internal input: %s\n',CFG.input_file);
OUT = struct('input_file',CFG.input_file,'audit',Audit);
end
