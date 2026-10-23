function OUT = BUILD_ANOR363_RATIOS_INPUT_V2(data_file,labels_file,data_banks_file,names_file)
%BUILD_ANOR363_RATIOS_INPUT_V2 Build the locked 363x382 ANOR ratios-only input.
% data_ECB_ANOR.mat must contain DataECB (365x383). Column 1 is log assets
% and is EXCLUDED from fitting. Columns 2:383 are the 382 ratios.
% Exact duplicate pairs must be 93==99 and 108==126. Raw rows 93 and 108
% are dropped simultaneously, yielding 363 unique observations.

root=fileparts(mfilename('fullpath'));CFG=anor_tfkm_config_v2(root);
if ~isfolder(CFG.output_root),mkdir(CFG.output_root);end
if nargin<3,data_banks_file="";end
if nargin<4,names_file="";end
assert(isfile(data_file),'data_ECB_ANOR.mat not found: %s',data_file);
assert(isfile(labels_file),'labels_ECB_ANOR.mat not found: %s',labels_file);
D=load(data_file,'DataECB');L=load(labels_file,'labelsPlus');
assert(isfield(D,'DataECB'),'data_ECB_ANOR.mat must contain DataECB.');
assert(isfield(L,'labelsPlus'),'labels_ECB_ANOR.mat must contain labelsPlus.');
Xfull=double(D.DataECB);
assert(isequal(size(Xfull),[CFG.n_raw CFG.p_raw]),'Expected DataECB %dx%d.',CFG.n_raw,CFG.p_raw);
assert(all(isfinite(Xfull(:))),'DataECB contains nonfinite values.');
assert(numel(L.labelsPlus)==CFG.p,'Expected %d ratio labels.',CFG.p);

% Verify the two exact source duplicate pairs before any deletion.
for k=1:size(CFG.duplicate_pairs_raw,1)
    a=CFG.duplicate_pairs_raw(k,1);b=CFG.duplicate_pairs_raw(k,2);
    assert(isequal(Xfull(a,:),Xfull(b,:)),'Expected exact duplicate raw rows %d and %d.',a,b);
end
assert(size(unique(Xfull,'rows'),1)==CFG.n,'Source must contain exactly %d unique rows.',CFG.n);

% Critical representation lock: column 1 is deliberately NOT part of X.
log_assets_raw=Xfull(:,CFG.drop_source_column); %#ok<NASGU>
Xratio_raw=Xfull(:,2:end);
assert(size(Xratio_raw,2)==CFG.p,'Ratios-only column count mismatch.');
keep=true(CFG.n_raw,1);keep(CFG.drop_raw_rows)=false;
X=Xratio_raw(keep,:);
original_row_index=find(keep);
Lab_vars=L.labelsPlus;
assert(isequal(size(X),[CFG.n CFG.p]),'Clean ratios-only input size mismatch.');
assert(size(unique(X,'rows'),1)==CFG.n,'Duplicate row remains after prescribed deletions.');

% Optional but strong historical row-order check against data_banks.mat.
alignment_verified=false;alignment_shift=NaN;alignment_max_dev=NaN;
if strlength(string(data_banks_file))>0
    assert(isfile(data_banks_file),'Optional data_banks.mat not found: %s',data_banks_file);
    B=load(data_banks_file,'DataECB_def');
    assert(isfield(B,'DataECB_def')&&size(B.DataECB_def,1)==CFG.n,'data_banks.mat row contract mismatch.');
    cleaned_log_assets=Xfull(keep,1);
    d=double(B.DataECB_def(:,1))-cleaned_log_assets;
    alignment_shift=median(d);alignment_max_dev=max(abs(d-alignment_shift));
    assert(alignment_max_dev<1e-10,'Cleaned ANOR row order does not align with historical 363-row banking file.');
    alignment_verified=true;
end
names_verified=false;
if strlength(string(names_file))>0
    assert(isfile(names_file),'Optional names_banks.mat not found: %s',names_file);
    N=load(names_file,'names_bank_top');
    assert(isfield(N,'names_bank_top')&&numel(N.names_bank_top)==CFG.n,'Expected %d final bank-name entries.',CFG.n);
    names_verified=true;
end

external_standardization="NONE_RATIOS_AS_SUPPLIED";
source_column_1_status="EXCLUDED_LOG_ASSETS";
protocol_version=CFG.version;pair_seed=CFG.pair_seed;alpha_grid=CFG.alpha_grid;pairs=CFG.pairs;options=CFG.options;
save(CFG.input_file,'X','Lab_vars','original_row_index','alignment_verified','alignment_shift','alignment_max_dev', ...
    'names_verified','external_standardization','source_column_1_status','protocol_version','pair_seed','alpha_grid','pairs','options','-v7.3');

Audit=table(CFG.n_raw,CFG.p_raw,CFG.drop_source_column,CFG.n,CFG.p,CFG.drop_raw_rows(1),CFG.drop_raw_rows(2), ...
    CFG.duplicate_pairs_raw(1,1),CFG.duplicate_pairs_raw(1,2),CFG.duplicate_pairs_raw(2,1),CFG.duplicate_pairs_raw(2,2), ...
    alignment_verified,alignment_shift,alignment_max_dev,CFG.p/CFG.n,string(external_standardization),double(pair_seed), ...
    'VariableNames',{'n_raw','p_raw','excluded_source_column','n_clean','p_ratios','dropped_raw_row_1','dropped_raw_row_2', ...
    'duplicate_pair1_a','duplicate_pair1_b','duplicate_pair2_a','duplicate_pair2_b','historical_alignment_verified', ...
    'historical_logasset_shift','historical_alignment_max_dev','p_over_n_clean','external_standardization','pair_seed'});
writetable(Audit,fullfile(CFG.output_root,'INPUT_CONTRACT_ANONYMOUS.csv'));

fprintf('\nANOR TFKM ratios-only input contract built successfully.\n');
fprintf('  source DataECB: %d x %d\n',size(Xfull,1),size(Xfull,2));
fprintf('  source column 1: EXCLUDED (log assets)\n');
fprintf('  ratio block used for fitting: columns 2:383\n');
fprintf('  raw duplicate pairs verified: 93==99 and 108==126 exactly\n');
fprintf('  dropped raw rows: 93 and 108\n');
fprintf('  X clean: %d x %d\n',size(X,1),size(X,2));
if alignment_verified,fprintf('  historical 363-row order check: PASS (constant log-asset shift; max dev %.3g)\n',alignment_max_dev);end
fprintf('  external standardization: NONE; ratios used as supplied\n');
fprintf('  pair seed: %d\n',CFG.pair_seed);
fprintf('  internal input: %s\n',CFG.input_file);
OUT=struct('input_file',CFG.input_file,'audit',Audit);
end
