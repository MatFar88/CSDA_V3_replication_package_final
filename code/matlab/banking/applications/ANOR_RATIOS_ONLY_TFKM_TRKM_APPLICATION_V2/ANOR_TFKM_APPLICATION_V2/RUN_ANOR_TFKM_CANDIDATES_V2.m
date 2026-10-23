function OUT = RUN_ANOR_TFKM_CANDIDATES_V2(worker_id,n_workers)
%RUN_ANOR_TFKM_CANDIDATES_V2 Fit the 150 locked TFKM candidates.
%
% Restartable and safe for several independent MATLAB sessions.
% Examples:
%   RUN_ANOR_TFKM_CANDIDATES_V2(1,1)   % serial
%   RUN_ANOR_TFKM_CANDIDATES_V2(1,4)   % worker 1 of 4
%   RUN_ANOR_TFKM_CANDIDATES_V2(2,4)   % worker 2 of 4, separate MATLAB
%
% Each candidate uses the same dataset-level pair seed and the exact 100-start
% HIER4 source bundled under vendor/. No bank names enter this function.

if nargin<1, worker_id=1; end
if nargin<2, n_workers=1; end
worker_id=double(worker_id); n_workers=double(n_workers);
assert(worker_id>=1 && worker_id<=n_workers && worker_id==floor(worker_id) && n_workers==floor(n_workers), ...
    'Require integer 1 <= worker_id <= n_workers.');

root = fileparts(mfilename('fullpath'));
CFG = anor_tfkm_config_v2(root);
anor_tfkm_path_lock_v2(root);
assert(isfile(CFG.input_file),'Build the input first with BUILD_ANOR363_RATIOS_INPUT_V2.');
S = load(CFG.input_file,'X','protocol_version','pair_seed','alpha_grid','pairs','options');
X = double(S.X);
assert(isequal(size(X),[CFG.n CFG.p]),'Input dimension drift.');
assert(string(S.protocol_version)==CFG.version,'Input protocol version drift.');
assert(double(S.pair_seed)==CFG.pair_seed,'Pair-seed drift.');
assert(isequal(double(S.alpha_grid),CFG.alpha_grid),'Alpha-grid drift.');
assert(isequal(double(S.pairs),CFG.pairs),'Structural-grid drift.');
options = S.options;
check_options_local(options,CFG);

if ~isfolder(CFG.candidate_root), mkdir(CFG.candidate_root); end

Jobs = build_jobs_local(CFG);
idx = (1:height(Jobs))';
Jobs = Jobs(mod(idx-1,n_workers)+1==worker_id,:);

fprintf('\n============================================================\n');
fprintf('ANOR TFKM CANDIDATES V2 | worker %d/%d | jobs=%d\n',worker_id,n_workers,height(Jobs));
fprintf('Input: 363 x 382 | ratios only; source column 1 (log assets) excluded\n');
fprintf('Grid: 15 (r,g) pairs x 10 alpha values = 150 candidates\n');
fprintf('Starts: 50/35/10/5 | HIER4 | pair seed=%d\n',CFG.pair_seed);
fprintf('============================================================\n');

Ndone=0; Nnew=0; Ninvalid=0;
for j=1:height(Jobs)
    k=double(Jobs.candidate_index(j));
    r=double(Jobs.r(j)); g=double(Jobs.g(j)); a=double(Jobs.alpha(j));
    outfile=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));
    donefile=fullfile(CFG.candidate_root,sprintf('candidate_%03d.COMPLETE.txt',k));
    if isfile(outfile) && isfile(donefile)
        Qold=load(outfile,'Row','CandidateMeta');
        assert(isfield(Qold,'CandidateMeta'),'Existing candidate %d predates the application contract; move/delete it before continuing.',k);
        assert(string(Qold.CandidateMeta.protocol_version)==CFG.version && ...
            double(Qold.CandidateMeta.pair_seed)==CFG.pair_seed && ...
            double(Qold.CandidateMeta.n)==CFG.n && double(Qold.CandidateMeta.p)==CFG.p && ...
            isequal(double(Qold.CandidateMeta.drop_raw_rows(:))',CFG.drop_raw_rows), ...
            'Existing candidate %d belongs to a different application contract.',k);
        assert(double(Qold.Row.candidate_index)==k && double(Qold.Row.r)==r && ...
            double(Qold.Row.g)==g && abs(double(Qold.Row.alpha)-a)<=CFG.tol && ...
            double(Qold.Row.seed)==CFG.pair_seed, ...
            'Existing candidate %d state/seed mismatch.',k);
        Ndone=Ndone+1;
        fprintf('[%3d/150] (%d,%d,%.3f) already COMPLETE\n',k,r,g,a);
        continue
    end

    fprintf('[%3d/150] TFKM (%d,%d,%.3f) ... ',k,r,g,a);
    t0=tic;
    fit = run_hd_multistart_top_hier_stable4_unseen_v1(X,r,g,a,options,CFG.pair_seed);
    runtime_seconds=toc(t0);
    Row = compute_operational_candidate_row_unseen_v1( ...
        X,fit,0,1,CFG.pair_seed,"REG_TFKM_HD",r,g,a,runtime_seconds);
    Row = addvars(Row,k,'Before',1,'NewVariableNames','candidate_index');
    Row.dataset = CFG.dataset;

    if fit.valid
        Labels = uint16(fit.labels(:));
        Trim = logical(fit.trim(:));
        U = double(fit.U);
        assert(numel(Labels)==CFG.n && numel(Trim)==CFG.n,'Fit size mismatch.');
        assert(size(U,1)==CFG.p && size(U,2)==r,'Subspace size mismatch.');
        fprintf('H=%.6f | A=%.6f | J=%.6f | %s | %.1fs\n', ...
            double(fit.hier_cross_H_score),double(fit.hier_cross_ari_score), ...
            double(fit.hier_cross_j_score),string(fit.selected_init_family),runtime_seconds);
    else
        Labels = zeros(CFG.n,1,'uint16');
        Trim = false(CFG.n,1);
        U = zeros(CFG.p,0);
        Ninvalid=Ninvalid+1;
        fprintf('INVALID | %s | %.1fs\n',string(fit.message),runtime_seconds);
    end

    % Keep candidate outputs anonymous. Bank names remain only in the internal
    % input contract MAT file. CandidateMeta binds restart files to this
    % pre-specified application contract without storing bank identifiers.
    CandidateMeta=struct('protocol_version',CFG.version,'pair_seed',CFG.pair_seed, ...
        'n',CFG.n,'p',CFG.p,'drop_raw_rows',CFG.drop_raw_rows, ...
        'excluded_source_column',CFG.drop_source_column, ...
        'external_standardization',"NONE_RATIOS_AS_SUPPLIED");
    save(outfile,'Row','Labels','Trim','U','fit','CandidateMeta','-v7.3');
    fid=fopen(donefile,'w'); assert(fid>0);
    fprintf(fid,'COMPLETE\ncandidate_index=%d\nr=%d\ng=%d\nalpha=%.12g\nvalid=%d\n', ...
        k,r,g,a,logical(fit.valid)); fclose(fid);
    Nnew=Nnew+1;
end

Summary = collect_worker_rows_local(Jobs,CFG);
summary_file=fullfile(CFG.candidate_root,sprintf('worker_%02d_of_%02d_summary.csv',worker_id,n_workers));
writetable(Summary,summary_file);

fprintf('\nWorker %d/%d finished: pre-existing=%d, newly fitted=%d, invalid_new=%d.\n', ...
    worker_id,n_workers,Ndone,Nnew,Ninvalid);
OUT=struct('worker_id',worker_id,'n_workers',n_workers,'summary_file',string(summary_file));
end

function Jobs=build_jobs_local(CFG)
Rows=zeros(150,4); z=0;
for ip=1:size(CFG.pairs,1)
    r=CFG.pairs(ip,1); g=CFG.pairs(ip,2);
    for ia=1:numel(CFG.alpha_grid)
        z=z+1; Rows(z,:)=[z r g CFG.alpha_grid(ia)];
    end
end
Jobs=array2table(Rows,'VariableNames',{'candidate_index','r','g','alpha'});
end

function T=collect_worker_rows_local(Jobs,CFG)
T=table();
for j=1:height(Jobs)
    k=double(Jobs.candidate_index(j));
    f=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));
    if ~isfile(f),continue,end
    Q=load(f,'Row');
    if isempty(T),T=Q.Row;else,T=[T;Q.Row];end %#ok<AGROW>
end
if ~isempty(T),T=sortrows(T,'candidate_index');end
end

function check_options_local(options,CFG)
assert(isfield(options,'hd_tfkm_allocation') && isfield(options,'hd_init_names') && ...
    isfield(options,'tauW') && isfield(options,'tauT') && ...
    isfield(options,'max_iterations') && isfield(options,'trim_stability_weight'), ...
    'Input options lack required fields.');
assert(isequal(double(options.hd_tfkm_allocation(:))',CFG.options.hd_tfkm_allocation), ...
    'TFKM start allocation drift.');
assert(isequal(string(options.hd_init_names(:))',string(CFG.options.hd_init_names(:))'), ...
    'Initialization-family drift.');
assert(double(options.max_iterations)==CFG.options.max_iterations,'max_iterations drift.');
assert(abs(double(options.tauW)-CFG.options.tauW)<=eps,'tauW drift.');
assert(abs(double(options.tauT)-CFG.options.tauT)<=eps,'tauT drift.');
assert(abs(double(options.trim_stability_weight)-CFG.options.trim_stability_weight)<=eps,'trim weight drift.');
end
