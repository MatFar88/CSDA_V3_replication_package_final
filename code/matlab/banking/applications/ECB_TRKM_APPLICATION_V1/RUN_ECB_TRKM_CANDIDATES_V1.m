function OUT = RUN_ECB_TRKM_CANDIDATES_V1(worker_id,n_workers)
%RUN_ECB_TRKM_CANDIDATES_V1 Fit the 150 locked TRKM candidates.
%
% Restartable and safe for independent MATLAB sessions.
% Each candidate uses the exact 100-start RV-CV production selector.

if nargin<1,worker_id=1;end
if nargin<2,n_workers=1;end
worker_id=double(worker_id);n_workers=double(n_workers);
assert(worker_id>=1 && worker_id<=n_workers && worker_id==floor(worker_id) && n_workers==floor(n_workers), ...
    'Require integer 1 <= worker_id <= n_workers.');

root=fileparts(mfilename('fullpath'));
CFG=ecb_trkm_config_v1(root);
ecb_trkm_path_lock_v1(root);
assert(isfile(CFG.input_file),'Build the input first with BUILD_ECB363_INPUT_TRKM_V1.');
S=load(CFG.input_file,'X','protocol_version','pair_seed','alpha_grid','pairs','options');
X=double(S.X);
assert(isequal(size(X),[CFG.n CFG.p]),'Input dimension drift.');
assert(string(S.protocol_version)==CFG.version,'Input protocol version drift.');
assert(double(S.pair_seed)==CFG.pair_seed,'Pair-seed drift.');
assert(isequal(double(S.alpha_grid),CFG.alpha_grid),'Alpha-grid drift.');
assert(isequal(double(S.pairs),CFG.pairs),'Structural-grid drift.');
options=S.options;check_options_local(options,CFG);

if ~isfolder(CFG.candidate_root),mkdir(CFG.candidate_root);end
Jobs=build_jobs_local(CFG);idx=(1:height(Jobs))';Jobs=Jobs(mod(idx-1,n_workers)+1==worker_id,:);

fprintf('\n============================================================\n');
fprintf('ECB TRKM CANDIDATES V1 | worker %d/%d | jobs=%d\n',worker_id,n_workers,height(Jobs));
fprintf('Input: 363 x 361 | original log-level scales retained\n');
fprintf('Grid: 15 (r,g) pairs x 10 alpha values = 150 candidates\n');
fprintf('Starts: 40/30/15/15 | RV-CV | pair seed=%d\n',CFG.pair_seed);
fprintf('============================================================\n');

Ndone=0;Nnew=0;Ninvalid=0;
for j=1:height(Jobs)
    k=double(Jobs.candidate_index(j));r=double(Jobs.r(j));g=double(Jobs.g(j));a=double(Jobs.alpha(j));
    outfile=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));
    donefile=fullfile(CFG.candidate_root,sprintf('candidate_%03d.COMPLETE.txt',k));
    if isfile(outfile)&&isfile(donefile)
        Qold=load(outfile,'Row','CandidateMeta');
        assert(isfield(Qold,'CandidateMeta'),'Existing candidate %d predates application contract; move/delete it.',k);
        assert(string(Qold.CandidateMeta.protocol_version)==CFG.version && double(Qold.CandidateMeta.pair_seed)==CFG.pair_seed && ...
            double(Qold.CandidateMeta.n)==CFG.n && double(Qold.CandidateMeta.p)==CFG.p && double(Qold.CandidateMeta.drop_raw_row)==CFG.drop_raw_row, ...
            'Existing candidate %d belongs to a different application contract.',k);
        assert(double(Qold.Row.candidate_index)==k && double(Qold.Row.r)==r && double(Qold.Row.g)==g && ...
            abs(double(Qold.Row.alpha)-a)<=CFG.tol && double(Qold.Row.seed)==CFG.pair_seed, ...
            'Existing candidate %d state/seed mismatch.',k);
        Ndone=Ndone+1;fprintf('[%3d/150] (%d,%d,%.3f) already COMPLETE\n',k,r,g,a);continue
    end

    fprintf('[%3d/150] TRKM (%d,%d,%.3f) ... ',k,r,g,a);
    t0=tic;
    fit=run_hd_multistart_top_v1(X,r,g,a,"REG_TRKM_HD",options,CFG.pair_seed);
    runtime_seconds=toc(t0);
    Row=compute_operational_candidate_row_unseen_v1(X,fit,0,1,CFG.pair_seed,"REG_TRKM_HD",r,g,a,runtime_seconds);
    Row=addvars(Row,k,'Before',1,'NewVariableNames','candidate_index');Row.dataset="ECB_Q4_2014_363x361";

    if fit.valid
        Labels=uint16(fit.labels(:));Trim=logical(fit.trim(:));U=double(fit.U);
        assert(numel(Labels)==CFG.n && numel(Trim)==CFG.n,'Fit size mismatch.');
        assert(size(U,1)==CFG.p && size(U,2)==r,'Subspace size mismatch.');
        fprintf('RV=%.6f | R2=%.6f | %s | start=%d | %.1fs\n', ...
            double(fit.selector_cv_score),double(Row.R2),string(fit.selected_init_family),double(fit.selected_start_id),runtime_seconds);
    else
        Labels=zeros(CFG.n,1,'uint16');Trim=false(CFG.n,1);U=zeros(CFG.p,0);Ninvalid=Ninvalid+1;
        fprintf('INVALID | %s | %.1fs\n',string(fit.message),runtime_seconds);
    end

    CandidateMeta=struct('protocol_version',CFG.version,'pair_seed',CFG.pair_seed,'n',CFG.n,'p',CFG.p, ...
        'drop_raw_row',CFG.drop_raw_row,'external_standardization',"NONE_LOG_LEVELS_RETAINED");
    save(outfile,'Row','Labels','Trim','U','fit','CandidateMeta','-v7.3');
    fid=fopen(donefile,'w');assert(fid>0);fprintf(fid,'COMPLETE\ncandidate_index=%d\nr=%d\ng=%d\nalpha=%.12g\nvalid=%d\n',k,r,g,a,logical(fit.valid));fclose(fid);
    Nnew=Nnew+1;
end
Summary=collect_worker_rows_local(Jobs,CFG);
summary_file=fullfile(CFG.candidate_root,sprintf('worker_%02d_of_%02d_summary.csv',worker_id,n_workers));writetable(Summary,summary_file);
fprintf('\nWorker %d/%d finished: pre-existing=%d, newly fitted=%d, invalid_new=%d.\n',worker_id,n_workers,Ndone,Nnew,Ninvalid);
OUT=struct('worker_id',worker_id,'n_workers',n_workers,'summary_file',string(summary_file));
end

function Jobs=build_jobs_local(CFG)
Rows=zeros(150,4);z=0;
for ip=1:size(CFG.pairs,1)
    r=CFG.pairs(ip,1);g=CFG.pairs(ip,2);
    for ia=1:numel(CFG.alpha_grid),z=z+1;Rows(z,:)=[z r g CFG.alpha_grid(ia)];end
end
Jobs=array2table(Rows,'VariableNames',{'candidate_index','r','g','alpha'});
end
function T=collect_worker_rows_local(Jobs,CFG)
T=table();for j=1:height(Jobs),k=double(Jobs.candidate_index(j));f=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));if ~isfile(f),continue,end;Q=load(f,'Row');if isempty(T),T=Q.Row;else,T=[T;Q.Row];end,end %#ok<AGROW>
if ~isempty(T),T=sortrows(T,'candidate_index');end
end
function check_options_local(options,CFG)
assert(isfield(options,'hd_trkm_allocation')&&isfield(options,'hd_init_names')&&isfield(options,'tauW')&&isfield(options,'tauT')&&isfield(options,'max_iterations')&&isfield(options,'trim_stability_weight'),'Input options lack required fields.');
assert(isequal(double(options.hd_trkm_allocation(:))',CFG.options.hd_trkm_allocation),'TRKM start allocation drift.');
assert(isequal(string(options.hd_init_names(:))',string(CFG.options.hd_init_names(:))'),'Initialization-family drift.');
assert(double(options.max_iterations)==CFG.options.max_iterations,'max_iterations drift.');
assert(abs(double(options.tauW)-CFG.options.tauW)<=eps,'tauW drift.');
assert(abs(double(options.tauT)-CFG.options.tauT)<=eps,'tauT drift.');
assert(abs(double(options.trim_stability_weight)-CFG.options.trim_stability_weight)<=eps,'trim weight drift.');
end
