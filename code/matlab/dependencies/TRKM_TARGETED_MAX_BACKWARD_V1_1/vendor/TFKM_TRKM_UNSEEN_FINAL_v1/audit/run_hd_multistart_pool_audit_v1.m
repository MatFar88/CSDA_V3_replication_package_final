function A = run_hd_multistart_pool_audit_v1(sim,r,c,alpha,method,options,seed)
%RUN_HD_MULTISTART_POOL_AUDIT_V1 Refit and expose the complete frozen pool.
%
% This is a forensic/audit wrapper around the frozen v11 multistart logic.
% It does NOT alter the production estimator. All 100 starts are generated
% with the exact frozen family allocations and deterministic start seeds.
%
% Returned object A contains only audit summaries plus the compact pool;
% truth is used only for post-fit diagnostics in the start table.

method = string(method);
X = sim.X;

switch method
    case "REG_TFKM_HD"
        alloc = options.hd_tfkm_allocation;
        coremethod = "REG_TFKM";
        selector_rule = "SIGN_ALL_STABLE";
    case "REG_TRKM_HD"
        alloc = options.hd_trkm_allocation;
        coremethod = "REG_TRKM";
        selector_rule = "RV_CV";
    otherwise
        error('basinAudit:Method','Only REG_TFKM_HD and REG_TRKM_HD are supported.');
end

names = string(options.hd_init_names);
N = sum(alloc);
candidates = cell(N,1);
valid = false(N,1);
z = 0;

for j = 1:numel(names)
    for h = 1:alloc(j)
        z = z + 1;
        start_seed = seed + 1000*j + h;
        try
            rng(start_seed,'twister');
            init = make_robust_initialization_v1(X,c,r,alpha,names(j),start_seed);
            q = regularized_trimmed_fit_corrected_v1( ...
                X,init,r,c,alpha,coremethod,options.max_iterations,options.tauW,options.tauT);
            q.init_family = names(j);
            q.start_id = h;
            q.start_seed = start_seed;
            % Literal frozen production validity definition.
            q.valid = all(isfinite(q.labels)) && isfinite(q.objective) && size(q.U,2)==r;
            q.message = "";
            candidates{z} = q;
            valid(z) = q.valid;
        catch ME
            candidates{z} = struct('valid',false,'message', ...
                string(ME.identifier)+": "+string(ME.message), ...
                'init_family',names(j),'start_id',h,'start_seed',start_seed);
        end
    end
end

valid_idx = find(valid);
if isempty(valid_idx)
    error('basinAudit:NoValidStarts','No valid starts for (%d,%d,%.6f), %s.',r,c,alpha,method);
end

% Exact compact pool used by the frozen selectors.
labels_pool = zeros(size(X,1),N,'uint16');
trim_pool = false(size(X,1),N);
objective_pool = nan(N,1);
family = strings(N,1);
start_id = nan(N,1);
start_seed = nan(N,1);
for i=1:N
    q=candidates{i};
    family(i)=string(q.init_family);
    start_id(i)=q.start_id;
    start_seed(i)=q.start_seed;
    if ~valid(i), continue; end
    labels_pool(:,i)=uint16(q.labels(:));
    trim_pool(:,i)=logical(q.trim(:));
    objective_pool(i)=q.objective;
end
[~,ord] = sort(objective_pool(valid_idx),'ascend');
ordered = valid_idx(ord);
objective_rank = nan(N,1);
objective_rank(ordered) = (1:numel(ordered))';

pool=struct();
pool.valid=valid;
pool.family=family;
pool.start_id=start_id;
pool.start_seed=start_seed;
pool.objective=objective_pool;
pool.objective_rank=objective_rank;
pool.labels=labels_pool;
pool.trim=trim_pool;
pool.trim_stability_weight=options.trim_stability_weight;
pool.current_candidate_index=legacy_candidate_local(pool,valid_idx,options.trim_stability_weight);

% Frozen current selector: this MUST reproduce the stored candidate row.
if method=="REG_TFKM_HD"
    current_sel = select_tfkm_sign_v4(pool,X,c,r,selector_rule);
else
    current_sel = select_robust_view_v3(pool,X,c,r,selector_rule);
end
current_idx = current_sel.chosen;

% Common truth-blind scores for audit challengers.
rv = score_candidates_robust_views_v3(pool,X,c,r);
S = pairwise_similarity(pool);
[within_family, all_stability, family_balanced] = stability_scores(pool,S,names);
[family_winner_idx, family_winner_score] = choose_family_winners(pool,within_family,names);

% Truth-open + geometric diagnostics per terminal fit.
truth_ari = nan(N,1);
trim_precision = nan(N,1); trim_recall=nan(N,1); trim_f1=nan(N,1);
W_raw=nan(N,1); W_ratio=nan(N,1); R2=nan(N,1); explained=nan(N,1);
iterations=nan(N,1); converged=false(N,1);
labels_changed_first=nan(N,1); trim_changed_first=nan(N,1);
truth = double(sim.labels(:)); outlier = logical(sim.outlier(:)); intrue=~outlier;
for i=valid_idx(:)'
    q=candidates{i};
    truth_ari(i)=adjusted_rand_index_v1(truth(intrue),double(q.labels(intrue)));
    tr=logical(q.trim(:));
    if any(outlier)
        tp=sum(tr&outlier); fp=sum(tr&~outlier); fn=sum(~tr&outlier);
        trim_precision(i)=tp/max(tp+fp,1);
        trim_recall(i)=tp/max(tp+fn,1);
        trim_f1(i)=2*trim_precision(i)*trim_recall(i)/max(trim_precision(i)+trim_recall(i),eps);
    end
    st=compute_external_rg_fit_stats_v1(X,q,method);
    W_raw(i)=st.W_raw; W_ratio(i)=st.W_ratio; R2(i)=st.R2;
    if method=="REG_TFKM_HD", explained(i)=st.W_ratio*st.R2; else, explained(i)=st.R2; end
    if isfield(q,'iterations'), iterations(i)=q.iterations; end
    if isfield(q,'converged'), converged(i)=logical(q.converged); end
    if isfield(q,'labels_changed_first'), labels_changed_first(i)=q.labels_changed_first; end
    if isfield(q,'trim_changed_first'), trim_changed_first(i)=q.trim_changed_first; end
end

sign_score=nan(N,1);
if method=="REG_TFKM_HD" && isfield(current_sel,'score')
    sign_score=current_sel.score;
end
current_native=nan(N,1);
if method=="REG_TFKM_HD", current_native=sign_score; else, current_native=rv.cv_score; end

is_family_winner=false(N,1);
for j=1:numel(family_winner_idx)
    if isfinite(family_winner_idx(j)), is_family_winner(family_winner_idx(j))=true; end
end

start_index=(1:N)';
current_eligible=valid;
if method=="REG_TFKM_HD", current_eligible=valid & family=="SIGN_PCA_TKM"; end
is_current_winner=false(N,1); is_current_winner(current_idx)=true;

starts=table(start_index,family,start_id,start_seed,valid,current_eligible, ...
    objective_pool,objective_rank,iterations,converged,labels_changed_first,trim_changed_first, ...
    sign_score,rv.cv_score,rv.geom_score,within_family,all_stability,family_balanced, ...
    truth_ari,trim_precision,trim_recall,trim_f1,W_raw,W_ratio,R2,explained, ...
    is_current_winner,is_family_winner, ...
    'VariableNames',{'start_index','family','start_id','start_seed','valid','current_eligible', ...
    'objective','objective_rank','iterations','converged','labels_changed_first','trim_changed_first', ...
    'sign_score','rv_cv','rv_geom','within_family_stability','all_stability','family_balanced_stability', ...
    'truth_ari_inlier','trim_precision','trim_recall','trim_f1','W_raw','W_ratio','R2','explained', ...
    'is_current_winner','is_family_winner'});

% Family summaries and winners.
fam_rows=table();
for f=1:numel(names)
    idx=find(valid & family==names(f));
    if isempty(idx)
        row=table(names(f),0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            'VariableNames',{'family','n_valid','mean_ari','best_ari','p_ari_ge_090','p_ari_ge_095', ...
            'winner_index','winner_ari','winner_within_stability','winner_objective_rank'});
    else
        widx=family_winner_idx(f);
        a=truth_ari(idx);
        row=table(names(f),numel(idx),mean(a,'omitnan'),max(a,[],'omitnan'), ...
            mean(a>=.90,'omitnan'),mean(a>=.95,'omitnan'),double(widx),truth_ari(widx), ...
            family_winner_score(f),objective_rank(widx), ...
            'VariableNames',{'family','n_valid','mean_ari','best_ari','p_ari_ge_090','p_ari_ge_095', ...
            'winner_index','winner_ari','winner_within_stability','winner_objective_rank'});
    end
    fam_rows=[fam_rows;row]; %#ok<AGROW>
end

% Truth-blind rule winners for TFKM. TRKM remains frozen and is only audited.
rules=table();
if method=="REG_TFKM_HD"
    rules=[rules; rule_row("CURRENT_SIGN_ALL_STABLE",current_idx,current_native(current_idx),true,starts)]; %#ok<AGROW>

    rv_idx=choose_rv(pool,rv.cv_score,rv.geom_score,valid_idx);
    rules=[rules; rule_row("RV_CV_ALL",rv_idx,rv.cv_score(rv_idx),true,starts)]; %#ok<AGROW>

    fb_idx=choose_max_then_objective(family_balanced,valid,objective_pool);
    rules=[rules; rule_row("FB_STABLE_ALL",fb_idx,family_balanced(fb_idx),true,starts)]; %#ok<AGROW>

    all_idx=choose_max_then_objective(all_stability,valid,objective_pool);
    rules=[rules; rule_row("ALL_STABLE_NAIVE_DIAGNOSTIC",all_idx,all_stability(all_idx),true,starts)]; %#ok<AGROW>

    [hstab_idx,hstab_score]=choose_hier_stable(pool,S,family_winner_idx);
    rules=[rules; rule_row("HIER_STABLE4",hstab_idx,hstab_score,true,starts)]; %#ok<AGROW>

    hrv_idx=choose_rv(pool,rv.cv_score,rv.geom_score,family_winner_idx(isfinite(family_winner_idx)));
    rules=[rules; rule_row("HIER_RVCV4",hrv_idx,rv.cv_score(hrv_idx),true,starts)]; %#ok<AGROW>

    [~,mobj]=min(objective_pool(valid_idx)); minobj_idx=valid_idx(mobj);
    rules=[rules; rule_row("MINOBJ_ALL_DIAGNOSTIC",minobj_idx,-objective_pool(minobj_idx),true,starts)]; %#ok<AGROW>
else
    rules=[rules; rule_row("CURRENT_RV_CV",current_idx,rv.cv_score(current_idx),true,starts)]; %#ok<AGROW>
end

% Truth-open references: never selector challengers.
[~,m]=max(truth_ari(valid_idx)); oracle_idx=valid_idx(m);
rules=[rules; rule_row("ORACLE_BEST_ALL",oracle_idx,truth_ari(oracle_idx),false,starts)]; %#ok<AGROW>
sign_idx=find(valid & family=="SIGN_PCA_TKM");
if ~isempty(sign_idx)
    [~,m]=max(truth_ari(sign_idx)); oracle_sign=sign_idx(m);
    rules=[rules; rule_row("ORACLE_BEST_SIGN",oracle_sign,truth_ari(oracle_sign),false,starts)]; %#ok<AGROW>
end

% Pool-level truth-open quantities repeated later in every rules row by caller.
best_all_ari=max(truth_ari(valid_idx),[],'omitnan');
best_sign_ari=NaN;
if ~isempty(sign_idx), best_sign_ari=max(truth_ari(sign_idx),[],'omitnan'); end
non_sign_idx=find(valid & family~="SIGN_PCA_TKM");
best_non_sign_ari=NaN;
if ~isempty(non_sign_idx), best_non_sign_ari=max(truth_ari(non_sign_idx),[],'omitnan'); end

A=struct();
A.method=method; A.r=r; A.g=c; A.alpha=alpha; A.seed=seed;
A.pool=pool; A.starts=starts; A.family_summary=fam_rows; A.rules=rules;
A.current_idx=current_idx; A.family_winner_idx=family_winner_idx;
A.best_all_ari=best_all_ari; A.best_sign_ari=best_sign_ari; A.best_non_sign_ari=best_non_sign_ari;
A.oracle_best_idx=oracle_idx; A.oracle_best_family=family(oracle_idx);
A.pairwise_similarity=S;
end

function row=rule_row(rule,idx,score,truth_blind,starts)
row=table(string(rule),logical(truth_blind),double(idx),string(starts.family(idx)), ...
    starts.start_id(idx),starts.start_seed(idx),double(score),starts.truth_ari_inlier(idx), ...
    starts.objective(idx),starts.objective_rank(idx),starts.W_raw(idx),starts.W_ratio(idx), ...
    starts.R2(idx),starts.explained(idx),starts.rv_cv(idx),starts.rv_geom(idx), ...
    starts.within_family_stability(idx),starts.family_balanced_stability(idx), ...
    'VariableNames',{'rule','truth_blind','chosen_index','chosen_family','chosen_start_id','chosen_start_seed', ...
    'rule_score','selected_ari','objective','objective_rank','W_raw','W_ratio','R2','explained','rv_cv','rv_geom', ...
    'within_family_stability','family_balanced_stability'});
end

function S=pairwise_similarity(pool)
valid=logical(pool.valid(:)); N=numel(valid); S=nan(N); w=pool.trim_stability_weight;
if isempty(w)||~isfinite(w),w=.2;end
idx=find(valid);
for a=1:numel(idx)
    ia=idx(a); S(ia,ia)=0;
    for b=a+1:numel(idx)
        ib=idx(b); ta=logical(pool.trim(:,ia)); tb=logical(pool.trim(:,ib));
        in=~ta & ~tb;
        if sum(in)>1
            ar=adjusted_rand_index_v1(double(pool.labels(in,ia)),double(pool.labels(in,ib)));
            if ~isfinite(ar),ar=0;end
        else
            ar=0;
        end
        un=ta|tb; if any(un), jac=sum(ta&tb)/sum(un); else, jac=1; end
        S(ia,ib)=ar+w*jac; S(ib,ia)=S(ia,ib);
    end
end
end

function [wf,alls,fb]=stability_scores(pool,S,names)
valid=logical(pool.valid(:)); fam=string(pool.family(:)); N=numel(valid);
wf=nan(N,1); alls=nan(N,1); fb=nan(N,1);
vidx=find(valid);
if ~isempty(vidx), alls(vidx)=mean(S(vidx,vidx),2,'omitnan'); end
for f=1:numel(names)
    idx=find(valid & fam==names(f));
    if isempty(idx),continue;end
    % Exact frozen diagonal-zero convention within a family.
    wf(idx)=mean(S(idx,idx),2,'omitnan');
end
for j=vidx(:)'
    fm=nan(numel(names),1);
    for f=1:numel(names)
        idx=find(valid & fam==names(f)); idx(idx==j)=[];
        if ~isempty(idx),fm(f)=mean(S(j,idx),'omitnan');end
    end
    fb(j)=mean(fm,'omitnan');
end
end

function [widx,wscore]=choose_family_winners(pool,wf,names)
valid=logical(pool.valid(:)); fam=string(pool.family(:)); obj=double(pool.objective(:));
widx=nan(numel(names),1); wscore=nan(numel(names),1);
for f=1:numel(names)
    idx=find(valid & fam==names(f)); if isempty(idx),continue;end
    mx=max(wf(idx),[],'omitnan'); cand=idx(abs(wf(idx)-mx)<=1e-12);
    if numel(cand)>1,[~,m]=min(obj(cand)); cand=cand(m);end
    widx(f)=cand; wscore(f)=wf(cand);
end
end

function idx=choose_max_then_objective(score,valid,obj)
cand0=find(valid & isfinite(score));
if isempty(cand0), cand0=find(valid); [~,m]=min(obj(cand0)); idx=cand0(m); return; end
mx=max(score(cand0)); cand=cand0(abs(score(cand0)-mx)<=1e-12);
if numel(cand)>1,[~,m]=min(obj(cand)); idx=cand(m);else,idx=cand;end
end

function idx=choose_rv(pool,cv,gm,eligible)
obj=double(pool.objective(:)); eligible=eligible(:); eligible=eligible(isfinite(eligible));
eligible=eligible(logical(pool.valid(eligible)));
if isempty(eligible), error('basinAudit:NoEligibleRV','No eligible starts for RV selection.'); end
mx=max(cv(eligible),[],'omitnan'); cand=eligible(abs(cv(eligible)-mx)<=1e-12);
if isempty(cand),[~,m]=min(obj(eligible));idx=eligible(m);return;end
if numel(cand)>1
    a=max(cv(cand),[],'omitnan'); cand=cand(abs(cv(cand)-a)<=1e-12);
end
if numel(cand)>1
    b=max(gm(cand),[],'omitnan'); cand=cand(abs(gm(cand)-b)<=1e-12);
end
if numel(cand)>1,[~,m]=min(obj(cand));idx=cand(m);else,idx=cand;end
end

function [idx,score]=choose_hier_stable(pool,S,winners)
winners=winners(isfinite(winners)); winners=double(winners(:));
if isempty(winners),error('basinAudit:NoFamilyWinners','No valid family winners.');end
if numel(winners)==1,idx=winners(1);score=0;return;end
Q=S(winners,winners); qs=mean(Q,2,'omitnan'); mx=max(qs); cand=winners(abs(qs-mx)<=1e-12);
obj=double(pool.objective(:));
if numel(cand)>1,[~,m]=min(obj(cand));idx=cand(m);else,idx=cand;end
score=qs(winners==idx);
end

function chosen=legacy_candidate_local(pool,valid_idx,w)
obj=pool.objective(valid_idx); [~,ord]=sort(obj,'ascend'); ordered=valid_idx(ord);
top_n=min(numel(ordered),max(3,ceil(numel(ordered)/4))); top=ordered(1:top_n);
if numel(top)==1,chosen=top(1);return;end
S=zeros(numel(top));
for a=1:numel(top)
    for b=a+1:numel(top)
        ta=logical(pool.trim(:,top(a)));tb=logical(pool.trim(:,top(b)));in=~ta&~tb;
        if sum(in)>1
            ca=adjusted_rand_index_v1(double(pool.labels(in,top(a))),double(pool.labels(in,top(b))));
            if ~isfinite(ca),ca=0;end
        else,ca=0;end
        un=ta|tb;if any(un),tj=sum(ta&tb)/sum(un);else,tj=1;end
        S(a,b)=ca+w*tj;S(b,a)=S(a,b);
    end
end
[~,m]=max(mean(S,2));chosen=top(m);
end
