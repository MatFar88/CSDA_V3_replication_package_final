function fit = run_hd_multistart_top_hier_stable4_v1(X,r,c,alpha,options,seed)
%RUN_HD_MULTISTART_TOP_HIER_STABLE4_V1 TFKM challenger candidate fit.
%
% Candidate fitting, initialization allocation, deterministic seeds, ridge,
% ALS core, validity contract, and objective are identical to frozen v11.
% The ONLY methodological change is multistart arbitration:
% SIGN_ALL_STABLE -> HIER_STABLE4.

method = "REG_TFKM_HD";
alloc = options.hd_tfkm_allocation;
coremethod = "REG_TFKM";
selector_rule = "HIER_STABLE4";
tauW = options.tauW;
tauT = options.tauT;
names = string(options.hd_init_names);
N = sum(alloc);
candidates = cell(N,1);
valid = false(N,1);
z = 0;

for j=1:numel(names)
    for h=1:alloc(j)
        z=z+1;
        start_seed = seed + 1000*j + h;
        try
            rng(start_seed,'twister');
            init = make_robust_initialization_v1(X,c,r,alpha,names(j),start_seed);
            q = regularized_trimmed_fit_corrected_v1( ...
                X,init,r,c,alpha,coremethod,options.max_iterations,tauW,tauT);
            q.init_family = names(j);
            q.start_id = h;
            q.start_seed = start_seed;
            q.valid = all(isfinite(q.labels)) && isfinite(q.objective) && size(q.U,2)==r;
            q.message = "";
            candidates{z}=q;
            valid(z)=q.valid;
        catch ME
            candidates{z}=struct('valid',false,'message', ...
                string(ME.identifier)+": "+string(ME.message), ...
                'init_family',names(j),'start_id',h,'start_seed',start_seed);
        end
    end
end

pool_id = coremethod + "|" + string(seed) + "|" + strjoin(names + "=" + string(alloc),",");
valid_idx=find(valid);
if isempty(valid_idx)
    fit=invalid_fit(size(X,1),method,N,tauW,tauT,pool_id,selector_rule,candidates);
    return;
end

labels_pool=zeros(size(X,1),N,'uint16');
trim_pool=false(size(X,1),N);
objective_pool=nan(N,1);
family=strings(N,1);
start_id=nan(N,1);
start_seed=nan(N,1);
for i=1:N
    q=candidates{i};
    family(i)=string(q.init_family);
    start_id(i)=q.start_id;
    start_seed(i)=q.start_seed;
    if ~valid(i),continue;end
    labels_pool(:,i)=uint16(q.labels(:));
    trim_pool(:,i)=logical(q.trim(:));
    objective_pool(i)=q.objective;
end
[~,ord]=sort(objective_pool(valid_idx),'ascend');
ordered=valid_idx(ord);
objective_rank=nan(N,1);
objective_rank(ordered)=(1:numel(ordered))';

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

sel=select_tfkm_hier_stable4_v1(pool);
chosen=sel.chosen;
fit=candidates{chosen};
fit.selector=selector_rule;
fit.N_starts=N;
fit.valid_starts=numel(valid_idx);
fit.penalty_mode="RIDGE_PENALTY";
fit.tauW=tauW;
fit.tauT=tauT;
fit.valid=true;
fit.message="";
fit.init_pool_id=pool_id;
fit.selected_init_family=string(fit.init_family);
fit.selected_start_id=fit.start_id;
fit.selected_start_seed=fit.start_seed;
fit.selected_objective_rank=objective_rank(chosen);
fit.selector_score=sel.score(chosen);
fit.selector_cv_score=NaN;
fit.hier_family_names=sel.family_names;
fit.hier_family_winner_idx=sel.family_winner_idx;
fit.hier_family_winner_score=sel.family_winner_score;
fit.hier_within_family_score=sel.within_family_score(chosen);
end

function fit=invalid_fit(n,method,N,tauW,tauT,pool_id,selector_rule,candidates)
messages=strings(numel(candidates),1);
for j=1:numel(candidates)
    if isstruct(candidates{j}) && isfield(candidates{j},'message')
        messages(j)=string(candidates{j}.message);
    end
end
messages=unique(messages(strlength(messages)>0),'stable');
if isempty(messages),msg="No valid multistart candidate.";else,msg="No valid multistart candidate. First error: "+messages(1);end
fit=struct('labels',ones(n,1),'trim',false(n,1),'U',[], ...
    'objective',NaN,'iterations',NaN,'valid',false,'message',msg, ...
    'selector',selector_rule,'N_starts',N,'valid_starts',0, ...
    'penalty_mode',"RIDGE_PENALTY",'tauW',tauW,'tauT',tauT, ...
    'method',method,'init_pool_id',pool_id,'selected_init_family',"", ...
    'selected_start_id',NaN,'selected_start_seed',NaN, ...
    'selected_objective_rank',NaN,'selector_score',NaN,'selector_cv_score',NaN);
end
