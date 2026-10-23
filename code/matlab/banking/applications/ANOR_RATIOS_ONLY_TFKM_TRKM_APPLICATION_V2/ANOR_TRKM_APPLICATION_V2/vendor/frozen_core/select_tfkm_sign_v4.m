function q = select_tfkm_sign_v4(pool,X,c,r,rule)
%SELECT_TFKM_SIGN_V4 Truth-blind TFKM-specific selector candidates.
% Restricts arbitration to SIGN_PCA_TKM starts. Objective/stability/CV are
% used without access to truth labels or outlier indicators.

rule=string(rule); valid=logical(pool.valid(:)); fam=string(pool.family(:));
idx=find(valid & fam=="SIGN_PCA_TKM");
if isempty(idx)
    chosen=pool.current_candidate_index;
    q=struct('chosen',chosen,'score',nan(numel(valid),1),'cv_score',nan(numel(valid),1), ...
        'rule',rule,'fallback_legacy',true);
    return;
end
obj=double(pool.objective(:)); N=numel(valid); score=nan(N,1); cv=nan(N,1);

switch rule
    case "SIGN_MINOBJ"
        [~,m]=min(obj(idx)); chosen=idx(m); score(chosen)=1;
    case {"SIGN_Q25_STABLE","SIGN_Q50_STABLE","SIGN_ALL_STABLE"}
        [~,ord]=sort(obj(idx),'ascend'); ordered=idx(ord);
        if rule=="SIGN_Q25_STABLE"
            keep_n=max(3,ceil(numel(ordered)/4));
        elseif rule=="SIGN_Q50_STABLE"
            keep_n=max(3,ceil(numel(ordered)/2));
        else
            keep_n=numel(ordered);
        end
        keep=ordered(1:keep_n);
        s=within_set_stability(pool,keep);
        score(keep)=s;
        mx=max(s); cand=keep(abs(s-mx)<=1e-12);
        if numel(cand)>1, [~,m]=min(obj(cand)); chosen=cand(m); else, chosen=cand; end
    case {"SIGN_RV_CV","SIGN_Q50_RV_CV","SIGN_OBJ_CV_RANK"}
        sc=score_candidates_robust_views_v3(pool,X,c,r); cv=sc.cv_score;
        if rule=="SIGN_RV_CV"
            keep=idx; score=cv;
        elseif rule=="SIGN_Q50_RV_CV"
            [~,ord]=sort(obj(idx),'ascend'); ordered=idx(ord); keep=ordered(1:max(3,ceil(numel(ordered)/2))); score(keep)=cv(keep);
        else
            keep=idx;
            ro=ascending_rank(obj,keep,N); rc=descending_rank(cv,keep,N);
            score(keep)=-(ro(keep)+rc(keep))/2;
        end
        mx=max(score(keep),[],'omitnan'); cand=keep(abs(score(keep)-mx)<=1e-12);
        if isempty(cand), [~,m]=min(obj(keep)); chosen=keep(m);
        elseif numel(cand)>1, [~,m]=min(obj(cand)); chosen=cand(m);
        else, chosen=cand; end
    otherwise
        error('selectorv4:Rule','Unknown TFKM rule %s.',rule);
end
q=struct('chosen',chosen,'score',score,'cv_score',cv,'rule',rule,'fallback_legacy',false);
end

function s=within_set_stability(pool,keep)
m=numel(keep); S=zeros(m); w=pool.trim_stability_weight;
if isempty(w) || ~isfinite(w), w=0.2; end
for a=1:m
    for b=a+1:m
        ia=keep(a); ib=keep(b); ta=logical(pool.trim(:,ia)); tb=logical(pool.trim(:,ib));
        in=~ta & ~tb;
        if sum(in)>1
            ca=adjusted_rand_index_v1(double(pool.labels(in,ia)),double(pool.labels(in,ib)));
            if ~isfinite(ca),ca=0;end
        else, ca=0; end
        un=ta|tb; if any(un),tj=sum(ta&tb)/sum(un);else,tj=1;end
        S(a,b)=ca+w*tj;S(b,a)=S(a,b);
    end
end
s=mean(S,2); % same diagonal-zero convention as legacy
end

function r=ascending_rank(x,keep,N)
r=nan(N,1); [~,ord]=sort(x(keep),'ascend'); r(keep(ord))=(1:numel(keep))';
end
function r=descending_rank(x,keep,N)
r=nan(N,1); [~,ord]=sort(x(keep),'descend'); r(keep(ord))=(1:numel(keep))';
end
