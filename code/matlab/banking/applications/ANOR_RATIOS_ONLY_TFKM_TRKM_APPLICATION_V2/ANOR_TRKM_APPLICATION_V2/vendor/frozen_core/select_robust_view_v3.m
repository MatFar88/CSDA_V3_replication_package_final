function q = select_robust_view_v3(pool,X,c,r,rule)
%SELECT_ROBUST_VIEW_V3 Prespecified truth-blind robust-view selectors.
% Objective is used only as the last deterministic tie-breaker.

rule=string(rule); sc=score_candidates_robust_views_v3(pool,X,c,r);
valid=logical(pool.valid(:)); obj=double(pool.objective(:));
cv=sc.cv_score; gm=sc.geom_score;
N=numel(valid); score=nan(N,1);

switch rule
    case "RV_CV"
        score=cv;
    case "RV_GEOM"
        score=gm;
    case "RV_CV_GEOM"
        % Equal-rank aggregation avoids arbitrary scale mixing.
        rcv=descending_rank(cv,valid); rgm=descending_rank(gm,valid);
        score=-(rcv+rgm)/2;
    otherwise
        error('selectorv3:Rule','Unknown rule %s.',rule);
end

mx=max(score(valid),[],'omitnan'); cand=find(valid & abs(score-mx)<=1e-12);
if isempty(cand)
    cand=find(valid); [~,m]=min(obj(cand)); chosen=cand(m);
elseif numel(cand)==1
    chosen=cand;
else
    % Tie-break 1: higher CV, 2: higher geometry, 3: lower objective.
    a=max(cv(cand),[],'omitnan'); c2=cand(abs(cv(cand)-a)<=1e-12);
    if numel(c2)>1
        b=max(gm(c2),[],'omitnan'); c2=c2(abs(gm(c2)-b)<=1e-12);
    end
    if numel(c2)>1, [~,m]=min(obj(c2)); chosen=c2(m); else, chosen=c2; end
end
q=struct('chosen',chosen,'score',score,'cv_score',cv,'geom_score',gm,'rule',rule, ...
    'view_cv',sc.view_cv,'view_geom',sc.view_geom);
end

function r=descending_rank(x,valid)
r=nan(size(x)); idx=find(valid & isfinite(x));
if isempty(idx), return; end
[~,ord]=sort(x(idx),'descend'); r(idx(ord))=(1:numel(idx))';
end
