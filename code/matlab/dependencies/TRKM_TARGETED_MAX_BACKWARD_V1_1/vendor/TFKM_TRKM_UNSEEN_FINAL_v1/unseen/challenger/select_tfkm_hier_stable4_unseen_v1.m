function q = select_tfkm_hier_stable4_unseen_v1(pool)
%SELECT_TFKM_HIER_STABLE4_UNSEEN_V1 Exact frozen HIER4 + component bookkeeping.
% Operational choice is IDENTICAL to HIER_STABLE4 development definition.
% Extra outputs only decompose the chosen cross-family score into ARI and J.

valid=logical(pool.valid(:));
fam=string(pool.family(:));
obj=double(pool.objective(:));
N=numel(valid);
family_names=["WINSOR_PCA_TKM";"SIGN_PCA_TKM";"TKM_ORIGINAL";"RP_TKM"];
within_score=nan(N,1);
family_winner_idx=nan(numel(family_names),1);
family_winner_score=nan(numel(family_names),1);

for f=1:numel(family_names)
    idx=find(valid & fam==family_names(f));
    if isempty(idx), continue; end
    [~,~,H]=component_matrices(pool,idx);
    s=mean(H,2,'omitnan');
    within_score(idx)=s;
    mx=max(s,[],'omitnan');
    cand=idx(abs(s-mx)<=1e-12);
    if numel(cand)>1
        [~,m]=min(obj(cand)); cand=cand(m);
    end
    family_winner_idx(f)=cand;
    family_winner_score(f)=within_score(cand);
end

winners=double(family_winner_idx(isfinite(family_winner_idx)));
winners=winners(:);
assert(~isempty(winners),'hierStable4:NoValidStarts');

final_score=nan(N,1);
cross_ari=nan(N,1);
cross_j=nan(N,1);

if numel(winners)==1
    chosen=winners(1);
    final_score(chosen)=within_score(chosen);
    cross_ari(chosen)=0;
    cross_j(chosen)=0;
else
    [A,J,H]=component_matrices(pool,winners);
    sH=mean(H,2,'omitnan');
    sA=mean(A,2,'omitnan');
    sJ=mean(J,2,'omitnan');
    final_score(winners)=sH;
    cross_ari(winners)=sA;
    cross_j(winners)=sJ;
    mx=max(sH,[],'omitnan');
    cand=winners(abs(sH-mx)<=1e-12);
    if numel(cand)>1
        [~,m]=min(obj(cand)); chosen=cand(m);
    else
        chosen=cand;
    end
end

q=struct();
q.chosen=chosen;
q.score=final_score;
q.within_family_score=within_score;
q.family_names=family_names;
q.family_winner_idx=family_winner_idx;
q.family_winner_score=family_winner_score;
q.cross_ari_score=cross_ari;
q.cross_j_score=cross_j;
q.rule="HIER_STABLE4";
end

function [A,J,H]=component_matrices(pool,keep)
keep=double(keep(:)); m=numel(keep);
A=zeros(m); J=zeros(m);
w=pool.trim_stability_weight;
if isempty(w)||~isfinite(w), w=.2; end
for a=1:m
    ia=keep(a);
    for b=a+1:m
        ib=keep(b);
        ta=logical(pool.trim(:,ia)); tb=logical(pool.trim(:,ib));
        retained=~ta & ~tb;
        if sum(retained)>1
            ar=adjusted_rand_index_v1(double(pool.labels(retained,ia)),double(pool.labels(retained,ib)));
            if ~isfinite(ar), ar=0; end
        else
            ar=0;
        end
        un=ta|tb;
        if any(un), jac=sum(ta&tb)/sum(un); else, jac=1; end
        A(a,b)=ar; A(b,a)=ar;
        J(a,b)=jac; J(b,a)=jac;
    end
end
H=A+w*J;
end
