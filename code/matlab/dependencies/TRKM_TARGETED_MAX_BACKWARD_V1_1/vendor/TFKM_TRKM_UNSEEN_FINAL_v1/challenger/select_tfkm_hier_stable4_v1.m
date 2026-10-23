function q = select_tfkm_hier_stable4_v1(pool)
%SELECT_TFKM_HIER_STABLE4_V1 Truth-blind hierarchical TFKM arbitration.
%
% Frozen challenger definition (chosen after the fixed-triplet basin audit):
%   1. Within each initialization family, choose the most stable terminal fit.
%      Stability is the same TFKM similarity used by SIGN_ALL_STABLE:
%          ARI on mutually retained observations + w * Jaccard(trim sets),
%      with zero diagonal and ordinary row mean. Exact ties use lower objective.
%   2. Compare the available family winners with the same similarity.
%      Choose the most stable family representative; exact ties use lower objective.
%
% The unequal start allocation (50/35/10/5) therefore controls exploration only.
% Each family contributes at most one representative to the final competition.
% No truth labels, contamination indicators, oracle subspaces, or performance
% outcomes are accessed.

valid = logical(pool.valid(:));
fam = string(pool.family(:));
obj = double(pool.objective(:));
N = numel(valid);

family_names = ["WINSOR_PCA_TKM";"SIGN_PCA_TKM";"TKM_ORIGINAL";"RP_TKM"];
within_score = nan(N,1);
family_winner_idx = nan(numel(family_names),1);
family_winner_score = nan(numel(family_names),1);

for f = 1:numel(family_names)
    idx = find(valid & fam==family_names(f));
    if isempty(idx), continue; end
    s = within_set_stability(pool,idx);
    within_score(idx) = s;
    mx = max(s,[],'omitnan');
    cand = idx(abs(s-mx)<=1e-12);
    if numel(cand)>1
        [~,m] = min(obj(cand));
        cand = cand(m);
    end
    family_winner_idx(f) = cand;
    family_winner_score(f) = within_score(cand);
end

winners = family_winner_idx(isfinite(family_winner_idx));
winners = double(winners(:));
if isempty(winners)
    error('hierStable4:NoValidStarts','No valid starts in any initialization family.');
end

final_score = nan(N,1);
if numel(winners)==1
    chosen = winners(1);
    final_score(chosen) = within_score(chosen);
else
    S = similarity_matrix(pool,winners);
    s4 = mean(S,2,'omitnan'); % zero diagonal, same convention as SIGN_ALL_STABLE
    final_score(winners) = s4;
    mx = max(s4,[],'omitnan');
    cand = winners(abs(s4-mx)<=1e-12);
    if numel(cand)>1
        [~,m] = min(obj(cand));
        chosen = cand(m);
    else
        chosen = cand;
    end
end

q = struct();
q.chosen = chosen;
q.score = final_score;
q.within_family_score = within_score;
q.family_names = family_names;
q.family_winner_idx = family_winner_idx;
q.family_winner_score = family_winner_score;
q.rule = "HIER_STABLE4";
q.fallback_legacy = false;
end

function s = within_set_stability(pool,keep)
% Exact pairwise similarity and diagonal-zero normalization used by TFKM.
m = numel(keep);
if m==1
    s = 0;
    return;
end
S = similarity_matrix(pool,keep);
s = mean(S,2,'omitnan');
end

function S = similarity_matrix(pool,keep)
keep = double(keep(:));
m = numel(keep);
S = zeros(m);
w = pool.trim_stability_weight;
if isempty(w) || ~isfinite(w), w = 0.2; end
for a=1:m
    ia=keep(a);
    for b=a+1:m
        ib=keep(b);
        ta=logical(pool.trim(:,ia));
        tb=logical(pool.trim(:,ib));
        retained = ~ta & ~tb;
        if sum(retained)>1
            ar=adjusted_rand_index_v1(double(pool.labels(retained,ia)),double(pool.labels(retained,ib)));
            if ~isfinite(ar), ar=0; end
        else
            ar=0;
        end
        un=ta|tb;
        if any(un), jac=sum(ta&tb)/sum(un); else, jac=1; end
        S(a,b)=ar+w*jac;
        S(b,a)=S(a,b);
    end
end
end
