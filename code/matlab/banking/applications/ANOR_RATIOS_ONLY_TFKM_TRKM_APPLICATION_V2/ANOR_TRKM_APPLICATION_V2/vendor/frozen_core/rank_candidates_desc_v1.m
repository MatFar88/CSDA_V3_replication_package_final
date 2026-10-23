function T=rank_candidates_desc_v1(T,score_name,nkeep)
%RANK_CANDIDATES_DESC_V1 Descending score, deterministic parsimonious ties.
s=double(T.(score_name));q=isfinite(s);T=T(q,:);s=s(q);if isempty(T),return;end
K=double(T.r)+double(T.g);A=[-s K double(T.g) double(T.r) double(T.alpha)];[~,o]=sortrows(A,[1 2 3 4 5]);T=T(o,:);
if nargin>=3&&height(T)>nkeep,T=T(1:nkeep,:);end
end
