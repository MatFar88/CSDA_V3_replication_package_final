function sel=select_rg_frozen_v1(P,method)
%SELECT_RG_FROZEN_V1 Frozen external (r,g) selectors.
% TFKM: RANK_EN. TRKM: EXPL_NATIVE_010.
% Truth-blind: uses only fitted candidate diagnostics and method identity.
method=string(method); P=sortrows(P,{'r','g'});
assert(all(string(P.method)==method),'Mixed methods in RG path.');
E=nan(height(P),1); N=nan(height(P),1); K=double(P.r)+double(P.g);
if method=="REG_TFKM_HD"
 E=double(P.W_ratio).*double(P.R2); N=double(P.selector_score); rule="RANK_EN";
elseif method=="REG_TRKM_HD"
 E=double(P.R2); N=double(P.selector_cv_score); rule="EXPL_NATIVE_010";
else
 error('jointRG:Method','Unsupported method %s.',method);
end
N(~isfinite(N))=-Inf;
valid=logical(P.valid) & isfinite(E); idx=find(valid); assert(~isempty(idx),'No valid RG candidates.');
if rule=="RANK_EN"
 re=rank_desc(E,idx,height(P)); rn=rank_desc(N,idx,height(P)); score=re+rn;
 v=score(idx); mn=min(v); cand=idx(abs(v-mn)<=1e-12); k=parsimonious(P,K,cand);
else
 mx=max(E(idx)); keep=idx(E(idx)>=mx-.01-1e-12);
 vals=N(keep); mnative=max(vals);
 if ~isfinite(mnative), k=parsimonious(P,K,keep);
 else, cand=keep(abs(vals-mnative)<=1e-12); k=parsimonious(P,K,cand); end
end
sel=struct('r',double(P.r(k)),'g',double(P.g(k)),'row_index',k,'rule',rule, ...
 'explained',E(k),'native_score',N(k),'complexity',K(k));
end

function k=parsimonious(P,K,cand)
A=[K(cand) double(P.g(cand)) double(P.r(cand))]; [~,o]=sortrows(A,[1 2 3]); k=cand(o(1));
end
function r=rank_desc(x,idx,N)
r=nan(N,1); good=idx(isfinite(x(idx)));
if isempty(good),r(idx)=numel(idx);return;end
[~,o]=sort(x(good),'descend'); r(good(o))=(1:numel(good))'; missing=setdiff(idx,good); r(missing)=numel(idx)+1;
end
