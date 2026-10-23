function sel=select_alpha_from_path_v1(P,rule)
%SELECT_ALPHA_FROM_PATH_V1 Truth-blind alpha selection from a completed path.
% Truth/alpha0/ARI/contamination labels are never used by this function.
rule=string(rule);
valid=logical(P.valid) & isfinite(P.W_raw) & isfinite(P.explained);
assert(any(valid),'No valid alpha candidate.');
idx=find(valid); A=double(P.alpha); E=double(P.explained); N=double(P.native_score); G=double(P.trim_gap_log10); W=double(P.W_raw);
N(~isfinite(N))=-Inf; G(~isfinite(G))=-Inf;

switch rule
 case {"SAFEDEF_BASEZERO_1","SAFEDEF_BASEZERO_2","SAFEDEF_BASEZERO_3","SAFEDEF_TRI_1","SAFEDEF_TRI_2","SAFEDEF_TRI_3","SAFEDEF_TRI_4","ALWAYS_SAFE_MIN"}
   ss=select_alpha_safedefault_rule_v4(P,rule); k=ss.row_index;
 case {"LEGACY_ALPHA_005","LEGACY_ALPHA_010","LEGACY_ALPHA_020","LEGACY_ALPHA_050"}
   tol=rule_tol(rule); k=legacy_alpha(P,idx,tol);
 case "MIN_WRAW"
   mn=min(W(idx)); cand=idx(abs(W(idx)-mn)<=1e-12); k=smallest_alpha(P,cand);
 case "EXPL_MAX"
   mx=max(E(idx)); cand=idx(abs(E(idx)-mx)<=1e-12); k=smallest_alpha(P,cand);
 case {"EXPL_PLAT_005","EXPL_PLAT_010","EXPL_PLAT_020"}
   tol=rule_tol(rule); mx=max(E(idx)); cand=idx(E(idx)>=mx-tol); k=smallest_alpha(P,cand);
 case {"EXPL_NATIVE_005","EXPL_NATIVE_010","EXPL_NATIVE_020"}
   tol=rule_tol(rule); mx=max(E(idx)); keep=idx(E(idx)>=mx-tol); k=best_native_then_alpha(P,N,keep);
 case "GAP_MAX"
   mx=max(G(idx)); cand=idx(abs(G(idx)-mx)<=1e-12); k=smallest_alpha(P,cand);
 case {"RANK_EN_A05","RANK_EG_A05","RANK_ENG_A05","RANK_WNG_A05"}
   re=rank_desc(E,idx,height(P)); rn=rank_desc(N,idx,height(P)); rg=rank_desc(G,idx,height(P)); rw=rank_asc(W,idx,height(P)); ra=rank_asc(A,idx,height(P));
   if rule=="RANK_EN_A05"
      score=re+rn+0.5*ra;
   elseif rule=="RANK_EG_A05"
      score=re+rg+0.5*ra;
   elseif rule=="RANK_ENG_A05"
      score=re+rn+rg+0.5*ra;
   else
      score=rw+rn+rg+0.5*ra;
   end
   v=score(idx); mn=min(v); cand=idx(abs(v-mn)<=1e-12); k=smallest_alpha(P,cand);
 case "FIRSTDROP_R30_SUFFIXN010"
   k=firstdrop_suffix_native(A,E,N,idx,.30,.01);
   if isempty(k)
      mx=max(E(idx)); keep=idx(E(idx)>=mx-.02); k=best_native_then_alpha(P,N,keep);
   end
 otherwise
   error('selectAlpha:Rule','Unknown alpha rule %s.',rule);
end
sel=struct('alpha',A(k),'row_index',k,'rule',rule,'explained',E(k),'native_score',N(k),'gap_log10',G(k),'W_raw',W(k));
end

function k=firstdrop_suffix_native(A,E,N,idx,thr,native_tol)
k=[]; ord=idx(:); [~,o]=sort(A(ord),'ascend'); ord=ord(o);
if numel(ord)<2,return;end
Er=E(ord); D=Er(1:end-1)-Er(2:end); pr=max(Er)-min(Er);
if ~(isfinite(pr)&&pr>0),return;end
j=find(D>=thr*pr & D>0,1,'first');
if isempty(j),return;end
pos=j+1; suf=ord(pos:end); vals=N(suf); nm=max(vals);
if ~isfinite(nm),k=suf(1);return;end
cand=suf(vals>=nm-native_tol); [~,jj]=min(A(cand)); k=cand(jj);
end
function k=legacy_alpha(P,idx,tol)
W=double(P.W_raw); method=string(P.method(1)); mn=min(W(idx));
if method=="REG_TFKM_HD", wmax=.05; else, wmax=Inf; end
eligible=idx(W(idx)<=wmax & W(idx)-mn<=tol);
if isempty(eligible), cand=idx(abs(W(idx)-mn)<=1e-12); k=smallest_alpha(P,cand); else, k=smallest_alpha(P,eligible); end
end
function k=best_native_then_alpha(P,N,keep)
vals=N(keep); mx=max(vals); if ~isfinite(mx),k=smallest_alpha(P,keep);return;end
cand=keep(abs(vals-mx)<=1e-12);k=smallest_alpha(P,cand);
end
function k=smallest_alpha(P,cand)
A=double(P.alpha(cand));mn=min(A);cand=cand(abs(A-mn)<=1e-12);k=cand(1);
end
function t=rule_tol(rule)
if contains(rule,"005"),t=.005;elseif contains(rule,"010"),t=.01;elseif contains(rule,"020"),t=.02;elseif contains(rule,"050"),t=.05;else,error('No tolerance in %s.',rule);end
end
function r=rank_desc(x,idx,Ntot)
r=nan(Ntot,1);good=idx(isfinite(x(idx)));if isempty(good),r(idx)=numel(idx);return;end
[~,ord]=sort(x(good),'descend');rr=zeros(numel(good),1);rr(ord)=1:numel(good);r(good)=rr;r(idx(~isfinite(x(idx))))=numel(good)+1;
end
function r=rank_asc(x,idx,Ntot)
r=nan(Ntot,1);good=idx(isfinite(x(idx)));if isempty(good),r(idx)=numel(idx);return;end
[~,ord]=sort(x(good),'ascend');rr=zeros(numel(good),1);rr(ord)=1:numel(good);r(good)=rr;r(idx(~isfinite(x(idx))))=numel(good)+1;
end
