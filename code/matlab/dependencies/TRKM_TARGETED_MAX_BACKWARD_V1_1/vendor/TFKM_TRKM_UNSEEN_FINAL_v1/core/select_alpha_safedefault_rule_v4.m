function sel=select_alpha_safedefault_rule_v4(P,rule)
% Truth-blind default-safe alpha selector. Positive branch is alpha>=0.10.
% Rule never uses alpha0, contamination labels, ARI or truth masks.
P=sortrows(P,'alpha'); rule=string(rule);
valid=logical(P.valid)&isfinite(P.alpha)&isfinite(P.explained);
A=double(P.alpha); E=double(P.explained); N=double(P.native_score); G=double(P.trim_gap_log10);
N(~isfinite(N))=-Inf; G(~isfinite(G))=-Inf; idx=find(valid);
zero=find(valid & abs(A)<=1e-12,1,'first'); assert(~isempty(zero),'alpha=0 missing');
safe=idx(A(idx)>=.10-1e-12); assert(~isempty(safe),'alpha>=.10 missing');
ks=min_alpha_idx(A,safe); % exact-safe branch: smallest admissible positive alpha
if rule=="ALWAYS_SAFE_MIN"
    k=ks; why="ALWAYS_SAFE";
else
    ga=max(N(safe))-N(zero);
    dr=boundary_drop_rel(A,E,idx,safe);
    mg=max(G(safe));
    bz=base_rule_idx(P); baseZero=A(bz)<=1e-12;
    % Default is SAFE. ZERO is released only under strong clean-like evidence.
    switch rule
        case "SAFEDEF_NLOSS020"
            goZero=ga<=-.020;
        case "SAFEDEF_NLOSS010"
            goZero=ga<=-.010;
        case "SAFEDEF_NLOSS005"
            goZero=ga<=-.005;
        case "SAFEDEF_NLOSS0025"
            goZero=ga<=-.0025;
        case "SAFEDEF_NGAIN000"
            goZero=ga<=0;
        case "SAFEDEF_TRI_1"
            goZero=(ga<=.005)&&(dr<=.05)&&(mg<=.25);
        case "SAFEDEF_TRI_2"
            goZero=(ga<=.005)&&(dr<=.10)&&(mg<=.50);
        case "SAFEDEF_TRI_3"
            goZero=(ga<=.010)&&(dr<=.10)&&(mg<=.50);
        case "SAFEDEF_TRI_4"
            goZero=(ga<=.010)&&(dr<=.20)&&(mg<=1.00);
        case "SAFEDEF_BASEZERO_1"
            goZero=baseZero&&(ga<=.005)&&(dr<=.05)&&(mg<=.25);
        case "SAFEDEF_BASEZERO_2"
            goZero=baseZero&&(ga<=.010)&&(dr<=.10)&&(mg<=.50);
        case "SAFEDEF_BASEZERO_3"
            goZero=baseZero&&(ga<=.020)&&(dr<=.20)&&(mg<=1.00);
        case "SAFEDEF_ZEROALL_1"
            % Strongest clean release: zero must look better in native score,
            % with no marked explained transition and weak trim-gap evidence.
            goZero=(ga<=0)&&(dr<=.10)&&(mg<=.50);
        case "SAFEDEF_ZEROALL_2"
            goZero=(ga<=.0025)&&(dr<=.05)&&(mg<=.50);
        case "SAFEDEF_ZEROALL_3"
            goZero=(ga<=.005)&&(dr<=.05)&&(mg<=1.00);
        otherwise
            error('alphaSafeV4:Rule','Unknown rule %s.',rule);
    end
    if goZero, k=zero; why="ZERO_RELEASE"; else, k=ks; why="SAFE_DEFAULT"; end
end
sel=struct('alpha',A(k),'row_index',k,'rule',rule,'branch',why, ...
    'native_gain_safe_vs_zero',max(N(safe))-N(zero), ...
    'boundary_drop_rel',boundary_drop_rel(A,E,idx,safe),'max_safe_gap',max(G(safe)));
end
function d=boundary_drop_rel(A,E,idx,safe)
[~,j]=min(A(safe));ks=safe(j);pre=idx(A(idx)<A(ks)-1e-12);vals=E(idx);rg=max(vals)-min(vals);
if isempty(pre)||~isfinite(rg)||rg<=eps,d=0;return;end
d=max(0,max(E(pre))-E(ks))/rg;
end
function k=min_alpha_idx(A,cand)
mn=min(A(cand));c=cand(abs(A(cand)-mn)<=1e-12);k=c(1);
end
function k=base_rule_idx(P)
if string(P.method(1))=="REG_TFKM_HD",k=rank_en_a05(P);else,k=firstdrop_r30_suffixn010(P);end
end
function k=rank_en_a05(P)
valid=logical(P.valid)&isfinite(P.explained)&isfinite(P.alpha);idx=find(valid);A=double(P.alpha);E=double(P.explained);N=double(P.native_score);N(~isfinite(N))=-Inf;
re=rank_desc(E,idx,height(P));rn=rank_desc(N,idx,height(P));ra=rank_asc(A,idx,height(P));s=re+rn+.5*ra;v=s(idx);mn=min(v);cand=idx(abs(v-mn)<=1e-12);k=min_alpha_idx(A,cand);
end
function k=firstdrop_r30_suffixn010(P)
valid=logical(P.valid)&isfinite(P.explained)&isfinite(P.alpha);idx=find(valid);A=double(P.alpha);E=double(P.explained);N=double(P.native_score);N(~isfinite(N))=-Inf;
ord=idx(:);[~,o]=sort(A(ord),'ascend');ord=ord(o);Er=E(ord);D=Er(1:end-1)-Er(2:end);rg=max(Er)-min(Er);pos=[];
if isfinite(rg)&&rg>0,pos=find(D>=.30*rg & D>0,1,'first');end
if isempty(pos),mx=max(E(idx));keep=idx(E(idx)>=mx-.02);vals=N(keep);nm=max(vals);if isfinite(nm),cand=keep(vals>=nm-.01);else,cand=keep;end;k=min_alpha_idx(A,cand);return;end
suf=ord(pos+1:end);vals=N(suf);nm=max(vals);if ~isfinite(nm),k=suf(1);return;end;cand=suf(vals>=nm-.01);k=min_alpha_idx(A,cand);
end
function r=rank_desc(x,idx,Nn)
r=nan(Nn,1);good=idx(isfinite(x(idx)));[~,o]=sort(x(good),'descend');r(good(o))=(1:numel(good))';miss=setdiff(idx,good);r(miss)=numel(idx)+1;
end
function r=rank_asc(x,idx,Nn)
r=nan(Nn,1);good=idx(isfinite(x(idx)));[~,o]=sort(x(good),'ascend');r(good(o))=(1:numel(good))';miss=setdiff(idx,good);r(miss)=numel(idx)+1;
end
