function sel=apply_general_alpha_rule_v2(P,R)
%APPLY_GENERAL_ALPHA_RULE_V2 General truth-blind alpha-path selector.
% NEVER reads alpha0, ARI, contamination labels, role, truth masks, scenario
% difficulty, or a privileged alpha value.
P=sortrows(P,'alpha');
A=double(P.alpha(:)); valid=logical(P.valid) & isfinite(A);
assert(sum(valid)>=3,'Need at least 3 valid alpha candidates.');
idx=find(valid); ord=idx(:); [~,oo]=sort(A(ord)); ord=ord(oo);
W=finitevec(double(P.W_raw)); O=finitevec(double(P.objective_per_keep));
N=finitevec(double(P.native_score)); E=finitevec(double(P.explained));

% All fit components are path-normalized losses in [0,1].
lw=norm_low_loss(logpos(W),idx); lo=norm_low_loss(logpos(O),idx);
ln=norm_high_loss(N,idx); le=norm_high_loss(E,idx);
ar=zeros(size(A)); amin=min(A(idx)); amax=max(A(idx));
if amax>amin, ar(idx)=(A(idx)-amin)/(amax-amin); end

fam=string(R.family); l=double(R.p1); s=double(R.p2);
switch fam
 case "IC_OBJ"
  loss=lo + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "IC_W"
  loss=lw + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "IC_ON"
  loss=.5*(lo+ln) + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "IC_WN"
  loss=.5*(lw+ln) + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "IC_OWN"
  loss=(lo+lw+ln)/3 + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "IC_OWNE"
  loss=(lo+lw+ln+le)/4 + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "IC_ROB"
  M=[lo lw ln le]; rob=median(M,2); loss=rob + l*ar; k=argmin_then_alpha(loss,idx,A);
 case "GAIN_OBJ"
  k=gain_rule(A,lo,idx,l);
 case "GAIN_OWN"
  base=(lo+lw+ln)/3; k=gain_rule(A,base,idx,l);
 case "STOP_OBJ"
  k=stop_rule(A,lo,idx,l);
 case "STOP_W"
  k=stop_rule(A,lw,idx,l);
 case "STOP_OWN"
  base=(lo+lw+ln)/3; k=stop_rule(A,base,idx,l);
 case "IC_STAB"
  % Penalized fit loss, but reject unstable candidates by adding a smooth
  % excess penalty rather than a hard alpha cutoff.
  base=.5*(lo+lw); stab=.5*(ln+le);
  excess=max(stab-s,0); loss=base + l*ar + 2*excess;
  k=argmin_then_alpha(loss,idx,A);
 otherwise
  error('generalAlpha:Rule','Unknown family %s.',fam);
end
if isempty(k)||~isfinite(k)||k<1, k=min_alpha_idx(A,idx); end
sel=struct('alpha',A(k),'row_index',k,'rule',string(R.rule),'family',fam);
end

function k=gain_rule(A,loss,idx,lambda)
% Maximize improvement from the least-trimmed fit minus trimming cost.
ord=idx(:); [~,o]=sort(A(ord)); ord=ord(o);
a=A(ord); base=loss(ord(1)); gain=base-loss(ord);
ar=(a-a(1))/max(a(end)-a(1),eps); utility=gain-lambda*ar;
mx=max(utility); c=ord(abs(utility-mx)<=1e-12); k=min_alpha_idx(A,c);
end

function k=stop_rule(A,loss,idx,thr)
% Select alpha after the last materially beneficial adjacent step.
ord=idx(:); [~,o]=sort(A(ord)); ord=ord(o); m=numel(ord);
if m<2, k=ord(1); return; end
a=A(ord); y=loss(ord); da=max(diff(a),eps); imp=max(0,-diff(y)./da);
mx=max(imp); if ~isfinite(mx)||mx<=eps, k=ord(1); return; end
sig=find(imp>=thr*mx-1e-12);
if isempty(sig), k=ord(1); else, k=ord(max(sig)+1); end
end

function x=finitevec(x)
x=x(:); bad=~isfinite(x); if all(bad),x(:)=0;elseif any(bad),x(bad)=median(x(~bad));end
end
function y=logpos(x)
x=x(:); pos=x(isfinite(x)&x>0); if isempty(pos),s=eps;else,s=max(min(pos)*1e-6,eps);end
y=log(max(x,s));
end
function loss=norm_low_loss(x,idx)
loss=ones(size(x)); v=x(idx); lo=min(v); hi=max(v); rg=hi-lo;
if rg<=eps(max(abs(v))+1), loss(idx)=0; else, loss(idx)=(v-lo)/rg; end
end
function loss=norm_high_loss(x,idx)
loss=ones(size(x)); v=x(idx); lo=min(v); hi=max(v); rg=hi-lo;
if rg<=eps(max(abs(v))+1), loss(idx)=0; else, loss(idx)=(hi-v)/rg; end
end
function k=argmin_then_alpha(x,idx,A)
v=x(idx); m=min(v); c=idx(abs(v-m)<=1e-12); k=min_alpha_idx(A,c);
end
function k=min_alpha_idx(A,c)
[~,j]=min(A(c)); k=c(j);
end
