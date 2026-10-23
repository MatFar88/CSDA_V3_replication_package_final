function views = build_robust_views_v3(X,r,c)
%BUILD_ROBUST_VIEWS_V3 Two deterministic robust validation views.
% These views are independent of the candidate fit and never use truth.
% They mirror the robust geometry already used by the initialization bank:
% SIGN-PCA and WINSOR-PCA.

Xc = X - median(X,1);
s = median(abs(Xc),1) + eps;
Xr = Xc ./ s;
n = size(X,1); p = size(X,2);
d = min([10,p,n-1,max([r,c-1])*2]);
d = max(1,d);

% SIGN-PCA directions; validation coordinates use standardized Xr projected
% onto directions learned from row-normalized data.
Zs = Xr ./ max(vecnorm(Xr,2,2),eps);
[~,~,Vs] = svd(Zs,'econ');
Vs = Vs(:,1:min(d,size(Vs,2)));
S = Xr * Vs;

% WINSOR-PCA view.
Zw = max(min(Xr,4),-4);
[~,~,Vw] = svd(Zw,'econ');
Vw = Vw(:,1:min(d,size(Vw,2)));
W = Zw * Vw;

% Standardize each view coordinate robustly so no single component dominates.
S = robust_colscale(S);
W = robust_colscale(W);
views = {S,W};
end

function Z=robust_colscale(Z)
Z = Z - median(Z,1);
sc = median(abs(Z),1) + eps;
Z = Z ./ sc;
end
