function Z=build_bqs_common_view_v1(X)
%BUILD_BQS_COMMON_VIEW_V1 Candidate-independent WINSOR-PCA score space.
% Truth blind. Fit candidates remain in original X; only BQS scoring uses Z.
Xc=X-median(X,1);s=median(abs(Xc),1)+eps;Xr=Xc./s;
Zw=max(min(Xr,4),-4);n=size(X,1);p=size(X,2);d=max(1,min([10,p,n-1]));
[~,~,V]=svd(Zw,'econ');V=V(:,1:min(d,size(V,2)));Z=Zw*V;
Z=Z-median(Z,1);sc=median(abs(Z),1)+eps;Z=Z./sc;
end
