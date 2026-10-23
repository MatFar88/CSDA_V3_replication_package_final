function init=make_robust_initialization_v1(X,c,r,alpha,fam,seed)
rng(seed,'twister'); Xc=X-median(X,1); s=median(abs(Xc),1)+eps; Xr=Xc./s;
switch fam
 case "TKM_ORIGINAL", init=trimmed_kmeans_simple_v1(Xr,c,alpha,seed,20);
 case "SIGN_PCA_TKM", Z=Xr./max(vecnorm(Xr,2,2),eps); [~,~,V]=svd(Z,'econ'); init=trimmed_kmeans_simple_v1(Xr*V(:,1:r),c,alpha,seed,20);
 case "WINSOR_PCA_TKM", Z=max(min(Xr,4),-4); [~,~,V]=svd(Z,'econ'); init=trimmed_kmeans_simple_v1(Z*V(:,1:r),c,alpha,seed,20);
 case "RP_TKM", R=randn(size(X,2),max(r,5))/sqrt(max(r,5)); init=trimmed_kmeans_simple_v1(Xr*R,c,alpha,seed,20);
 case "PCA_KMEANS", [~,~,V]=svd(Xc,'econ'); g=kmeans(Xc*V(:,1:r),c,'Replicates',1,'Start','plus'); init=struct('labels',g,'trim',false(size(g)));
 otherwise, g=kmeans(Xr,c,'Replicates',1,'Start','plus'); init=struct('labels',g,'trim',false(size(g)));
end
end
