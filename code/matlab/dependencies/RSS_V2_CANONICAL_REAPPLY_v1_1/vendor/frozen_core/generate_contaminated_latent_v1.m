function sim=generate_contaminated_latent_v1(cfg,seed)
rng(seed,'twister'); n=cfg.n; p=cfg.p; r=cfg.r_true; c=cfg.c_true;
if isfield(cfg,'separation'), separation=cfg.separation; elseif isfield(cfg,'delta'), separation=cfg.delta; else, error('Missing separation/delta field.'); end
if isfield(cfg,'spectral_ratio'), spectral_ratio=cfg.spectral_ratio; elseif isfield(cfg,'rho'), spectral_ratio=cfg.rho; else, error('Missing spectral_ratio/rho field.'); end
[U,~]=qr(randn(p,r),0); [V,~]=qr(randn(p,p-r),0); V=V-U*(U'*V); [V,~]=qr(V,0);
E=eye(c)-ones(c)/c; [Q,~,~]=svd(E,'econ'); C=Q(:,1:r); C=C*(separation/mean(pdist(C)));
g=repelem((1:c)',ceil(n/c)); g=g(1:n); g=g(randperm(n)); Z=C(g,:)+randn(n,r); W=randn(n,p-r)*sqrt(spectral_ratio); X=Z*U'+W*V';
out=false(n,1); m=round(cfg.alpha*n);
if m>0
 idx=randperm(n,m); out(idx)=true;
 if cfg.contamination=="leverage"
   X(idx,:)=X(idx,:)+10*randn(m,p).*repmat(std(X,0,1)+eps,m,1);
 else
   for ii=1:m
     k=g(idx(ii)); kk=mod(k,c)+1; z=.5*(C(k,:)+C(kk,:))+2*randn(1,r); X(idx(ii),:)=z*U'+randn(1,p-r)*sqrt(spectral_ratio)*V';
   end
 end
end
X=X-mean(X,1); sim=struct('X',X,'labels',g,'outlier',out,'U',U);
end
