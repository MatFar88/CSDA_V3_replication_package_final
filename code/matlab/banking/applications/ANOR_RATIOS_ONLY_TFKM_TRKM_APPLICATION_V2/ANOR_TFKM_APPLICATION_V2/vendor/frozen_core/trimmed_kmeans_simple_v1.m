function out=trimmed_kmeans_simple_v1(Y,c,alpha,seed,maxit)
rng(seed,'twister'); n=size(Y,1); g=kmeans(Y,c,'Replicates',1,'MaxIter',100,'Start','plus'); trim=false(n,1);
for it=1:maxit
 C=zeros(c,size(Y,2)); for k=1:c, C(k,:)=mean(Y(g==k & ~trim,:),1,'omitnan'); if any(isnan(C(k,:))), C(k,:)=Y(randi(n),:); end, end
 D=zeros(n,c); for k=1:c, Z=Y-C(k,:); D(:,k)=sum(Z.^2,2); end
 [d,gnew]=min(D,[],2); m=round(alpha*n); trimnew=false(n,1); if m>0, [~,ix]=maxk(d,m); trimnew(ix)=true; end
 if isequal(gnew,g)&&isequal(trimnew,trim), break; end; g=gnew; trim=trimnew;
end
out=struct('labels',g,'trim',trim);
end
