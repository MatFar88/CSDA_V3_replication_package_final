function ari=adjusted_rand_index_v1(a,b)
a=a(:); b=b(:); [~,~,a]=unique(a); [~,~,b]=unique(b); N=accumarray([a b],1); n=sum(N(:)); c2=@(x) x.*(x-1)/2; s=sum(c2(N(:))); sa=sum(c2(sum(N,2))); sb=sum(c2(sum(N,1))); e=sa*sb/max(c2(n),eps); ari=(s-e)/max(.5*(sa+sb)-e,eps);
end
