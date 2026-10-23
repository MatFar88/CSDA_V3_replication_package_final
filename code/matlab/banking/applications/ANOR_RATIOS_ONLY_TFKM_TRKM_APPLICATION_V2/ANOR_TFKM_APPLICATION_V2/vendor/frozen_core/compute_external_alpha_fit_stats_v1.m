function stats = compute_external_alpha_fit_stats_v1(X,fit,method)
%COMPUTE_EXTERNAL_ALPHA_FIT_STATS Truth-blind scalar diagnostics for alpha path.
% No true labels/outlier indicators are used here.
assert(fit.valid,'Cannot score an invalid fit.');
method=string(method);
base=compute_external_rg_fit_stats_v1(X,fit,method);

in=~logical(fit.trim(:)); trim=~in;
labels=double(fit.labels(:)); U=fit.U;
mu=mean(X(in,:),1); Y=(X-mu)*U;
g=max(labels); r=size(U,2); centers=zeros(g,r);
for k=1:g
 q=in & labels==k;
 if any(q)
   centers(k,:)=mean(Y(q,:),1);
 else
   q=labels==k;
   if any(q), centers(k,:)=mean(Y(q,:),1); end
 end
end
resid=sum((Y-centers(labels,:)).^2,2);
keep_resid=resid(in);
scale=max(median(keep_resid),eps);
if any(trim)
 max_keep=max(keep_resid);
 min_trim=min(resid(trim));
 gap_ratio=(min_trim+eps(scale))/(max_keep+eps(scale));
 gap_log10=log10(max(gap_ratio,realmin));
 medk=median(keep_resid);
 madk=median(abs(keep_resid-medk));
 gap_robust=(min_trim-max_keep)/max(1.4826*madk,eps(scale));
else
 max_keep=max(keep_resid); min_trim=NaN;
 gap_ratio=1; gap_log10=0; gap_robust=0;
end
if method=="REG_TFKM_HD"
 explained=base.W_ratio*base.R2;
 native=field_or(fit,'selector_score',NaN);
else
 explained=base.R2;
 native=field_or(fit,'selector_cv_score',NaN);
end
stats=base;
stats.explained=explained;
stats.native_score=native;
stats.objective_per_keep=fit.objective/max(sum(in),1);
stats.trim_rate=mean(trim);
stats.trim_cut_max_keep=max_keep;
stats.trim_cut_min_trim=min_trim;
stats.trim_gap_ratio=gap_ratio;
stats.trim_gap_log10=gap_log10;
stats.trim_gap_robust=gap_robust;
end
function v=field_or(s,name,d)
if isfield(s,name)&&~isempty(s.(name)),v=s.(name);else,v=d;end
end
