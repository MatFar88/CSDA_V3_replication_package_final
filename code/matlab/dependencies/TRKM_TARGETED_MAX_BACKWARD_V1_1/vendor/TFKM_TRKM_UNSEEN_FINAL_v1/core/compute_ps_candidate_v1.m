function out=compute_ps_candidate_v1(X,r,g,alpha,method,options,data_seed)
%COMPUTE_PS_CANDIDATE_V1 One deterministic 2-fold PS-style joint candidate.
method=string(method);n=size(X,1);rng(data_seed+700000,'twister');perm=randperm(n);n1=floor(n/2);ia=perm(1:n1);ib=perm(n1+1:end);
if method=="REG_TFKM_HD",off=110000;else,off=220000;end
sa=data_seed+off+710000;sb=data_seed+off+720000;t0=tic;
try
 fa=run_hd_multistart_top_v1(X(ia,:),r,g,alpha,method,options,sa);if fa.valid,fa=normalize_fit_contract_v1(fa,numel(ia),g,alpha,method);end
 fb=run_hd_multistart_top_v1(X(ib,:),r,g,alpha,method,options,sb);if fb.valid,fb=normalize_fit_contract_v1(fb,numel(ib),g,alpha,method);end
 if ~(fa.valid&&fb.valid),error('competitor:PSFit','Invalid half fit.');end
 [pb,ok1]=predict_labels_from_hd_fit_v1(X(ia,:),fa,X(ib,:),g);[pa,ok2]=predict_labels_from_hd_fit_v1(X(ib,:),fb,X(ia,:),g);
 [lb,ok3]=predict_labels_from_hd_fit_v1(X(ib,:),fb,X(ib,:),g);[la,ok4]=predict_labels_from_hd_fit_v1(X(ia,:),fa,X(ia,:),g);
 if ~(ok1&&ok2&&ok3&&ok4),error('competitor:PSPredict','Unable to construct complete fold partitions/predictions.');end
 ps_ab=prediction_strength_from_labels_v1(lb,pb,g);ps_ba=prediction_strength_from_labels_v1(la,pa,g);ps=mean([ps_ab ps_ba]);valid=isfinite(ps);msg="";
catch ME
 ps_ab=NaN;ps_ba=NaN;ps=NaN;valid=false;msg=string(ME.identifier)+": "+string(ME.message);
end
out=struct('ps',ps,'ps_ab',ps_ab,'ps_ba',ps_ba,'valid',valid,'message',msg,'runtime_seconds',toc(t0));
end
