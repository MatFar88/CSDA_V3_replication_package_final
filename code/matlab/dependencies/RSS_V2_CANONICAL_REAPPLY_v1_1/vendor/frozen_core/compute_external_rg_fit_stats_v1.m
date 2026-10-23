function stats=compute_external_rg_fit_stats_v1(X,fit,method)
%COMPUTE_EXTERNAL_RG_FIT_STATS_V1 Legacy-compatible truth-blind RG statistics.
% Uses only X and the fitted labels/trim/subspace; no truth information.
method=string(method);
assert(fit.valid,'Cannot score an invalid fit.');
in=~logical(fit.trim(:));
labels=double(fit.labels(in));
Xk=X(in,:); mu=mean(Xk,1); Xc=Xk-mu;
U=fit.U; Y=Xc*U;
g=max(labels); r=size(U,2);
cent=zeros(g,r);
for k=1:g
    q=labels==k;
    if any(q),cent(k,:)=mean(Y(q,:),1);end
end
TSS=sum(Xc(:).^2);
Yc=Y-mean(Y,1); subTSS=sum(Yc(:).^2);
if TSS>0, W_ratio=subTSS/TSS; else, W_ratio=NaN; end

if method=="REG_TFKM_HD"
    R=Y-cent(labels,:);
    within=sum(R(:).^2);
    if subTSS>0,W_raw=within/subTSS;else,W_raw=NaN;end
    R2=1-W_raw;
elseif method=="REG_TRKM_HD"
    Xhat=cent(labels,:)*U';
    R=Xc-Xhat;
    RSS=sum(R(:).^2);
    if TSS>0,W_raw=RSS/TSS;else,W_raw=NaN;end
    R2=1-W_raw;
else
    error('selectRG:Method','Unknown method %s.',method);
end
stats=struct('W_raw',W_raw,'W_ratio',W_ratio,'R2',R2,'n_keep',sum(in), ...
    'subspace_variance',subTSS,'total_variance',TSS);
end
