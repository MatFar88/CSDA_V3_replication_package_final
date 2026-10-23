function [lab,ok]=predict_labels_from_hd_fit_v1(Xtrain,fit,Xnew,g)
%PREDICT_LABELS_FROM_HD_FIT_V1 Complete nearest-centroid assignment for PS.
ok=false;lab=nan(size(Xnew,1),1);
if ~isfield(fit,'valid')||~fit.valid||~isfield(fit,'U')||isempty(fit.U),return;end
keep=~logical(fit.trim(:));labels=double(fit.labels(:));
if sum(keep)<=g,return;end
mu=mean(Xtrain(keep,:),1);Ytr=(Xtrain-mu)*fit.U;Ynew=(Xnew-mu)*fit.U;
cent=zeros(g,size(Ytr,2));
for k=1:g
 q=keep & labels==k;if ~any(q),return;end
 cent(k,:)=mean(Ytr(q,:),1);
end
D=zeros(size(Ynew,1),g);for k=1:g,Q=Ynew-cent(k,:);D(:,k)=sum(Q.^2,2);end
[~,lab]=min(D,[],2);ok=all(isfinite(lab));
end
