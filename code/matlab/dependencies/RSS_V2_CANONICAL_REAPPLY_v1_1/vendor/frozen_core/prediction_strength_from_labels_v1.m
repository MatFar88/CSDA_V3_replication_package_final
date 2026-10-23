function ps=prediction_strength_from_labels_v1(test_labels,pred_labels,g)
%PREDICTION_STRENGTH_FROM_LABELS_V1 Tibshirani-Walther pairwise co-membership score.
test_labels=double(test_labels(:));pred_labels=double(pred_labels(:));vals=zeros(g,1);
for k=1:g
 idx=find(test_labels==k);nk=numel(idx);
 if nk<2,vals(k)=0;continue;end
 pr=pred_labels(idx);u=unique(pr);num=0;
 for j=1:numel(u),m=sum(pr==u(j));num=num+m*(m-1);end
 vals(k)=num/(nk*(nk-1));
end
ps=min(vals);if ~isfinite(ps),ps=NaN;end
end
