function [score,ok,msg]=bqs_smooth_score_v1(Zfit,labels,trim,Zscore,g)
%BQS_SMOOTH_SCORE_V1 Smooth quadratic score in common view.
% Cluster parameters use retained points; all Zscore rows are evaluated.
score=NaN;ok=false;msg="";labels=double(labels(:));trim=logical(trim(:));d=size(Zfit,2);keep=~trim;h=sum(keep);
if h<=g*(d+1),msg="too_few_retained";return;end
K=g;mu=zeros(K,d);Rcell=cell(K,1);logdet=zeros(K,1);pii=zeros(K,1);
for k=1:K
 q=keep & labels==k;nk=sum(q);if nk<=d,msg="cluster_too_small";return;end
 Zk=Zfit(q,:);mu(k,:)=mean(Zk,1);S=cov(Zk,1);S=(S+S')/2;
 [R,pd]=chol(S);if pd~=0||any(~isfinite(R(:))),msg="non_pd_cov";return;end
 Rcell{k}=R;logdet(k)=2*sum(log(diag(R)));pii(k)=nk/h;
end
if any(pii<=0),msg="zero_proportion";return;end
N=size(Zscore,1);QS=zeros(N,K);
for k=1:K
 D=Zscore-mu(k,:);A=D/Rcell{k};mah=sum(A.^2,2);
 QS(:,k)=log(pii(k))-.5*logdet(k)-.5*mah;
end
mx=max(QS,[],2);E=exp(QS-mx);den=sum(E,2);if any(~isfinite(den)|den<=0),msg="softmax_failure";return;end
T=E./den;score=mean(sum(T.*QS,2));ok=isfinite(score);if ~ok,msg="nonfinite_score";end
end
