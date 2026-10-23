function S=bqs_boot_summary_v1(B,required_total,boot_na_share)
%BQS_BOOT_SUMMARY_V1 Summarize candidate bootstrap scores and missing-score filter.
if isempty(B),S=table();return;end
[G,r,g,a]=findgroups(double(B.r),double(B.g),double(B.alpha));n_total=splitapply(@numel,B.score,G);n_obs=splitapply(@(x)sum(isfinite(x)),double(B.score),G);n_missing=n_total-n_obs;
mean_score=splitapply(@localmean,double(B.score),G);lower_qnt=splitapply(@localq,double(B.score),G);
eligible=n_total>=required_total & (n_missing./max(n_total,1))<=boot_na_share & isfinite(lower_qnt);
S=table(r,g,a,n_total,n_obs,n_missing,mean_score,lower_qnt,eligible,'VariableNames',{'r','g','alpha','n_total','n_obs','n_missing','mean_score','lower_qnt','eligible'});
end
function m=localmean(x)
x=x(isfinite(x));if isempty(x),m=NaN;else,m=mean(x);end
end
function q=localq(x)
x=x(isfinite(x));if isempty(x),q=NaN;else,q=prctile(x,2.5);end
end
