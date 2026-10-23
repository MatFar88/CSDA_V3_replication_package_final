function row=compute_sealed_truth_row_unseen_v1(sim,fit,cfg,method,rep,r,g,alpha)
%COMPUTE_SEALED_TRUTH_ROW_UNSEEN_V1 Evaluation-only row. NEVER read by selector.
method=string(method);
truth=sim.labels(:); out=logical(sim.outlier(:)); intrue=~out;
if isfield(fit,'valid') && fit.valid
    pred=fit.labels(:); trim=logical(fit.trim(:));
    ari=adjusted_rand_index_v1(truth(intrue),pred(intrue));
    tp=sum(trim&out); fp=sum(trim&~out); fn=sum(~trim&out);
    if any(out)
        pr=tp/max(tp+fp,1); re=tp/max(tp+fn,1); f1=2*pr*re/max(pr+re,eps);
    else
        pr=NaN; re=NaN; f1=NaN;
    end
else
    ari=NaN; pr=NaN; re=NaN; f1=NaN;
end
row=table(double(cfg.scenario_id),double(rep),method,r,g,alpha, ...
 double(cfg.r0),double(cfg.c0),double(cfg.alpha0),string(cfg.contamination), ...
 ari,pr,re,f1, ...
 'VariableNames',{'scenario_id','replicate','method','r','g','alpha', ...
 'r0','g0','alpha0','contamination','ari_inlier','trim_precision','trim_recall','trim_f1'});
end
