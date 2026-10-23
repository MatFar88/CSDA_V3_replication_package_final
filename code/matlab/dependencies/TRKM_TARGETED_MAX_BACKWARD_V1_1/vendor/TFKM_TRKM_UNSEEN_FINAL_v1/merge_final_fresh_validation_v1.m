function out=merge_final_fresh_validation_v1(output_root)
%MERGE_FINAL_FRESH_VALIDATION_V1 Merge fresh cube and apply frozen final pipeline.
pkg=fileparts(mfilename('fullpath')); addpath(genpath(pkg),'-begin'); output_root=char(string(output_root));
S=load(fullfile(output_root,'study.mat')); study=S.study; catalog=S.catalog; methods=S.methods; manifest=S.manifest;
All=table();
for i=1:height(manifest)
    f=fullfile(output_root,'results',sprintf('joint_candidates_%s.csv',char(manifest.session(i))));
    assert(isfile(f),'Missing session file: %s',f); T=readtable(f,'TextType','string');
    assert(height(T)==double(manifest.expected_candidate_fits(i)),'Incomplete %s: %d rows.',f,height(T)); All=[All;T]; %#ok<AGROW>
end
All=sortrows(All,{'scenario_id','replicate','method','r','g','alpha'});
expected=sum(double(manifest.expected_candidate_fits)); assert(height(All)==expected,'Merged cube %d/%d.',height(All),expected);
writetable(All,fullfile(output_root,'summary','all_joint_candidates.csv'));
Sel=table();
for ii=1:height(catalog)
    sid=double(catalog.scenario_id(ii));
    for rep=reshape(double(study.validation_reps),1,[])
        for mm=1:height(methods)
            method=string(methods.method(mm)); U=All(double(All.scenario_id)==sid & double(All.replicate)==rep & string(All.method)==method,:);
            assert(height(U)==150,'Expected 150 cube rows for sid=%d rep=%d %s.',sid,rep,method);
            s=select_final_pipeline_v1(U,method); row=s.selected_row; r0=double(row.r0); g0=double(row.g0); a0=double(row.alpha0);
            tr=U(double(U.r)==r0 & double(U.g)==g0 & abs(double(U.alpha)-a0)<1e-12,:); assert(height(tr)==1,'True-fixed row absent.');
            ari=double(row.ari_inlier); ari_true=double(tr.ari_inlier); preg=max(0,ari_true-ari); structural_over=(s.r>r0)||(s.g>g0); low_ari=ari<double(study.forensic_ari_threshold);
            Sel=[Sel;table(sid,rep,method,r0,g0,a0,s.r,s.g,s.alpha,ari,ari_true,preg,double(row.trim_precision),double(row.trim_recall),double(row.trim_f1), ...
                s.stage1_r,s.stage1_g,s.stage1_alpha,logical(s.stage1_guarded),s.stage1_anchor_r,s.stage2g_r,s.stage2g_g,s.stage2gr_r,s.stage2gr_g,logical(s.g_changed),logical(s.r_changed), ...
                logical(s.native010_one_step_triggered),s.native010_r,s.native010_g,s.native010_delta_R2,s.native010_delta_cv, ...
                s.r>=r0,s.g>=g0,(s.r>=r0&&s.g>=g0),(s.r==r0&&s.g==g0),abs(s.alpha-a0)<=1e-12,(s.r==r0&&s.g==g0&&abs(s.alpha-a0)<=1e-12), ...
                max(s.r-r0,0),max(s.g-g0,0),logical(structural_over),logical(low_ari),logical(structural_over&&ari>=study.forensic_ari_threshold), ...
                'VariableNames',{'scenario_id','replicate','method','r0','g0','alpha0','selected_r','selected_g','selected_alpha','ari_inlier','ari_true_fixed','positive_regret','trim_precision','trim_recall','trim_f1', ...
                'stage1_r','stage1_g','stage1_alpha','stage1_guarded','stage1_anchor_r','stage2g_r','stage2g_g','stage2gr_r','stage2gr_g','g_changed','r_changed', ...
                'native010_one_step_triggered','native010_r','native010_g','native010_delta_R2','native010_delta_cv', ...
                'no_under_r','no_under_g','no_under_rg','exact_rg','exact_alpha','exact_joint','r_overshoot','g_overshoot','structural_over','low_ari_lt_095','over_ari_ge_095'})]; %#ok<AGROW>
        end
    end
end
writetable(Sel,fullfile(output_root,'summary','final_selections.csv'));
Overall=summarize(Sel,{'method'}); ByScenario=summarize(Sel,{'scenario_id','method'});
writetable(Overall,fullfile(output_root,'summary','final_summary_by_method.csv')); writetable(ByScenario,fullfile(output_root,'summary','final_summary_by_scenario.csv'));
fprintf('\n============================================================\nFINAL FRESH VALIDATION MERGED | %d candidate rows | %d selected units\n============================================================\n',height(All),height(Sel));
fprintf('\nOVERALL FINAL SUMMARY\n'); disp(Overall);
fprintf('\nBY-SCENARIO FINAL SUMMARY\n'); disp(ByScenario);
fprintf('\nOVERSELECTION QUALITY (primary interpretive endpoint: ARI >= %.2f when structurally overselected)\n',study.forensic_ari_threshold);
Q=Overall(:,{'method','N','structural_over_rate','mean_ari_over','min_ari_over','over_ari_ge_095_rate','low_ari_rate','mean_positive_regret','q90_positive_regret'}); disp(Q);
out=struct('candidates',All,'selections',Sel,'overall',Overall,'by_scenario',ByScenario); save(fullfile(output_root,'summary','final_fresh_report.mat'),'out','-v7.3');
fprintf('NEXT: run_final_fresh_forensics_session_v1(A..H,output_root), then merge_final_fresh_forensics_v1(output_root).\n');
end

function S=summarize(T,groups)
[G,keys]=findgroups(T(:,groups)); S=table();
for k=1:max(G)
    q=G==k; Z=T(q,:); row=keys(k,:); row.N=height(Z); row.mean_ari=mean(double(Z.ari_inlier),'omitnan'); row.min_ari=min(double(Z.ari_inlier));
    row.mean_true_fixed_ari=mean(double(Z.ari_true_fixed),'omitnan'); row.mean_positive_regret=mean(double(Z.positive_regret),'omitnan'); row.q90_positive_regret=quantile_local(double(Z.positive_regret),.90); row.mean_trim_precision=mean(double(Z.trim_precision),'omitnan'); row.mean_trim_recall=mean(double(Z.trim_recall),'omitnan'); row.mean_trim_f1=mean(double(Z.trim_f1),'omitnan');
    row.no_under_rg_rate=mean(double(Z.no_under_rg)); row.exact_rg_rate=mean(double(Z.exact_rg)); row.exact_alpha_rate=mean(double(Z.exact_alpha)); row.exact_joint_rate=mean(double(Z.exact_joint));
    row.mean_r_overshoot=mean(double(Z.r_overshoot)); row.mean_g_overshoot=mean(double(Z.g_overshoot)); row.structural_over_rate=mean(double(Z.structural_over)); row.low_ari_rate=mean(double(Z.low_ari_lt_095));
    row.native010_one_step_trigger_rate=mean(double(Z.native010_one_step_triggered));
    ov=logical(Z.structural_over); if any(ov),row.mean_ari_over=mean(double(Z.ari_inlier(ov)),'omitnan');row.min_ari_over=min(double(Z.ari_inlier(ov)));row.over_ari_ge_095_rate=mean(double(Z.ari_inlier(ov))>=.95);else,row.mean_ari_over=NaN;row.min_ari_over=NaN;row.over_ari_ge_095_rate=NaN;end
    S=[S;row]; %#ok<AGROW>
end
end
function q=quantile_local(x,p),x=sort(x(isfinite(x)));if isempty(x),q=NaN;return;end;pos=1+(numel(x)-1)*p;lo=floor(pos);hi=ceil(pos);if lo==hi,q=x(lo);else,q=x(lo)+(pos-lo)*(x(hi)-x(lo));end,end
