function out=merge_hier_stable4_full_pipeline_v1(output_root,scope)
%MERGE_HIER_STABLE4_FULL_PIPELINE_V1 Apply the unchanged frozen pipeline to challenger cubes.
%
% This is DEVELOPMENT evidence on the consumed v11 datasets. It does not
% constitute fresh validation of HIER_STABLE4.

pkg=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(pkg),'-begin');
output_root=char(string(output_root));
if nargin<2 || strlength(string(scope))==0,scope="CLEAN_BRIDGE";end
scope=upper(string(scope));
S=load(fullfile(output_root,'study.mat')); study=S.study;
if scope=="CLEAN_BRIDGE",scenarios=[23001 23004]; tag="clean_bridge";
elseif scope=="ALL_TFKM",scenarios=[23001 23004 23005 23007]; tag="all_tfkm";
else,error('hierStable4:Scope','Unknown scope %s.',scope);end

root=fullfile(output_root,'hier_stable4_full_pipeline_v1'); resdir=fullfile(root,'candidates'); sumdir=fullfile(root,'summary');
if ~isfolder(sumdir),mkdir(sumdir);end
All=table();
for sid=scenarios
    for rep=1:8
        f=fullfile(resdir,sprintf('hier_s%d_rep%02d.csv',sid,rep));
        assert(isfile(f),'Missing challenger job: %s',f);
        T=readtable(f,'TextType','string');
        assert(height(T)==150 && all(logical(T.valid)),'Incomplete/invalid challenger job: %s',f);
        All=[All;T]; %#ok<AGROW>
    end
end
All=sortrows(All,{'scenario_id','replicate','r','g','alpha'});
writetable(All,fullfile(sumdir,"all_hier_candidates_"+tag+".csv"));

current_file=fullfile(output_root,'summary','final_selections.csv');
assert(isfile(current_file),'Missing frozen current final selections: %s',current_file);
Cur=readtable(current_file,'TextType','string');
Cur=Cur(string(Cur.method)=="REG_TFKM_HD",:);

Sel=table();
for sid=scenarios
    for rep=1:8
        U=All(double(All.scenario_id)==sid & double(All.replicate)==rep,:);
        assert(height(U)==150,'Expected 150 challenger rows for sid=%d rep=%d.',sid,rep);
        s=select_final_pipeline_v1(U,"REG_TFKM_HD");
        row=s.selected_row; r0=double(row.r0); g0=double(row.g0); a0=double(row.alpha0);
        tr=U(double(U.r)==r0 & double(U.g)==g0 & abs(double(U.alpha)-a0)<1e-12,:);
        assert(height(tr)==1,'Challenger true-fixed row absent.');
        C=Cur(double(Cur.scenario_id)==sid & double(Cur.replicate)==rep,:);
        assert(height(C)==1,'Frozen current selection not unique for sid=%d rep=%d.',sid,rep);
        contamination=string(row.contamination);
        hier_ari=double(row.ari_inlier); current_ari=double(C.ari_inlier);
        h_over_r=s.r>r0; h_over_g=s.g>g0; h_under_r=s.r<r0; h_under_g=s.g<g0;
        c_over_r=double(C.selected_r)>r0; c_over_g=double(C.selected_g)>g0; c_under_r=double(C.selected_r)<r0; c_under_g=double(C.selected_g)<g0;
        Sel=[Sel;table(sid,rep,contamination,r0,g0,a0, ...
            double(C.selected_r),double(C.selected_g),double(C.selected_alpha),current_ari, ...
            s.r,s.g,s.alpha,hier_ari,hier_ari-current_ari,double(tr.ari_inlier),max(double(tr.ari_inlier)-hier_ari,0), ...
            c_over_r,c_over_g,(c_over_r||c_over_g),c_under_r,c_under_g,(c_under_r||c_under_g), ...
            h_over_r,h_over_g,(h_over_r||h_over_g),h_under_r,h_under_g,(h_under_r||h_under_g), ...
            s.r==r0,s.g==g0,(s.r==r0&&s.g==g0),abs(s.alpha-a0)<=1e-12, ...
            s.stage1_r,s.stage1_g,s.stage1_alpha,logical(s.stage1_guarded),s.stage1_anchor_r, ...
            s.stage2g_r,s.stage2g_g,s.stage2gr_r,s.stage2gr_g,logical(s.g_changed),logical(s.r_changed), ...
            string(row.selected_init_family),double(row.selector_score),double(row.explained),double(row.R2),double(row.W_ratio), ...
            'VariableNames',{'scenario_id','replicate','contamination','r0','g0','alpha0', ...
            'current_r','current_g','current_alpha','current_ari', ...
            'hier_r','hier_g','hier_alpha','hier_ari','delta_ari','hier_true_fixed_ari','hier_positive_regret', ...
            'current_over_r','current_over_g','current_structural_over','current_under_r','current_under_g','current_any_under', ...
            'hier_over_r','hier_over_g','hier_structural_over','hier_under_r','hier_under_g','hier_any_under', ...
            'hier_exact_r','hier_exact_g','hier_exact_rg','hier_exact_alpha', ...
            'hier_stage1_r','hier_stage1_g','hier_stage1_alpha','hier_stage1_guarded','hier_stage1_anchor_r', ...
            'hier_stage2g_r','hier_stage2g_g','hier_stage2gr_r','hier_stage2gr_g','hier_g_changed','hier_r_changed', ...
            'hier_selected_init_family','hier_native_score','hier_explained','hier_R2','hier_W_ratio'})]; %#ok<AGROW>
    end
end
Sel=sortrows(Sel,{'scenario_id','replicate'});
writetable(Sel,fullfile(sumdir,"hier_final_selections_"+tag+".csv"));
Summary=summarize_compare(Sel);
writetable(Summary,fullfile(sumdir,"hier_vs_current_summary_"+tag+".csv"));
Family=summarize_family(All,Sel);
writetable(Family,fullfile(sumdir,"hier_family_usage_"+tag+".csv"));

fprintf('\n============================================================\n');
fprintf('HIER_STABLE4 FULL-PIPELINE DEVELOPMENT MERGE | scope=%s | units=%d\n',scope,height(Sel));
fprintf('Downstream RANK_EN, GMAX_R_E, Stage2, alpha and joint vote are unchanged.\n');
fprintf('============================================================\n');
disp(Summary);
fprintf('\nUNIT-LEVEL CHANGES\n');
disp(Sel(:,{'scenario_id','replicate','contamination','current_r','current_g','current_alpha','current_ari','hier_r','hier_g','hier_alpha','hier_ari','delta_ari','current_structural_over','hier_structural_over','current_any_under','hier_any_under','hier_selected_init_family'}));

out=struct('scope',scope,'candidates',All,'selections',Sel,'summary',Summary,'family_usage',Family, ...
    'development_only',true,'challenger',"HIER_STABLE4");
save(fullfile(sumdir,"hier_full_pipeline_report_"+tag+".mat"),'out','-v7.3');
end

function S=summarize_compare(T)
labels=["OVERALL";unique(string(T.contamination),'stable')]; S=table();
for i=1:numel(labels)
    if labels(i)=="OVERALL",Z=T;else,Z=T(string(T.contamination)==labels(i),:);end
    row=table(labels(i),height(Z), ...
        mean(double(Z.current_ari),'omitnan'),mean(double(Z.hier_ari),'omitnan'),mean(double(Z.delta_ari),'omitnan'), ...
        min(double(Z.current_ari)),min(double(Z.hier_ari)), ...
        mean(double(Z.current_structural_over)),mean(double(Z.hier_structural_over)), ...
        mean(double(Z.current_over_r)),mean(double(Z.hier_over_r)),mean(double(Z.current_over_g)),mean(double(Z.hier_over_g)), ...
        mean(double(Z.current_any_under)),mean(double(Z.hier_any_under)), ...
        mean(double(Z.hier_exact_rg)),mean(double(Z.hier_exact_alpha)), ...
        mean(double(Z.hier_r-r0col(Z)),'omitnan'),mean(double(Z.hier_g-g0col(Z)),'omitnan'), ...
        mean(double(Z.delta_ari)>0.01),mean(double(Z.delta_ari)<-0.01), ...
        'VariableNames',{'scope','N','mean_current_ari','mean_hier_ari','mean_delta_ari','min_current_ari','min_hier_ari', ...
        'current_structural_over_rate','hier_structural_over_rate','current_over_r_rate','hier_over_r_rate','current_over_g_rate','hier_over_g_rate', ...
        'current_any_under_rate','hier_any_under_rate','hier_exact_rg_rate','hier_exact_alpha_rate','mean_hier_r_error','mean_hier_g_error', ...
        'rate_ari_improve_gt_001','rate_ari_degrade_gt_001'});
    S=[S;row]; %#ok<AGROW>
end
end
function x=r0col(T),x=double(T.r0);end
function x=g0col(T),x=double(T.g0);end

function F=summarize_family(All,Sel)
fams=["WINSOR_PCA_TKM";"SIGN_PCA_TKM";"TKM_ORIGINAL";"RP_TKM"];
F=table();
for f=1:numel(fams)
    cand=mean(string(All.selected_init_family)==fams(f));
    fin=mean(string(Sel.hier_selected_init_family)==fams(f));
    F=[F;table(fams(f),cand,fin,'VariableNames',{'family','candidate_winner_rate','final_selected_family_rate'})]; %#ok<AGROW>
end
end
