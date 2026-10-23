function OUT = MERGE_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1(documents_root)
% Merge all completed worker files from the weak true-state pool audit.
if nargin<1,documents_root="";end
[root,outroot,~]=load_study_local(documents_root); %#ok<ASGLU>
dir0=fullfile(outroot,'fresh_v2_diagnostic_suite_v1','weak_true_state_pool_audit_v1');assert(isfolder(dir0),'Weak pool audit output directory missing.');
D=dir(fullfile(dir0,'worker_*_of_*_pool.csv'));assert(~isempty(D),'No worker pool files found.');Pool=table();
for k=1:numel(D),Q=readtable(fullfile(D(k).folder,D(k).name),'TextType','string');Pool=[Pool;Q];end %#ok<AGROW>
Pool=unique(Pool,'rows','stable');Pool=sortrows(Pool,{'method','scenario_id','replicate'});assert(height(Pool)==300,'Expected 300 unique method-unit pool rows, found %d.',height(Pool));
F=dir(fullfile(dir0,'worker_*_of_*_family.csv'));Fam=table();for k=1:numel(F),Q=readtable(fullfile(F(k).folder,F(k).name),'TextType','string');Fam=[Fam;Q];end %#ok<AGROW>
Fam=unique(Fam,'rows','stable');Fam=sortrows(Fam,{'method','scenario_id','replicate','family_index'});assert(height(Fam)==1200,'Expected 1200 family rows, found %d.',height(Fam));
writetable(Pool,fullfile(dir0,'WEAK_TRUESTATE_POOL_ALL.csv'));writetable(Fam,fullfile(dir0,'WEAK_TRUESTATE_FAMILY_ALL.csv'));
[G,m,pn,rho,mech]=findgroups(Pool.method,Pool.p_over_n,Pool.rho,Pool.mechanism_POSTHOC);N=splitapply(@numel,Pool.scenario_id,G);S=table(m,pn,rho,mech,N,'VariableNames',{'method','p_over_n','rho','mechanism_POSTHOC','N'});
[G2,m2,pn2,rho2]=findgroups(Pool.method,Pool.p_over_n,Pool.rho);Den=splitapply(@numel,Pool.scenario_id,G2);K=table(m2,pn2,rho2,Den,'VariableNames',{'method','p_over_n','rho','Den'});S=outerjoin(S,K,'Keys',{'method','p_over_n','rho'},'MergeKeys',true,'Type','left');S.rate_within_group=S.N./S.Den;S=sortrows(S,{'method','p_over_n','rho','mechanism_POSTHOC'});
[G3,m3]=findgroups(Pool.method);N3=splitapply(@numel,Pool.scenario_id,G3);sel=splitapply(@(x)mean(double(x),'omitnan'),Pool.selected_truth_ARI_POSTHOC,G3);best=splitapply(@(x)mean(double(x),'omitnan'),Pool.best_truth_ARI_POSTHOC,G3);r90=splitapply(@(x)mean(double(x)>=.90),Pool.selected_truth_ARI_POSTHOC,G3);b90=splitapply(@(x)mean(double(x)>=.90),Pool.best_truth_ARI_POSTHOC,G3);arb=splitapply(@(x)mean(string(x)=="ARBITRATION_FAILURE"),Pool.mechanism_POSTHOC,G3);bas=splitapply(@(x)mean(string(x)=="BASIN_GENERATION_FAILURE"),Pool.mechanism_POSTHOC,G3);good=splitapply(@(x)mean(string(x)=="SELECTED_GOOD"),Pool.mechanism_POSTHOC,G3);O=table(m3,N3,sel,best,r90,b90,good,arb,bas,'VariableNames',{'method','N','mean_selected_truth_ARI','mean_best_start_truth_ARI','selected_ge090_rate','best_start_ge090_rate','selected_good_rate','arbitration_failure_rate','basin_generation_failure_rate'});
writetable(S,fullfile(dir0,'WEAK_TRUESTATE_MECHANISM_BY_GEOMETRY.csv'));writetable(O,fullfile(dir0,'WEAK_TRUESTATE_POOL_OVERALL.csv'));
fid=fopen(fullfile(dir0,'WEAK_TRUESTATE_POOL_REPORT.txt'),'w');assert(fid>0);c=onCleanup(@()fclose(fid));fprintf(fid,'WEAK TRUE-STATE 100-START POOL AUDIT\n');fprintf(fid,'Truth open; development forensic refit only. V1 unchanged.\n\nOVERALL\n%s\nBY GEOMETRY\n%s\n',evalc('disp(O)'),evalc('disp(S)'));
fprintf('\n========== WEAK TRUE-STATE POOL OVERALL ==========\n');disp(O);fprintf('\n========== MECHANISM BY GEOMETRY =================\n');disp(S);fprintf('MERGE COMPLETE: %s\n',dir0);
OUT=struct('pool',Pool,'family',Fam,'overall',O,'mechanism_by_geometry',S,'output_root',string(dir0));
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end;assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
