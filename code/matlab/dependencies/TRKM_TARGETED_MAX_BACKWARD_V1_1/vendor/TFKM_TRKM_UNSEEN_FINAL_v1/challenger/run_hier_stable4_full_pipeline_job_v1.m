function T = run_hier_stable4_full_pipeline_job_v1(output_root,scenario_id,replicate)
%RUN_HIER_STABLE4_FULL_PIPELINE_JOB_V1 Refit one complete 150-candidate TFKM cube.
%
% Uses the ORIGINAL v11 output_root only as a source of frozen study design.
% Writes to a separate challenger directory. The frozen results are never edited.

pkg=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(pkg),'-begin');
output_root=char(string(output_root));
S=load(fullfile(output_root,'study.mat'));
study=S.study; options=S.options; catalog=S.catalog; pairs=S.pairs; alpha_grid=S.alpha_grid;
scenario_id=double(scenario_id); replicate=double(replicate);
assert(ismember(replicate,double(study.validation_reps)),'Replicate %d is not in frozen validation reps.',replicate);
ci=find(double(catalog.scenario_id)==scenario_id,1);
assert(~isempty(ci),'Scenario %d not found in frozen catalog.',scenario_id);
cfg=table2struct(catalog(ci,:));

root=fullfile(output_root,'hier_stable4_full_pipeline_v1');
resdir=fullfile(root,'candidates'); if ~isfolder(resdir),mkdir(resdir);end
sumdir=fullfile(root,'summary'); if ~isfolder(sumdir),mkdir(sumdir);end
outfile=fullfile(resdir,sprintf('hier_s%d_rep%02d.csv',scenario_id,replicate));
checkpoint=fullfile(resdir,sprintf('hier_s%d_rep%02d_checkpoint.mat',scenario_id,replicate));

T=table(); if isfile(outfile),T=readtable(outfile,'TextType','string');end
if isfile(checkpoint)
    Q=load(checkpoint,'T_checkpoint');
    if isfield(Q,'T_checkpoint'),T=merge_unique(T,Q.T_checkpoint);end
end
keys=make_keys(T);
expected=height(pairs)*height(alpha_grid);
seed=double(study.validation_base_seed)+scenario_id*100000+replicate;
pair_seed=seed+110000;
sim=generate_fixed_dataset_v1(cfg,seed);
primary=cfg.primary_tfkm; sentinel=cfg.sentinel_tfkm;
t0=tic;

fprintf('\n============================================================\n');
fprintf('HIER_STABLE4 FULL-PIPELINE JOB | sid=%d rep=%d | %s | candidates=%d | already=%d\n', ...
    scenario_id,replicate,char(string(cfg.contamination)),expected,height(T));
fprintf('Only TFKM multistart arbitration changes; ALS, seeds and downstream selectors are frozen.\n');
fprintf('============================================================\n');

for pp=1:height(pairs)
    r=double(pairs.r(pp)); g=double(pairs.g(pp));
    for aa=1:height(alpha_grid)
        alpha=double(alpha_grid.alpha(aa));
        key=sprintf('%d|%d|%d|%.6f',replicate,r,g,alpha);
        if ismember(string(key),keys),continue;end
        ft=tic;
        fprintf('  [%3d/%3d] HIER TFKM (%d,%d,%.3f) ... ',height(T)+1,expected,r,g,alpha);
        try
            fit=run_hd_multistart_top_hier_stable4_v1(sim.X,r,g,alpha,options,pair_seed);
            if fit.valid
                try
                    fit=normalize_fit_contract_v1(fit,size(sim.X,1),g,alpha,"REG_TFKM_HD");
                catch ME
                    fit.valid=false; fit.message=string(ME.identifier)+": "+string(ME.message);
                end
            end
        catch ME
            fit=struct('valid',false,'message',string(ME.identifier)+": "+string(ME.message), ...
                'objective',NaN,'selector',"HIER_STABLE4_UNCAUGHT_ERROR",'N_starts',100,'valid_starts',0, ...
                'selected_init_family',"",'selected_objective_rank',NaN,'selector_score',NaN, ...
                'selector_cv_score',NaN,'tauW',options.tauW,'tauT',options.tauT);
        end
        elapsed=toc(ft);
        row=compute_joint_candidate_row_v1(sim,fit,cfg,"REG_TFKM_HD",replicate,seed,r,g,alpha,elapsed,primary,sentinel);
        T=[T;row]; %#ok<AGROW>
        keys(end+1,1)=string(key); %#ok<AGROW>
        T_checkpoint=T; save(checkpoint,'T_checkpoint','-v7.3'); %#ok<NASGU>
        if row.valid
            fprintf('ARI=%.3f | family=%s | N=%.4f | %.1fs\n',row.ari_inlier,char(row.selected_init_family),row.selector_score,elapsed);
        else
            fprintf('INVALID | %.1fs | %s\n',elapsed,char(row.message));
        end
    end
end

T=sortrows(T,{'r','g','alpha'}); writetable(T,outfile);
assert(height(T)==expected,'Job ended with %d/%d rows.',height(T),expected);
assert(all(logical(T.valid)),'At least one challenger candidate is invalid. Inspect %s.',outfile);
assert(all(double(T.N_starts)==100),'Expected 100 starts per candidate.');
assert(all(string(T.internal_selector)=="HIER_STABLE4"),'Selector contract mismatch.');
if isfile(checkpoint),delete(checkpoint);end
fprintf('HIER JOB COMPLETE | sid=%d rep=%d | %d candidates | %.1f min\n',scenario_id,replicate,height(T),toc(t0)/60);
end

function keys=make_keys(T)
if isempty(T),keys=strings(0,1);return;end
keys=string(double(T.replicate))+"|"+string(double(T.r))+"|"+string(double(T.g))+"|"+compose('%.6f',double(T.alpha));
end
function T=merge_unique(A,B)
if isempty(A),T=B;return;elseif isempty(B),T=A;return;end
T=[A;B]; k=make_keys(T); [~,ia]=unique(k,'stable'); T=T(sort(ia),:);
end
