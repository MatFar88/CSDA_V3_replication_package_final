function run_final_fresh_session_v1(session,output_root)
%RUN_FINAL_FRESH_SESSION_V1 Run/resume one untouched validation replicate.
pkg=fileparts(mfilename('fullpath')); addpath(genpath(pkg),'-begin'); output_root=char(string(output_root));
S=load(fullfile(output_root,'study.mat')); study=S.study; options=S.options; catalog=S.catalog; methods=S.methods; pairs=S.pairs; alpha_grid=S.alpha_grid; manifest=S.manifest;
certfile=fullfile(output_root,'dgp_certification','DGP_CERTIFICATION_PASS.mat');
assert(isfile(certfile),'Run certify_fresh_dgp_catalog_v1 first.'); Cpass=load(certfile,'all_pass'); assert(Cpass.all_pass,'DGP certification did not pass: validation is locked.');
session=upper(string(session)); m=find(string(manifest.session)==session,1); assert(~isempty(m),'Unknown session %s.',session);
rep=double(manifest.replicate(m)); expected=double(manifest.expected_candidate_fits(m));
outfile=fullfile(output_root,'results',sprintf('joint_candidates_%s.csv',char(session)));
checkpoint=fullfile(output_root,'results',sprintf('joint_candidates_%s_checkpoint.mat',char(session)));
R=table(); if isfile(outfile),R=readtable(outfile,'TextType','string');end
if isfile(checkpoint),Q=load(checkpoint,'R_checkpoint'); if isfield(Q,'R_checkpoint'),R=merge_unique(R,Q.R_checkpoint);end,end
keys=make_keys(R); done=numel(keys); t0=tic;
fprintf('\n============================================================\nFINAL FRESH session %s | rep=%d | candidates=%d | already=%d\n============================================================\n',session,rep,expected,done);
for ii=1:height(catalog)
    cfg=table2struct(catalog(ii,:)); seed=double(study.validation_base_seed)+double(cfg.scenario_id)*100000+rep; sim=[]; changed=false;
    fprintf('\n--- sid=%d source=%d | %s | p/n=%.2f | true=(%d,%d,%.3f) | delta=%.4g rho=%.4g\n', ...
        cfg.scenario_id,cfg.source_confirm_id,string(cfg.contamination),cfg.p_over_n,cfg.r0,cfg.c0,cfg.alpha0,cfg.delta,cfg.rho);
    for mm=1:height(methods)
        method=string(methods.method(mm));
        if method=="REG_TFKM_HD", pair_seed=seed+110000; primary=cfg.primary_tfkm; sentinel=cfg.sentinel_tfkm; else, pair_seed=seed+220000; primary=cfg.primary_trkm; sentinel=cfg.sentinel_trkm; end
        for pp=1:height(pairs)
            r=double(pairs.r(pp)); g=double(pairs.g(pp));
            for aa=1:height(alpha_grid)
                alpha=double(alpha_grid.alpha(aa)); key=sprintf('%d|%d|%s|%d|%d|%.6f',cfg.scenario_id,rep,char(method),r,g,alpha);
                if ismember(string(key),keys),continue;end
                if isempty(sim),sim=generate_fixed_dataset_v1(cfg,seed);end
                ft=tic; fprintf('  [%4d/%4d] %-13s (%d,%d,%.3f) ... ',done+1,expected,char(method),r,g,alpha);
                try
                    fit=run_hd_multistart_top_v1(sim.X,r,g,alpha,method,options,pair_seed);
                    if fit.valid
                        try, fit=normalize_fit_contract_v1(fit,size(sim.X,1),g,alpha,method); catch ME, fit.valid=false; fit.message=string(ME.identifier)+": "+string(ME.message); end
                    end
                catch ME
                    fit=struct('valid',false,'message',string(ME.identifier)+": "+string(ME.message),'objective',NaN,'selector',"UNCAUGHT_ERROR", ...
                        'N_starts',100,'valid_starts',0,'selected_init_family',"",'selected_objective_rank',NaN,'selector_score',NaN,'selector_cv_score',NaN,'tauW',options.tauW,'tauT',options.tauT);
                end
                elapsed=toc(ft); row=compute_joint_candidate_row_v1(sim,fit,cfg,method,rep,seed,r,g,alpha,elapsed,primary,sentinel);
                R=[R;row]; keys(end+1,1)=string(key); done=done+1; changed=true; R_checkpoint=R; save(checkpoint,'R_checkpoint','-v7.3'); %#ok<AGROW,NASGU>
                if row.valid,fprintf('ARI=%.3f | %.1fs\n',row.ari_inlier,elapsed);else,fprintf('INVALID | %.1fs | %s\n',elapsed,char(row.message));end
            end
        end
    end
    if changed,writetable(R,outfile);end
    fprintf('Scenario checkpoint | completed=%d/%d | elapsed %.1f min\n',done,expected,toc(t0)/60);
end
R=sortrows(R,{'scenario_id','replicate','method','r','g','alpha'}); writetable(R,outfile);
assert(height(R)==expected,'Session %s ended with %d/%d rows.',session,height(R),expected); if isfile(checkpoint),delete(checkpoint);end
fprintf('SESSION %s COMPLETE | %d candidates | %.1f min\n',session,height(R),toc(t0)/60);
end
function keys=make_keys(T)
if isempty(T),keys=strings(0,1);return;end
keys=string(double(T.scenario_id))+"|"+string(double(T.replicate))+"|"+string(T.method)+"|"+string(double(T.r))+"|"+string(double(T.g))+"|"+compose('%.6f',double(T.alpha));
end
function T=merge_unique(A,B)
if isempty(A),T=B;return;elseif isempty(B),T=A;return;end
T=[A;B]; k=make_keys(T); [~,ia]=unique(k,'stable'); T=T(sort(ia),:);
end
