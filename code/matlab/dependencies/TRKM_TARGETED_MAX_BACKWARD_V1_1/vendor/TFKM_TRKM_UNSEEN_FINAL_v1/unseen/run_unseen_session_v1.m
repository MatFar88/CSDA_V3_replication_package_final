function run_unseen_session_v1(session,output_root)
%RUN_UNSEEN_SESSION_V1 Run/resume one unseen replicate for BOTH TFKM and TRKM.
% Operational candidate table contains no truth columns.
% Truth metrics are written separately into a SEALED folder and are never
% read by operational selection.

pkg=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(pkg),'-begin');
output_root=char(string(output_root));
S=load(fullfile(output_root,'study.mat'));
study=S.study; options=S.options; catalog=S.catalog; methods=S.methods;
pairs=S.pairs; alpha_grid=S.alpha_grid; manifest=S.manifest;

assert(~study.truth_open_allowed, ...
    'Truth-open flag unexpectedly enabled before unseen operational selection.');

session=upper(string(session));
m=find(string(manifest.session)==session,1);
assert(~isempty(m),'Unknown session %s.',session);
rep=double(manifest.replicate(m));
expected=double(manifest.expected_total(m));

opfile=fullfile(output_root,'candidates',sprintf('unseen_candidates_%s.csv',char(session)));
truthfile=fullfile(output_root,'_SEALED_TRUTH_DO_NOT_OPEN_BEFORE_SELECTION', ...
    sprintf('sealed_truth_%s.csv',char(session)));
checkpoint=fullfile(output_root,'candidates',sprintf('unseen_%s_checkpoint.mat',char(session)));

R=table(); Z=table();
if isfile(opfile), R=readtable(opfile,'TextType','string'); end
if isfile(truthfile), Z=readtable(truthfile,'TextType','string'); end
if isfile(checkpoint)
    Q=load(checkpoint,'R_checkpoint','Z_checkpoint');
    if isfield(Q,'R_checkpoint'), R=merge_unique(R,Q.R_checkpoint); end
    if isfield(Q,'Z_checkpoint'), Z=merge_unique(Z,Q.Z_checkpoint); end
end

opkeys=make_keys(R); truthkeys=make_keys(Z);
donekeys=intersect(opkeys,truthkeys,'stable');
t0=tic;

fprintf('\n============================================================\n');
fprintf('UNSEEN FINAL session %s | replicate=%d | total candidates=%d | complete=%d\n', ...
    session,rep,expected,numel(donekeys));
fprintf('TFKM = HIER_STABLE4 (+ J saved) | TRKM = unchanged RV_CV\n');
fprintf('Truth is SEALED and never printed or used for selection.\n');
fprintf('============================================================\n');

for ii=1:height(catalog)
    cfg=table2struct(catalog(ii,:));
    sid=double(cfg.scenario_id);
    seed=double(study.validation_base_seed)+sid*100000+rep;
    sim=[];
    fprintf('\n--- unseen scenario %d | rep=%d | n=%d p=%d ---\n',sid,rep,cfg.n,cfg.p);

    for mm=1:height(methods)
        method=string(methods.method(mm));
        if method=="REG_TFKM_HD"
            pair_seed=seed+110000;
        else
            pair_seed=seed+220000;
        end

        for pp=1:height(pairs)
            r=double(pairs.r(pp)); g=double(pairs.g(pp));
            for aa=1:height(alpha_grid)
                alpha=double(alpha_grid.alpha(aa));
                key=string(sprintf('%d|%d|%s|%d|%d|%.6f',sid,rep,char(method),r,g,alpha));
                if ismember(key,donekeys), continue; end
                if isempty(sim), sim=generate_fixed_dataset_v1(cfg,seed); end

                ft=tic;
                fprintf('  %-13s (%d,%d,%.3f) ... ',char(method),r,g,alpha);
                try
                    if method=="REG_TFKM_HD"
                        fit=run_hd_multistart_top_hier_stable4_unseen_v1(sim.X,r,g,alpha,options,pair_seed);
                    else
                        fit=run_hd_multistart_top_v1(sim.X,r,g,alpha,method,options,pair_seed);
                    end
                    if fit.valid
                        try
                            fit=normalize_fit_contract_v1(fit,size(sim.X,1),g,alpha,method);
                        catch ME
                            fit.valid=false;
                            fit.message=string(ME.identifier)+": "+string(ME.message);
                        end
                    end
                catch ME
                    fit=struct('valid',false,'message',string(ME.identifier)+": "+string(ME.message), ...
                        'objective',NaN,'selector',"UNCAUGHT_ERROR",'N_starts',100,'valid_starts',0, ...
                        'selected_init_family',"",'selected_objective_rank',NaN,'selector_score',NaN, ...
                        'selector_cv_score',NaN,'tauW',options.tauW,'tauT',options.tauT);
                end
                elapsed=toc(ft);

                op=compute_operational_candidate_row_unseen_v1( ...
                    sim.X,fit,sid,rep,seed,method,r,g,alpha,elapsed);
                tr=compute_sealed_truth_row_unseen_v1(sim,fit,cfg,method,rep,r,g,alpha);

                R=[R;op]; %#ok<AGROW>
                Z=[Z;tr]; %#ok<AGROW>
                R=merge_unique(table(),R);
                Z=merge_unique(table(),Z);
                opkeys=make_keys(R); truthkeys=make_keys(Z);
                donekeys=intersect(opkeys,truthkeys,'stable');

                R_checkpoint=R; Z_checkpoint=Z; %#ok<NASGU>
                save(checkpoint,'R_checkpoint','Z_checkpoint','-v7.3');

                if op.valid
                    if method=="REG_TFKM_HD"
                        fprintf('H=%.4f J=%.4f | %s | %.1fs\n', ...
                            op.selector_score,op.hier_cross_j_score,char(op.selected_init_family),elapsed);
                    else
                        fprintf('RV=%.4f | %s | %.1fs\n', ...
                            op.native_score,char(op.selected_init_family),elapsed);
                    end
                else
                    fprintf('INVALID | %.1fs | %s\n',elapsed,char(op.message));
                end
            end
        end
    end

    writetable(sortrows(R,{'scenario_id','method','r','g','alpha'}),opfile);
    writetable(sortrows(Z,{'scenario_id','method','r','g','alpha'}),truthfile);
    fprintf('Scenario checkpoint | complete=%d/%d | elapsed %.1f min\n', ...
        numel(donekeys),expected,toc(t0)/60);
end

R=sortrows(R,{'scenario_id','replicate','method','r','g','alpha'});
Z=sortrows(Z,{'scenario_id','replicate','method','r','g','alpha'});
writetable(R,opfile); writetable(Z,truthfile);
assert(height(R)==expected && height(Z)==expected, ...
    'Session %s ended with operational=%d truth=%d expected=%d.', ...
    session,height(R),height(Z),expected);
assert(all(logical(R.valid)),'At least one unseen candidate is invalid. Inspect %s.',opfile);
assert(all(double(R.N_starts)==100),'Expected 100 starts per candidate.');
if isfile(checkpoint), delete(checkpoint); end
fprintf('SESSION %s COMPLETE | %d candidates | %.1f min\n',session,height(R),toc(t0)/60);
end

function keys=make_keys(T)
if isempty(T), keys=strings(0,1); return; end
keys=string(double(T.scenario_id))+"|"+string(double(T.replicate))+"|"+string(T.method)+"|"+ ...
    string(double(T.r))+"|"+string(double(T.g))+"|"+compose('%.6f',double(T.alpha));
end

function T=merge_unique(A,B)
if isempty(A), T=B; elseif isempty(B), T=A; else, T=[A;B]; end
if isempty(T), return; end
k=make_keys(T); [~,ia]=unique(k,'stable'); T=T(sort(ia),:);
end
