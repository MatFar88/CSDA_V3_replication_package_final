function OUT = RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1(worker_id,n_workers,documents_root)
%RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1
% Truth-open development audit of the exact 100-start pool at (r0,g0,alpha0)
% for ALL 150 WEAK_DELTA7 units, both TFKM and TRKM (300 pools total).
%
% This is a POSTHOC REFIT/REPLAY diagnostic on development data. It does not
% alter V1 and must not be reported as holdout validation.
%
% Usage:
%   RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1("SELFTEST",[],[]);
%   RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1(1,6);
%   ...
%   RUN_FRESH_V2_WEAK_TRUESTATE_POOL_AUDIT_V1(6,6);

if nargin>=1 && isstring(worker_id) && strcmpi(worker_id,"SELFTEST")
    assert(mechanism_local(.95,.99)=="SELECTED_GOOD");
    assert(mechanism_local(.50,.99)=="ARBITRATION_FAILURE");
    assert(mechanism_local(.50,.89)=="BASIN_GENERATION_FAILURE");
    fprintf('SELFTEST PASS: weak true-state pool audit contracts.\n');
    OUT=struct('selftest_pass',true);return
end
if nargin<1||isempty(worker_id),worker_id=1;end
if nargin<2||isempty(n_workers),n_workers=1;end
if nargin<3,documents_root="";end
worker_id=double(worker_id);n_workers=double(n_workers);
assert(worker_id>=1&&worker_id<=n_workers&&worker_id==floor(worker_id)&&n_workers==floor(n_workers),'Require integer 1 <= worker_id <= n_workers.');

[root,outroot,S]=load_study_local(documents_root);catalog=S.catalog;manifest=S.manifest;options=S.options;
weak=catalog(string(catalog.panel)=="WEAK_DELTA7",:);assert(height(weak)==5,'Expected 5 weak cells.');
M=manifest(ismember(double(manifest.scenario_id),double(weak.scenario_id)),:);assert(height(M)==150,'Expected 150 weak units.');
lock_vendor_local(root);

names=string(options.hd_init_names(:)');expected=["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"];
assert(isequal(names,expected),'Initialization-family names drifted.');
allocT=double(options.hd_tfkm_allocation(:)');allocR=double(options.hd_trkm_allocation(:)');
assert(isequal(allocT,[50 35 10 5]),'TFKM allocation drifted.');assert(isequal(allocR,[40 30 15 15]),'TRKM allocation drifted.');

hierRoot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
tfCube=fullfile(hierRoot,'fits','TFKM');trCube=fullfile(outroot,'fits','TRKM');
assert(isfolder(tfCube)&&isfolder(trCube),'Frozen candidate roots missing.');
outdiag=fullfile(outroot,'fresh_v2_diagnostic_suite_v1','weak_true_state_pool_audit_v1');if ~isfolder(outdiag),mkdir(outdiag);end
poolfile=fullfile(outdiag,sprintf('worker_%02d_of_%02d_pool.csv',worker_id,n_workers));
famfile=fullfile(outdiag,sprintf('worker_%02d_of_%02d_family.csv',worker_id,n_workers));
if isfile(poolfile),Pool=readtable(poolfile,'TextType','string');else,Pool=table();end
if isfile(famfile),Fam=readtable(famfile,'TextType','string');else,Fam=table();end

Jobs=table();
for im=1:2
    method="TFKM";if im==2,method="TRKM";end
    Q=M;Q.method=repmat(method,height(Q),1);Jobs=[Jobs;Q]; %#ok<AGROW>
end
Jobs=sortrows(Jobs,{'replicate','scenario_id','method'});assigned=mod((1:height(Jobs))'-1,n_workers)+1==worker_id;Jobs=Jobs(assigned,:);

fprintf('\n============================================================\n');
fprintf('WEAK TRUE-STATE 100-START POOL AUDIT | worker %d/%d\n',worker_id,n_workers);
fprintf('TRUTH OPEN | DEVELOPMENT FORENSICS | EXACT FROZEN SEEDS\n');
fprintf('============================================================\n');

for ij=1:height(Jobs)
    sid=double(Jobs.scenario_id(ij));rep=double(Jobs.replicate(ij));method=string(Jobs.method(ij));
    already=~isempty(Pool)&&any(double(Pool.scenario_id)==sid&double(Pool.replicate)==rep&string(Pool.method)==method);
    if already,fprintf('sid=%d rep=%02d %s already complete\n',sid,rep,method);continue,end
    ic=find(double(catalog.scenario_id)==sid,1);assert(~isempty(ic));c=catalog(ic,:);cfg=table2struct(c);
    dataSeed=double(Jobs.data_seed(ij));if method=="TFKM",pairSeed=double(Jobs.tfkm_pair_seed(ij));cubeBase=tfCube;alloc=allocT;core="REG_TFKM";else,pairSeed=double(Jobs.trkm_pair_seed(ij));cubeBase=trCube;alloc=allocR;core="REG_TRKM";end
    sim=generate_fixed_dataset_v1(cfg,dataSeed);X=sim.X;[y,out]=truth_local(sim);in=~out;r0=double(c.r0);g0=double(c.g0);a0=double(c.alpha0);
    cf=fullfile(cubeBase,sprintf('s%d',sid),sprintf('rep%02d',rep),'cube.mat');assert(isfile(cf));Z=load(cf,'Candidates','Labels','Trim');C=Z.Candidates;L=Z.Labels;T=Z.Trim;
    old=C(double(C.r)==r0&double(C.g)==g0&abs(double(C.alpha)-a0)<=1e-12,:);assert(height(old)==1,'Stored true state missing/duplicated.');
    t0=tic;P=generate_pool_local(X,r0,g0,a0,options,pairSeed,names,alloc,core);
    truthARI=start_truth_ari_local(P,y,in);
    if method=="TFKM"
        H=hier4_decompose_local(P,.2,1e-12);chosen=double(H.chosen_index);selectorScore=double(H.cross_H);
        reproduce_tfkm_local(P,H,old,L,T,1e-8,5e-8);
    else
        sc=score_candidates_robust_views_v3(P,X,g0,r0);rv=double(sc.cv_score(:));sel=select_robust_view_v3(P,X,g0,r0,"RV_CV");chosen=double(sel.chosen);selectorScore=rv(chosen);
        reproduce_trkm_local(P,chosen,selectorScore,old,L,T,5e-8,5e-8);
    end
    valid=find(P.valid);vals=truthARI(valid);selected=truthARI(chosen);[best,mx]=max(vals);bestIdx=valid(mx);
    q10=qlinear_local(vals,.10);med=qlinear_local(vals,.50);q90=qlinear_local(vals,.90);q95=qlinear_local(vals,.95);
    n90=sum(vals>=.90-1e-12);n95=sum(vals>=.95-1e-12);n99=sum(vals>=.99-1e-12);rankSel=1+sum(vals>selected+1e-12);mech=mechanism_local(selected,best);
    row=table(method,sid,rep,string(c.panel),double(c.p_over_n),double(c.rho),double(c.delta),r0,g0,a0,dataSeed,pairSeed,numel(valid), ...
        chosen,string(P.family(chosen)),double(P.start_id(chosen)),selectorScore,selected,rankSel,best,string(P.family(bestIdx)),double(P.start_id(bestIdx)), ...
        best-selected,q10,med,q90,q95,n90,n95,n99,n90/numel(valid),n95/numel(valid),n99/numel(valid),mech,toc(t0), ...
        'VariableNames',{'method','scenario_id','replicate','panel','p_over_n','rho','delta','r0_POSTHOC','g0_POSTHOC','alpha0_POSTHOC','data_seed','pair_seed','n_valid_starts', ...
        'selected_pool_index','selected_family','selected_start_id','selected_native_score','selected_truth_ARI_POSTHOC','selected_truth_rank_POSTHOC','best_truth_ARI_POSTHOC','best_family_POSTHOC','best_start_id_POSTHOC', ...
        'best_minus_selected_ARI_POSTHOC','truth_ARI_q10_POSTHOC','truth_ARI_median_POSTHOC','truth_ARI_q90_POSTHOC','truth_ARI_q95_POSTHOC','n_truthlike_090_POSTHOC','n_truthlike_095_POSTHOC','n_truthlike_099_POSTHOC', ...
        'frac_truthlike_090_POSTHOC','frac_truthlike_095_POSTHOC','frac_truthlike_099_POSTHOC','mechanism_POSTHOC','runtime_seconds'});
    Pool=[Pool;row]; %#ok<AGROW>
    for jf=1:numel(names)
        idx=find(P.valid&P.family==names(jf));vv=truthARI(idx);famMax=NaN;famMean=NaN;famMed=NaN;nFam90=0;bestStart=NaN;
        if ~isempty(vv),[famMax,zz]=max(vv);famMean=mean(vv);famMed=median(vv);nFam90=sum(vv>=.90-1e-12);bestStart=P.start_id(idx(zz));end
        Fam=[Fam;table(method,sid,rep,r0,g0,a0,jf,names(jf),numel(idx),famMax,famMean,famMed,nFam90,bestStart,chosen>=1&&P.family(chosen)==names(jf), ...
            'VariableNames',{'method','scenario_id','replicate','r0_POSTHOC','g0_POSTHOC','alpha0_POSTHOC','family_index','family','n_valid_family','family_max_truth_ARI_POSTHOC','family_mean_truth_ARI_POSTHOC','family_median_truth_ARI_POSTHOC','family_n_truthlike_090_POSTHOC','family_best_start_id_POSTHOC','selected_family_flag'})]; %#ok<AGROW>
    end
    Pool=sortrows(Pool,{'method','scenario_id','replicate'});Fam=sortrows(Fam,{'method','scenario_id','replicate','family_index'});writetable(Pool,poolfile);writetable(Fam,famfile);
    fprintf('sid=%d rep=%02d %s | selected=%.4f best=%.4f n90=%d/100 | %s | %.1fs\n',sid,rep,method,selected,best,n90,mech,row.runtime_seconds);
end
fprintf('Worker %d/%d complete | %d pool rows.\n',worker_id,n_workers,height(Pool));
OUT=struct('pool',Pool,'family',Fam,'output_root',string(outdiag));
end

% =========================================================================
function P=generate_pool_local(X,r,g,a,options,pairSeed,names,alloc,core)
N=sum(alloc);n=size(X,1);valid=false(N,1);labels=zeros(n,N,'uint16');trim=false(n,N);objective=nan(N,1);family=strings(N,1);start_id=nan(N,1);start_seed=nan(N,1);z=0;
for j=1:numel(names)
    for h=1:alloc(j)
        z=z+1;ss=pairSeed+1000*j+h;family(z)=names(j);start_id(z)=h;start_seed(z)=ss;
        try
            rng(ss,'twister');init=make_robust_initialization_v1(X,g,r,a,names(j),ss);q=regularized_trimmed_fit_corrected_v1(X,init,r,g,a,core,options.max_iterations,options.tauW,options.tauT);
            ok=all(isfinite(q.labels))&&isfinite(q.objective)&&size(q.U,2)==r;valid(z)=ok;if ok,labels(:,z)=uint16(q.labels(:));trim(:,z)=logical(q.trim(:));objective(z)=double(q.objective);end
        catch,valid(z)=false;end
    end
end
assert(any(valid),'No valid starts at true state.');
P=struct('valid',valid,'labels',labels,'trim',trim,'objective',objective,'family',family,'start_id',start_id,'start_seed',start_seed);
if core=="REG_TRKM"
    valid_idx=find(valid);[~,ord]=sort(objective(valid_idx),'ascend');ordered=valid_idx(ord);objective_rank=nan(N,1);objective_rank(ordered)=(1:numel(ordered))';P.objective_rank=objective_rank;
    if isfield(options,'trim_stability_weight'),w=options.trim_stability_weight;else,w=.2;end;P.trim_stability_weight=w;P.current_candidate_index=legacy_candidate_local(P,valid_idx,w);
end
end

function H=hier4_decompose_local(P,w,tol)
valid=logical(P.valid(:));fam=string(P.family(:));obj=double(P.objective(:));N=numel(valid);fn=["WINSOR_PCA_TKM";"SIGN_PCA_TKM";"TKM_ORIGINAL";"RP_TKM"];
within=nan(N,1);fw=nan(numel(fn),1);
for f=1:numel(fn)
    idx=find(valid&fam==fn(f));if isempty(idx),continue,end;M=pair_H_local(P,idx,w);s=mean(M,2,'omitnan');within(idx)=s;mx=max(s,[],'omitnan');cand=idx(abs(s-mx)<=tol);if numel(cand)>1,[~,m]=min(obj(cand));cand=cand(m);end;fw(f)=cand;
end
win=double(fw(isfinite(fw)));assert(~isempty(win));final=nan(N,1);
if numel(win)==1,chosen=win(1);final(chosen)=within(chosen);else,M=pair_H_local(P,win,w);s=mean(M,2,'omitnan');final(win)=s;mx=max(s,[],'omitnan');cand=win(abs(s-mx)<=tol);if numel(cand)>1,[~,m]=min(obj(cand));chosen=cand(m);else,chosen=cand;end,end
H=struct('chosen_index',chosen,'cross_H',final(chosen),'within_H',within(chosen),'objective',obj(chosen),'family_winner_idx',fw);
end
function M=pair_H_local(P,idx,w)
idx=double(idx(:));m=numel(idx);M=zeros(m);
for a=1:m
    for b=a+1:m
        ia=idx(a);ib=idx(b);ta=logical(P.trim(:,ia));tb=logical(P.trim(:,ib));keep=~ta&~tb;ar=0;if sum(keep)>1,ar=adjusted_rand_index_v1(double(P.labels(keep,ia)),double(P.labels(keep,ib)));if ~isfinite(ar),ar=0;end,end;un=ta|tb;if any(un),jac=sum(ta&tb)/sum(un);else,jac=1;end;M(a,b)=ar+w*jac;M(b,a)=M(a,b);
    end
end
end
function chosen=legacy_candidate_local(P,valid_idx,w)
obj=P.objective(valid_idx);[~,ord]=sort(obj,'ascend');ordered=valid_idx(ord);top=ordered(1:min(numel(ordered),max(3,ceil(numel(ordered)/4))));if numel(top)==1,chosen=top(1);return,end
S=zeros(numel(top));for a=1:numel(top),for b=a+1:numel(top),ta=P.trim(:,top(a));tb=P.trim(:,top(b));in=~ta&~tb;ar=0;if sum(in)>1,ar=adjusted_rand_index_v1(double(P.labels(in,top(a))),double(P.labels(in,top(b))));if ~isfinite(ar),ar=0;end,end;un=ta|tb;if any(un),j=sum(ta&tb)/sum(un);else,j=1;end;S(a,b)=ar+w*j;S(b,a)=S(a,b);end,end;[~,m]=max(mean(S,2));chosen=top(m);
end
function a=start_truth_ari_local(P,y,in)
a=nan(numel(P.valid),1);for i=find(P.valid(:))',a(i)=adjusted_rand_index_v1(y(in),double(P.labels(in,i)));end
end
function reproduce_tfkm_local(P,H,old,L,T,scoreTol,objTol)
assert(abs(double(H.cross_H)-double(old.selector_score))<=scoreTol,'TFKM HIER4 score reproduction failure.');assert(string(P.family(H.chosen_index))==string(old.selected_init_family),'TFKM family reproduction failure.');
if ismember('selected_start_id',old.Properties.VariableNames)&&isfinite(double(old.selected_start_id)),assert(double(P.start_id(H.chosen_index))==double(old.selected_start_id),'TFKM start-id reproduction failure.');end
assert(abs(double(P.objective(H.chosen_index))-double(old.objective))<=objTol*max(1,abs(double(old.objective))),'TFKM objective reproduction failure.');idx=double(old.row_index);tr=logical(P.trim(:,H.chosen_index));assert(isequal(tr,logical(T(:,idx))),'TFKM trim reproduction failure.');keep=~tr;ar=adjusted_rand_index_v1(double(P.labels(keep,H.chosen_index)),double(L(keep,idx)));assert(isfinite(ar)&&abs(ar-1)<=1e-12,'TFKM labels reproduction failure.');
end
function reproduce_trkm_local(P,chosen,score,old,L,T,scoreTol,objTol)
assert(abs(score-double(old.selector_cv_score))<=scoreTol,'TRKM RV score reproduction failure.');assert(string(P.family(chosen))==string(old.selected_init_family),'TRKM family reproduction failure.');
if ismember('selected_start_id',old.Properties.VariableNames)&&isfinite(double(old.selected_start_id)),assert(double(P.start_id(chosen))==double(old.selected_start_id),'TRKM start-id reproduction failure.');end
assert(abs(double(P.objective(chosen))-double(old.objective))<=objTol*max(1,abs(double(old.objective))),'TRKM objective reproduction failure.');idx=double(old.row_index);tr=logical(P.trim(:,chosen));assert(isequal(tr,logical(T(:,idx))),'TRKM trim reproduction failure.');keep=~tr;ar=adjusted_rand_index_v1(double(P.labels(keep,chosen)),double(L(keep,idx)));assert(isfinite(ar)&&abs(ar-1)<=1e-12,'TRKM labels reproduction failure.');
end
function s=mechanism_local(selected,best)
if isfinite(selected)&&selected>=.90-1e-12,s="SELECTED_GOOD";elseif isfinite(best)&&best>=.90-1e-12,s="ARBITRATION_FAILURE";else,s="BASIN_GENERATION_FAILURE";end
end
function q=qlinear_local(x,p),x=sort(double(x(:)));x=x(isfinite(x));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x(1);return,end;t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end,end
function [y,out]=truth_local(sim),if isfield(sim,'labels'),y=double(sim.labels(:));elseif isfield(sim,'gtrue'),y=double(sim.gtrue(:));else,error('Truth labels missing.');end;if isfield(sim,'outlier'),out=logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),out=logical(sim.outlier_mask(:));else,error('Outlier mask missing.');end,end
function lock_vendor_local(root)
v=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));assert(isfolder(v));addpath(genpath(char(v)),'-begin');rehash path;req=["generate_fixed_dataset_v1","make_robust_initialization_v1","regularized_trimmed_fit_corrected_v1","adjusted_rand_index_v1","score_candidates_robust_views_v3","select_robust_view_v3"];
for q=req,p=string(which(char(q)));assert(strlength(p)>0,'Missing %s.',q);assert(startsWith(p,v,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,p);end
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end;assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
