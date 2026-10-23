function OUT = RUN_TFKM_TRKM_FRESH_INVARIANTS_V1(worker_id,n_workers,documents_root)
%RUN_TFKM_TRKM_FRESH_INVARIANTS_V1
% Truth-open, theorem-oriented invariant audit for the 1230-unit fresh study.
%
% PURPOSE
% -------
% This is NOT a new selector and does NOT refit TFKM/TRKM. It studies
% method-independent quantities that can serve as identification targets for
% contamination mass alpha, cluster count g, latent rank r, and the true
% partition. The only computations are DGP regeneration, oracle geometry,
% low-dimensional PCA/SVD, and split-sample between-scatter diagnostics.
%
% Main invariant families:
%   Q : robust radial-energy separation (diffuse leverage / alpha)
%   P : true-subspace classification margin (partition identifiability)
%   G : underfit merge loss vs Gaussian within-cluster split gain (g)
%   R : rank/eigengap of between scatter; raw vs cross-fitted (r)
%   U : oracle-inlier PCA accessibility of the generating subspace
%
% The cross-fitted between-scatter operator uses independent within-cluster
% half-sample centroids. Conditional on the partition and fixed weights,
% its expectation is the population between-scatter B0. This removes the
% positive high-dimensional mean-noise bias of the raw between scatter.
%
% Usage
% -----
%   RUN_TFKM_TRKM_FRESH_INVARIANTS_V1("SELFTEST",[],[]);
%   RUN_TFKM_TRKM_FRESH_INVARIANTS_V1(1,6);
%   ...
%   RUN_TFKM_TRKM_FRESH_INVARIANTS_V1(6,6);
%
% Then run:
%   MERGE_TFKM_TRKM_FRESH_INVARIANTS_V1;
%
% No method outcome is changed. Truth is already open in this development
% phase; the purpose is to discover stable mathematical targets.

if nargin>=1 && (ischar(worker_id) || isstring(worker_id)) && strcmpi(string(worker_id),"SELFTEST")
    selftest_local();
    OUT=struct('selftest_pass',true);
    return
end
if nargin<1 || isempty(worker_id), worker_id=1; end
if nargin<2 || isempty(n_workers), n_workers=1; end
if nargin<3, documents_root=""; end
worker_id=double(worker_id); n_workers=double(n_workers);
assert(worker_id>=1 && worker_id<=n_workers && worker_id==floor(worker_id) && n_workers==floor(n_workers), ...
    'Require integer 1 <= worker_id <= n_workers.');

[root,outroot,S]=load_study_local(documents_root);
assert(isfield(S,'catalog') && isfield(S,'manifest') && isfield(S,'options'),'study.mat schema mismatch.');
catalog=S.catalog; manifest=S.manifest;
assert(height(catalog)==41 && height(manifest)==1230,'Fresh study dimensions drifted.');
truthroot=fullfile(outroot,'fresh_truth_open_v1');
assert(isfile(fullfile(truthroot,'TRUTH_OPENED_AT.txt')),'Truth-open marker missing.');
lock_vendor_local(root);

outdir=fullfile(outroot,'fresh_invariants_v1');
if ~isfolder(outdir),mkdir(outdir);end
outfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_invariants.csv',worker_id,n_workers));
matfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_checkpoint.mat',worker_id,n_workers));
if isfile(outfile)
    Rows=readtable(outfile,'TextType','string');
else
    Rows=table();
end

unit_index=(1:height(manifest))';
assigned=mod(unit_index-1,n_workers)+1==worker_id;
Jobs=manifest(assigned,:);

CFG=struct();
CFG.cross_splits=12;
CFG.numeric_tol=1e-12;
CFG.pi_const=pi;
CFG.radial_grid=[.025 .05 .075 .10 .125 .15 .20 .25 .30];
CFG.radial_gmm_repro_tol=5e-8;

fprintf('\n============================================================\n');
fprintf('TFKM/TRKM FRESH INVARIANTS V1 | worker %d/%d\n',worker_id,n_workers);
fprintf('THEOREM-ORIENTED | TRUTH OPEN | NO TFKM/TRKM REFITS\n');
fprintf('Jobs: %d | cross-fitted B splits per unit: %d\n',height(Jobs),CFG.cross_splits);
fprintf('============================================================\n\n');

for iu=1:height(Jobs)
    sid=double(Jobs.scenario_id(iu)); rep=double(Jobs.replicate(iu));
    if ~isempty(Rows) && any(double(Rows.scenario_id)==sid & double(Rows.replicate)==rep)
        fprintf('sid=%d rep=%02d already complete\n',sid,rep);
        continue
    end
    ic=find(double(catalog.scenario_id)==sid,1);
    assert(~isempty(ic),'scenario_id %d absent from catalog.',sid);
    c=catalog(ic,:); cfg=table2struct(c);
    data_seed=double(Jobs.data_seed(iu));
    t0=tic;
    sim=generate_fixed_dataset_v1(cfg,data_seed);
    [y,out,U0]=truth_fields_local(sim,double(c.r0));
    X=double(sim.X); in=~out; Xin=X(in,:); yin=y(in);
    assert(size(X,1)==double(c.n) && size(X,2)==double(c.p),'DGP dimension mismatch.');
    assert(sum(in)>double(c.g0)+2,'Too few true inliers.');
    U0=orth_local(U0);
    assert(size(U0,2)==double(c.r0),'U0 rank mismatch.');

    % ------------------------------------------------------------------
    % Q: robust radial-energy invariant and exact frozen radial estimator.
    % ------------------------------------------------------------------
    Q=radial_invariant_local(X,out,double(c.alpha0),CFG);
    [re,ra,rw,rb,rs]=radial_gmm_bic_local(Q.q,CFG);
    [frep,frozen_ok]=frozen_radial_reproduction_local(outroot,sid,rep,re,ra,rw,rb,rs,CFG);

    % ------------------------------------------------------------------
    % P/G/R: oracle inlier geometry in the generating signal subspace.
    % ------------------------------------------------------------------
    z=relabel_local(yin);
    g0=double(c.g0); r0=double(c.r0);
    assert(numel(unique(z))==g0,'Truth group count mismatch.');
    [Ctarget,centroid_source]=target_centroids_local(sim,Xin,z,U0,g0,r0);
    Geo=geometry_invariants_local(Xin,z,U0,Ctarget,double(c.delta),CFG);

    % ------------------------------------------------------------------
    % U: can the generating signal subspace be seen from oracle inliers?
    % ------------------------------------------------------------------
    PCA=pca_accessibility_local(Xin,U0,r0,Geo.beta_min,double(c.rho));

    % ------------------------------------------------------------------
    % R: raw and cross-fitted between-scatter rank diagnostics.
    % Cross-fitting removes same-sample centroid-noise bias in expectation.
    % ------------------------------------------------------------------
    Rank=rank_invariant_local(Xin,z,U0,Ctarget,data_seed,CFG);

    R=empty_row_local();
    R.scenario_id=sid; R.replicate=rep; R.data_seed=data_seed;
    R.panel=string(c.panel); R.contamination=string(c.contamination);
    R.n=double(c.n); R.p=double(c.p); R.n_inlier=sum(in); R.n_outlier=sum(out);
    R.p_over_n=double(c.p_over_n); R.rho=double(c.rho); R.delta=double(c.delta);
    R.r0_POSTHOC=r0; R.g0_POSTHOC=g0; R.alpha0_POSTHOC=double(c.alpha0);
    R.centroid_source=centroid_source;

    R.radial_logq_median_in=Q.med_in; R.radial_logq_sd_in=Q.sd_in; R.radial_logq_mad_in=Q.mad_in;
    R.radial_logq_median_out=Q.med_out; R.radial_truth_gap=Q.gap;
    R.radial_truth_sep_sd=Q.sep_sd; R.radial_truth_sep_mad=Q.sep_mad; R.radial_truth_AUC_out_vs_in=Q.auc;
    R.radial_topalpha_precision=Q.top_precision; R.radial_topalpha_recall=Q.top_recall;
    R.radial_GMM_emit=re; R.radial_GMM_alpha_hat=ra; R.radial_GMM_high_weight=rw;
    R.radial_GMM_delta_BIC=rb; R.radial_GMM_separation_sd=rs;
    R.radial_GMM_alpha_abs_error=abs_nan_local(ra-double(c.alpha0));
    R.radial_frozen_reproduction_pass=frozen_ok;
    R.radial_frozen_source=frep;

    R.centroid_min_euclid=Geo.min_euclid; R.centroid_min_mahal=Geo.min_mahal;
    R.partition_margin_ratio_8logn=Geo.partition_margin_ratio;
    R.oracle_nearest_centroid_ARI=Geo.oracle_ari;
    R.oracle_nearest_centroid_exact=Geo.oracle_exact;
    R.within_signal_lambda_min=Geo.within_lam_min; R.within_signal_lambda_max=Geo.within_lam_max;

    R.g_merge_loss_min=Geo.merge_min; R.g_merge_loss_median=Geo.merge_median;
    R.g_gaussian_split_gain_max=Geo.split_gain_max;
    R.G_g_merge_over_split=Geo.Gg; R.G_g_margin=Geo.merge_min-Geo.split_gain_max;
    R.G_g_spherical_simplex_proxy=Geo.Gg_simplex_proxy;

    R.beta1_between=Geo.beta1; R.beta_min_between=Geo.beta_min; R.beta_condition=Geo.beta_condition;
    R.population_centroid_rank=Geo.centroid_rank;

    R.pca_overlap_trueU=PCA.overlap; R.pca_lambda_r=PCA.lam_r; R.pca_lambda_next=PCA.lam_next;
    R.pca_gap=PCA.gap; R.pca_ratio=PCA.ratio; R.pca_BBP_proxy=PCA.bbp_proxy;

    R.rawB_lambda_r=Rank.raw_lam_r; R.rawB_lambda_next=Rank.raw_lam_next; R.rawB_ratio=Rank.raw_ratio;
    R.crossB_lambda_r_median=Rank.cross_lam_r_med; R.crossB_lambda_r_q10=Rank.cross_lam_r_q10;
    R.crossB_lambda_next_median=Rank.cross_lam_next_med; R.crossB_lambda_next_abs_q90=Rank.cross_lam_next_abs_q90;
    R.crossB_rank_SNR=Rank.cross_snr; R.crossB_gap_median=Rank.cross_gap_med;
    R.crossB_op_error_median=Rank.err_med; R.crossB_op_error_q90=Rank.err_q90;
    R.crossB_Weyl_ratio_beta_over_2err90=Rank.weyl_ratio;
    R.crossB_Weyl_condition_rate=Rank.weyl_pass_rate;
    R.crossB_splits=CFG.cross_splits;

    R.runtime_seconds=toc(t0);
    row=struct2table(R,'AsArray',true);
    if isempty(Rows),Rows=row;else,Rows=[Rows;row];end %#ok<AGROW>
    Rows=sortrows(Rows,{'scenario_id','replicate'});
    writetable(Rows,outfile);
    save(matfile,'Rows','Jobs','CFG','worker_id','n_workers','-v7.3');
    fprintf(['sid=%d rep=%02d %-11s %-8s | Gq=%.2f Gg=%.2f Gp=%.2f ' ...
        'PCAov=%.3f Rweyl=%.2f | %.2fs\n'], ...
        sid,rep,char(string(c.panel)),char(string(c.contamination)), ...
        Q.sep_mad,Geo.Gg,Geo.partition_margin_ratio,PCA.overlap,Rank.weyl_ratio,R.runtime_seconds);
end

fprintf('\nWorker %d/%d complete: %d rows saved.\n%s\n',worker_id,n_workers,height(Rows),outfile);
OUT=struct('rows',Rows,'output_root',string(outdir),'worker_id',worker_id,'n_workers',n_workers);
end

% =========================================================================
function R=empty_row_local()
R=struct();
R.scenario_id=NaN;R.replicate=NaN;R.data_seed=NaN;R.panel="";R.contamination="";
R.n=NaN;R.p=NaN;R.n_inlier=NaN;R.n_outlier=NaN;R.p_over_n=NaN;R.rho=NaN;R.delta=NaN;
R.r0_POSTHOC=NaN;R.g0_POSTHOC=NaN;R.alpha0_POSTHOC=NaN;R.centroid_source="";
R.radial_logq_median_in=NaN;R.radial_logq_sd_in=NaN;R.radial_logq_mad_in=NaN;R.radial_logq_median_out=NaN;
R.radial_truth_gap=NaN;R.radial_truth_sep_sd=NaN;R.radial_truth_sep_mad=NaN;R.radial_truth_AUC_out_vs_in=NaN;
R.radial_topalpha_precision=NaN;R.radial_topalpha_recall=NaN;R.radial_GMM_emit=false;R.radial_GMM_alpha_hat=NaN;
R.radial_GMM_high_weight=NaN;R.radial_GMM_delta_BIC=NaN;R.radial_GMM_separation_sd=NaN;R.radial_GMM_alpha_abs_error=NaN;
R.radial_frozen_reproduction_pass=false;R.radial_frozen_source="";
R.centroid_min_euclid=NaN;R.centroid_min_mahal=NaN;R.partition_margin_ratio_8logn=NaN;
R.oracle_nearest_centroid_ARI=NaN;R.oracle_nearest_centroid_exact=false;
R.within_signal_lambda_min=NaN;R.within_signal_lambda_max=NaN;
R.g_merge_loss_min=NaN;R.g_merge_loss_median=NaN;R.g_gaussian_split_gain_max=NaN;
R.G_g_merge_over_split=NaN;R.G_g_margin=NaN;R.G_g_spherical_simplex_proxy=NaN;
R.beta1_between=NaN;R.beta_min_between=NaN;R.beta_condition=NaN;R.population_centroid_rank=NaN;
R.pca_overlap_trueU=NaN;R.pca_lambda_r=NaN;R.pca_lambda_next=NaN;R.pca_gap=NaN;R.pca_ratio=NaN;R.pca_BBP_proxy=NaN;
R.rawB_lambda_r=NaN;R.rawB_lambda_next=NaN;R.rawB_ratio=NaN;
R.crossB_lambda_r_median=NaN;R.crossB_lambda_r_q10=NaN;R.crossB_lambda_next_median=NaN;R.crossB_lambda_next_abs_q90=NaN;
R.crossB_rank_SNR=NaN;R.crossB_gap_median=NaN;R.crossB_op_error_median=NaN;R.crossB_op_error_q90=NaN;
R.crossB_Weyl_ratio_beta_over_2err90=NaN;R.crossB_Weyl_condition_rate=NaN;R.crossB_splits=NaN;R.runtime_seconds=NaN;
end

% =========================================================================
function Q=radial_invariant_local(X,out,alpha0,CFG)
q=radial_q_local(X);x=log(max(q,eps));in=~out;
Q=struct();Q.q=q;Q.med_in=median(x(in));Q.sd_in=std(x(in));Q.mad_in=1.4826*median(abs(x(in)-median(x(in))));
if ~isfinite(Q.mad_in)||Q.mad_in<=eps,Q.mad_in=max(Q.sd_in,eps);end
Q.med_out=NaN;Q.gap=NaN;Q.sep_sd=NaN;Q.sep_mad=NaN;Q.auc=NaN;Q.top_precision=NaN;Q.top_recall=NaN;
if any(out)
    Q.med_out=median(x(out));Q.gap=Q.med_out-Q.med_in;
    Q.sep_sd=Q.gap/max(Q.sd_in,eps);Q.sep_mad=Q.gap/max(Q.mad_in,eps);
    Q.auc=rank_auc_local(x,double(out));
end
if alpha0>0 && any(out)
    m=max(1,round(alpha0*numel(q)));[~,ord]=sort(q,'descend');sel=false(numel(q),1);sel(ord(1:min(m,numel(q))))=true;
    Q.top_precision=sum(sel&out)/max(sum(sel),1);Q.top_recall=sum(sel&out)/max(sum(out),1);
end
end

function q=radial_q_local(X)
med=median(X,1);sc=median(abs(X-med),1);pos=sc(isfinite(sc)&sc>0);
if isempty(pos),floor_sc=1;else,floor_sc=max(eps,1e-12*median(pos));end
sc(~isfinite(sc)|sc<floor_sc)=floor_sc;Z=(X-med)./sc;q=mean(Z.^2,2);
end

function [emit,aGrid,wHi,dBIC,sep]=radial_gmm_bic_local(q,CFG)
x=log(max(double(q(:)),eps));x=x(isfinite(x));n=numel(x);assert(n>=20,'Too few radial scores.');
mu=mean(x);v=max(mean((x-mu).^2),1e-8);ll1=sum(lognormpdf_local(x,mu,v));bic1=-2*ll1+2*log(n);
starts=[.30 .70 .50;.20 .80 .25;.40 .90 .20];bestLL=-Inf;best=[];
for s=1:size(starts,1)
    m1=quantile_simple_local(x,starts(s,1));m2=quantile_simple_local(x,starts(s,2));vv=max(v,1e-6);v1=vv;v2=vv;p2=starts(s,3);p1=1-p2;last=-Inf;
    for it=1:250
        l1=log(max(p1,eps))+lognormpdf_local(x,m1,v1);l2=log(max(p2,eps))+lognormpdf_local(x,m2,v2);mx=max([l1 l2],[],2);den=mx+log(exp(l1-mx)+exp(l2-mx));
        r2=exp(l2-den);r1=1-r2;n1=sum(r1);n2=sum(r2);if n1<2||n2<2,break,end
        p1=n1/n;p2=n2/n;m1=sum(r1.*x)/n1;m2=sum(r2.*x)/n2;v1=max(sum(r1.*(x-m1).^2)/n1,1e-6);v2=max(sum(r2.*(x-m2).^2)/n2,1e-6);ll=sum(den);
        if abs(ll-last)<=1e-10*(1+abs(ll)),break,end;last=ll;
    end
    l1=log(max(p1,eps))+lognormpdf_local(x,m1,v1);l2=log(max(p2,eps))+lognormpdf_local(x,m2,v2);mx=max([l1 l2],[],2);ll=sum(mx+log(exp(l1-mx)+exp(l2-mx)));
    if ll>bestLL,bestLL=ll;best=[p1 p2 m1 m2 v1 v2];end
end
assert(~isempty(best),'GMM fit failed.');p1=best(1);p2=best(2);m1=best(3);m2=best(4);v1=best(5);v2=best(6);
if m1>m2,tmp=p1;p1=p2;p2=tmp;tmp=m1;m1=m2;m2=tmp;tmp=v1;v1=v2;v2=tmp;end
wHi=p2;bic2=-2*bestLL+5*log(n);dBIC=bic1-bic2;sep=(m2-m1)/sqrt(max(.5*(v1+v2),eps));
emit=isfinite(dBIC)&&dBIC>0&&isfinite(wHi)&&wHi>0&&wHi<=max(CFG.radial_grid)+CFG.numeric_tol;
if emit,aGrid=project_grid_local(wHi,CFG.radial_grid);else,aGrid=NaN;end
end

function [source,pass]=frozen_radial_reproduction_local(outroot,sid,rep,emit,a,w,b,s,CFG)
source="RECOMPUTED_ONLY";pass=true;
f=fullfile(outroot,'blind_finalization_v1','TRKM',sprintf('s%d',sid),sprintf('rep%02d',rep),'blind_final.mat');
if ~isfile(f),return,end
Z=load(f,'Prediction');if ~isfield(Z,'Prediction')||height(Z.Prediction)~=1,return,end
P=Z.Prediction;req={'radial_emit','radial_alpha_hat','radial_high_weight','radial_delta_BIC','radial_separation_sd'};
if ~all(ismember(req,P.Properties.VariableNames)),return,end
pass=(logical(P.radial_emit)==logical(emit));
pass=pass && same_nan_num_local(double(P.radial_alpha_hat),a,CFG.radial_gmm_repro_tol);
pass=pass && same_nan_num_local(double(P.radial_high_weight),w,CFG.radial_gmm_repro_tol);
pass=pass && same_nan_num_local(double(P.radial_delta_BIC),b,CFG.radial_gmm_repro_tol*max(1,abs(b)));
pass=pass && same_nan_num_local(double(P.radial_separation_sd),s,CFG.radial_gmm_repro_tol*max(1,abs(s)));
assert(pass,'Frozen radial reproduction mismatch sid=%d rep=%d.',sid,rep);source="FROZEN_REPRODUCED";
end

% =========================================================================
function [Ctarget,source]=target_centroids_local(sim,Xin,z,U0,g0,r0)
source="REALIZED_SAMPLE_CENTROIDS";
Ctarget=[];
if isfield(sim,'centroids')
    C0=double(sim.centroids);
    if isequal(size(C0),[g0 r0]) && all(isfinite(C0(:)))
        Ctarget=C0;source="DGP_POPULATION_CENTROIDS";
    end
end
if isempty(Ctarget)
    Xc=Xin-mean(Xin,1);Y=Xc*U0;Ctarget=zeros(g0,r0);
    for k=1:g0,Ctarget(k,:)=mean(Y(z==k,:),1);end
end
% Center in the coordinate system used below; pairwise geometry is unchanged.
nk=accumarray(z,1,[g0 1]);piw=nk/sum(nk);Ctarget=Ctarget-(piw'*Ctarget);
end

function G=geometry_invariants_local(Xin,z,U0,Ctarget,delta,CFG)
g=max(z);r=size(U0,2);n=size(Xin,1);Xc=Xin-mean(Xin,1);Y=Xc*U0;
nk=accumarray(z,1,[g 1]);piw=nk/sum(nk);assert(all(nk>=2),'Tiny truth cluster.');
Csamp=zeros(g,r);Wk=cell(g,1);lamMax=zeros(g,1);
for k=1:g
    Z=Y(z==k,:);Csamp(k,:)=mean(Z,1);E=Z-mean(Z,1);W=(E'*E)/size(E,1);W=(W+W')/2;Wk{k}=W;lamMax(k)=max(real(eig(W)));
end
% Align target centroids to realized global centering. Only differences matter.
C=Ctarget-(piw'*Ctarget);
% Pooled within covariance in U0.
Wpool=zeros(r);den=0;
for k=1:g,E=Y(z==k,:)-Csamp(k,:);Wpool=Wpool+E'*E;den=den+size(E,1);end
Wpool=Wpool/max(den,1);Wpool=(Wpool+Wpool')/2;
we=sort(real(eig(Wpool)),'ascend');G.within_lam_min=max(we(1),0);G.within_lam_max=max(we(end),0);
ridge=1e-10*max(trace(Wpool)/max(r,1),1);Wi=pinv(Wpool+ridge*eye(r));
D=[];DM=[];Merge=[];
for a=1:g
    for b=a+1:g
        d=C(a,:)-C(b,:);D(end+1,1)=norm(d);DM(end+1,1)=sqrt(max(d*Wi*d',0)); %#ok<AGROW>
        Merge(end+1,1)=(piw(a)*piw(b)/(piw(a)+piw(b)))*sum(d.^2); %#ok<AGROW>
    end
end
G.min_euclid=min(D);G.min_mahal=min(DM);G.partition_margin_ratio=(G.min_mahal^2)/(8*log(max(n,3)));
% Nearest target centroid partition in true U0.
DD=zeros(n,g);for k=1:g,E=Y-C(k,:);DD(:,k)=sum(E.^2,2);end
[~,pred]=min(DD,[],2);G.oracle_ari=adjusted_rand_index_v1(z,pred);G.oracle_exact=abs(G.oracle_ari-1)<=1e-12;
% g-invariant: cost of forced merge vs Gaussian optimal 2-means split gain.
G.merge_min=min(Merge);G.merge_median=median(Merge);splitG=piw.*((2/CFG.pi_const)*lamMax);G.split_gain_max=max(splitG);
G.Gg=G.merge_min/max(G.split_gain_max,eps);
% Balanced spherical-simplex special-case proxy: pi*delta^2/4 when r=g-1.
if r==g-1,G.Gg_simplex_proxy=CFG.pi_const*delta^2/4;else,G.Gg_simplex_proxy=NaN;end
% Population/target between scatter in U0.
B=zeros(r);for k=1:g,B=B+piw(k)*(C(k,:)'*C(k,:));end;B=(B+B')/2;
beta=sort(max(real(eig(B)),0),'descend');G.beta1=beta(1);G.beta_min=beta(min(r,numel(beta)));G.beta_condition=G.beta_min/max(G.beta1,eps);
tol=max(1e-10*max(beta(1),1),1e-12);G.centroid_rank=sum(beta>tol);G.Btarget=B;G.piw=piw;
end

% =========================================================================
function P=pca_accessibility_local(Xin,U0,r0,beta_min,rho)
Xc=Xin-mean(Xin,1);n=size(Xc,1);[~,S,V]=svd(Xc,'econ');sv=diag(S);lam=(sv.^2)/max(n,1);
P.lam_r=lam(r0);if numel(lam)>=r0+1,P.lam_next=lam(r0+1);else,P.lam_next=0;end
P.gap=P.lam_r-P.lam_next;P.ratio=P.lam_r/max(P.lam_next,eps);Up=V(:,1:r0);s=svd(U0'*Up);s=min(max(real(s),0),1);P.overlap=sum(s.^2)/r0;
% BBP-style proxy for the weakest population spike. This is diagnostic,
% not a theorem claim for the full robust estimator.
gamma=size(Xin,2)/size(Xin,1);spike=(1+max(beta_min,0))/max(rho,eps);P.bbp_proxy=spike/(1+sqrt(gamma));
end

% =========================================================================
function R=rank_invariant_local(Xin,z,U0,Ctarget,data_seed,CFG)
g=max(z);r0=size(U0,2);n=size(Xin,1);p=size(Xin,2);nk=accumarray(z,1,[g 1]);piw=nk/sum(nk);
% Target B0 in the generating subspace, with realized mixture weights.
Ct=Ctarget-(piw'*Ctarget);B0=zeros(r0);for k=1:g,B0=B0+piw(k)*(Ct(k,:)'*Ct(k,:));end;B0=(B0+B0')/2;
beta=sort(max(real(eig(B0)),0),'descend');beta_min=beta(r0);
% Raw same-sample between scatter from truth labels in ambient space.
mu=mean(Xin,1);A=zeros(p,g);
for k=1:g,mk=mean(Xin(z==k,:),1)-mu;A(:,k)=sqrt(piw(k))*mk';end
raw=sort(max(real(eig(A'*A)),0),'descend');raw=pad_local(raw,r0+1,0);R.raw_lam_r=raw(r0);R.raw_lam_next=raw(r0+1);R.raw_ratio=R.raw_lam_r/max(R.raw_lam_next,eps);
% Cross-fitted B: independent half-sample centroids within each true group.
ls=zeros(CFG.cross_splits,1);ln=zeros(CFG.cross_splits,1);er=zeros(CFG.cross_splits,1);gp=zeros(CFG.cross_splits,1);
for s=1:CFG.cross_splits
    rng(data_seed+880000+s,'twister');MA=zeros(g,p);MB=zeros(g,p);
    for k=1:g
        idx=find(z==k);idx=idx(randperm(numel(idx)));na=floor(numel(idx)/2);assert(na>=2 && numel(idx)-na>=2,'Cluster too small for cross-fit.');
        ia=idx(1:na);ib=idx(na+1:end);MA(k,:)=mean(Xin(ia,:),1);MB(k,:)=mean(Xin(ib,:),1);
    end
    ga=piw'*MA;gb=piw'*MB;AA=zeros(p,g);BB=zeros(p,g);
    for k=1:g,AA(:,k)=sqrt(piw(k))*(MA(k,:)-ga)';BB(:,k)=sqrt(piw(k))*(MB(k,:)-gb)';end
    ev=cross_eigs_local(AA,BB);ev=sort(real(ev),'descend');ev=pad_local(ev,r0+1,0);ls(s)=ev(r0);ln(s)=ev(r0+1);gp(s)=ev(r0)-max(ev(r0+1),0);
    er(s)=cross_target_error_local(AA,BB,U0,B0);
end
R.cross_lam_r_med=median(ls);R.cross_lam_r_q10=qlinear_local(ls,.10);R.cross_lam_next_med=median(ln);R.cross_lam_next_abs_q90=qlinear_local(abs(ln),.90);
R.cross_snr=R.cross_lam_r_med/max(R.cross_lam_next_abs_q90,eps);R.cross_gap_med=median(gp);R.err_med=median(er);R.err_q90=qlinear_local(er,.90);
R.weyl_ratio=beta_min/max(2*R.err_q90,eps);R.weyl_pass_rate=mean(er<beta_min/2);
end

function ev=cross_eigs_local(A,B)
g=size(A,2);C=[A B];K=[zeros(g) .5*eye(g);.5*eye(g) zeros(g)];[~,R]=qr(C,0);ev=eig((R*K*R'+(R*K*R')')/2);
end
function e=cross_target_error_local(A,B,U0,B0)
g=size(A,2);r=size(U0,2);C=[A B U0];K=zeros(2*g+r);K(1:g,g+(1:g))=.5*eye(g);K(g+(1:g),1:g)=.5*eye(g);K(2*g+(1:r),2*g+(1:r))=-B0;
[~,R]=qr(C,0);D=R*K*R';D=(D+D')/2;ev=eig(D);e=max(abs(real(ev)));
end

% =========================================================================
function [y,out,U0]=truth_fields_local(sim,r0)
if isfield(sim,'labels'),y=double(sim.labels(:));elseif isfield(sim,'gtrue'),y=double(sim.gtrue(:));else,error('Truth labels missing.');end
if isfield(sim,'outlier'),out=logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),out=logical(sim.outlier_mask(:));else,error('Outlier mask missing.');end
assert(numel(y)==numel(out),'Truth field length mismatch.');assert(isfield(sim,'U')&&size(sim.U,2)>=r0,'Generating U missing.');U0=double(sim.U(:,1:r0));
end
function U=orth_local(U),[U,~]=qr(double(U),0);end
function z=relabel_local(y),[~,~,z]=unique(double(y(:)));z=double(z);end
function x=pad_local(x,n,val),x=x(:);if numel(x)<n,x(end+1:n,1)=val;end,end
function v=abs_nan_local(x),if isfinite(x),v=abs(x);else,v=NaN;end,end
function tf=same_nan_num_local(a,b,tol),if isnan(a)&&isnan(b),tf=true;elseif isfinite(a)&&isfinite(b),tf=abs(a-b)<=tol;else,tf=false;end,end

% =========================================================================
function a=rank_auc_local(score,y)
score=double(score(:));y=double(y(:));ok=isfinite(score)&isfinite(y);score=score(ok);y=y(ok);pos=(y==1);neg=(y==0);np=sum(pos);nn=sum(neg);if np==0||nn==0,a=NaN;return,end
[ss,ord]=sort(score);r=zeros(size(score));i=1;while i<=numel(ss),j=i;while j<numel(ss)&&ss(j+1)==ss(i),j=j+1;end;rr=(i+j)/2;r(ord(i:j))=rr;i=j+1;end;a=(sum(r(pos))-np*(np+1)/2)/(np*nn);
end
function y=lognormpdf_local(x,m,v),y=-.5*(log(2*pi*v)+(x-m).^2./v);end
function q=quantile_simple_local(x,p),q=qlinear_local(x,p);end
function q=qlinear_local(x,p),x=sort(double(x(:)));x=x(isfinite(x));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x(1);return,end;t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end,end
function a=project_grid_local(x,grid),[~,i]=min(abs(grid-x));a=grid(i);end

% =========================================================================
function lock_vendor_local(root)
v=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));assert(isfolder(v),'Canonical vendor root missing: %s',v);addpath(genpath(char(v)),'-begin');rehash path;
req=["generate_fixed_dataset_v1","adjusted_rand_index_v1"];
for q=req,p=string(which(char(q)));assert(strlength(p)>0,'Missing %s.',q);assert(startsWith(p,v,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,p);end
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found. Pass Documents explicitly.');S=load(fullfile(outroot,'study.mat'));
end

% =========================================================================
function selftest_local()
% Radial AUC.
assert(abs(rank_auc_local([0;1;2;3],[0;0;1;1])-1)<=1e-12);
% Merge vs Gaussian split special case: balanced g=4, delta=10, sigma2=1.
g=4;d=10;merge=(1/g*1/g)/(2/g)*d^2;split=(1/g)*(2/pi);assert(abs((merge/split)-(pi*d^2/4))<1e-10);
% Cross-fit compact eigs equals explicit matrix eigs.
rng(7);A=randn(12,4);B=randn(12,4);E=.5*(A*B'+B*A');e1=sort(real(eig(E)),'descend');e2=sort(real(cross_eigs_local(A,B)),'descend');nz=e1(abs(e1)>1e-10);nz2=e2(abs(e2)>1e-10);assert(numel(nz)==numel(nz2)&&max(abs(nz-nz2))<1e-9);
% Cross-target error compact representation.
U=orth_local(randn(12,3));B0=diag([3 2 1]);D=E-U*B0*U';ee=max(abs(eig((D+D')/2)));ec=cross_target_error_local(A,B,U,B0);assert(abs(ee-ec)<1e-9);
fprintf('SELFTEST PASS: radial, g-gap, cross-fitted rank, and compact-spectrum contracts.\n');
end
