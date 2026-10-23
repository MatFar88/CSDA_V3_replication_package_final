function run_final_fresh_forensics_session_v1(session,output_root,mode)
%RUN_FINAL_FRESH_FORENSICS_SESSION_V1 Truth-open explanation layer after frozen selection.
% Default target: every structurally overselected unit OR ARI<.95.
% No diagnostic feeds back into selection.
if nargin<3,mode="targeted";end
pkg=fileparts(mfilename('fullpath'));addpath(genpath(pkg),'-begin');output_root=char(string(output_root));mode=lower(string(mode));
S=load(fullfile(output_root,'study.mat'));study=S.study;options=S.options;catalog=S.catalog;manifest=S.manifest;
selfile=fullfile(output_root,'summary','final_selections.csv');assert(isfile(selfile),'Run merge_final_fresh_validation_v1 first.');Sel=readtable(selfile,'TextType','string');
session=upper(string(session));m=find(string(manifest.session)==session,1);assert(~isempty(m),'Unknown session %s.',session);rep=double(manifest.replicate(m));
Z=Sel(double(Sel.replicate)==rep,:);if mode=="targeted",Z=Z(logical(Z.structural_over)|logical(Z.low_ari_lt_095),:);elseif mode~="all",error('Mode must be targeted or all.');end
if isempty(Z), fid=fopen(fullfile(output_root,'forensics',sprintf('forensic_empty_%s.txt',char(session))),'w');fprintf(fid,'No targeted forensic units for session %s.\n',session);fclose(fid);fprintf('Forensic session %s has zero targets.\n',session);return;end
ufile=fullfile(output_root,'forensics',sprintf('forensic_units_%s.csv',char(session)));cfile=fullfile(output_root,'forensics',sprintf('forensic_clusters_%s.csv',char(session)));dfile=fullfile(output_root,'forensics',sprintf('forensic_directions_%s.csv',char(session)));xfile=fullfile(output_root,'forensics',sprintf('forensic_contingency_%s.csv',char(session)));
if isfile(ufile),E=readtable(ufile,'TextType','string');if height(E)==height(Z),fprintf('Forensic session %s already complete (%d targets).\n',session,height(E));return;end,end
Urows=table();Crows=table();Drows=table();Xrows=table();
fprintf('\n============================================================\nFORENSICS %s | rep=%d | targets=%d | mode=%s\nSelection is frozen; truth-open diagnostics only.\n============================================================\n',session,rep,height(Z),mode);
for i=1:height(Z)
    sr=Z(i,:);sid=double(sr.scenario_id);ic=find(double(catalog.scenario_id)==sid,1);cfg=table2struct(catalog(ic,:));method=string(sr.method);r=double(sr.selected_r);g=double(sr.selected_g);a=double(sr.selected_alpha);
    seed=double(study.validation_base_seed)+sid*100000+rep;if method=="REG_TFKM_HD",pair_seed=seed+110000;else,pair_seed=seed+220000;end
    sim=generate_fixed_dataset_v1(cfg,seed);t0=tic;fit=run_hd_multistart_top_v1(sim.X,r,g,a,method,options,pair_seed);assert(fit.valid,'Forensic refit invalid sid=%d rep=%d %s.',sid,rep,method);fit=normalize_fit_contract_v1(fit,size(sim.X,1),g,a,method);
    in=truth_inlier_mask(sim);ari_now=ari_local(double(sim.labels(in)),double(fit.labels(in)));adiff=abs(ari_now-double(sr.ari_inlier));assert(adiff<=1e-8,'Forensic reproduction mismatch %.3g.',adiff);
    [fm,CR,DR,XR]=forensic_metrics(sim,fit,cfg,sid,rep,method,seed,g,double(study.random_reference_B),double(study.noise_pc_count));
    sig=failure_signature(fm,r,g,double(cfg.r0),double(cfg.c0),ari_now,double(study.forensic_ari_threshold));
    row=table(sid,rep,method,double(cfg.r0),double(cfg.c0),double(cfg.alpha0),r,g,a,double(sr.ari_inlier),ari_now,adiff,toc(t0), ...
        logical(sr.structural_over),logical(sr.low_ari_lt_095),fm.sample_cluster_cert,fm.sample_rank_cert,fm.subspace_coverage,fm.estimate_truth_fraction,fm.min_principal_cos,fm.extra_dim, ...
        fm.n_pred_labels_on_inliers,fm.fitted_cluster_purity,fm.ari_after_truth_open_merge,fm.merge_ari_gain,fm.extra_labels_only_noninliers,fm.pure_fragmentation,fm.n_extra_clusters,fm.n_extra_clusters_with_inliers,fm.n_extra_clusters_outlier_only,string(sig), ...
        'VariableNames',{'scenario_id','replicate','method','r0','g0','alpha0','selected_r','selected_g','selected_alpha','stored_ari','recomputed_ari','ari_abs_diff','runtime_seconds', ...
        'structural_over','low_ari_lt_095','validation_sample_cluster_cert','validation_sample_rank_cert','subspace_coverage','estimate_truth_fraction','min_principal_cos','orthogonal_extra_dim', ...
        'n_pred_labels_on_inliers','fitted_cluster_purity','ari_after_truth_open_merge','merge_ari_gain','extra_labels_only_noninliers','pure_fragmentation','n_extra_clusters','n_extra_clusters_with_inliers','n_extra_clusters_outlier_only','failure_signature'});
    Urows=[Urows;row];Crows=[Crows;CR];Drows=[Drows;DR];Xrows=[Xrows;XR]; %#ok<AGROW>
    fprintf('  sid=%d %-13s final=(%d,%d,%.3f) ARI=%.3f merge=%.3f cover=%.3f | %s\n',sid,method,r,g,a,ari_now,fm.ari_after_truth_open_merge,fm.subspace_coverage,sig);
end
writetable(Urows,ufile);writetable(Crows,cfile);if ~isempty(Drows),writetable(Drows,dfile);end;writetable(Xrows,xfile);
fprintf('FORENSIC session %s COMPLETE.\n',session);
end

function [fm,CR,DR,XR]=forensic_metrics(sim,fit,cfg,sid,rep,method,seed,gh,B,Lmax)
Xall=double(sim.X);zall=double(sim.labels(:));out=~truth_inlier_mask(sim);in=~out;pred=double(fit.labels(:));trim=logical(fit.trim(:));U0=true_U(sim);[Uh,~]=qr(double(fit.U),0);r0=size(U0,2);rh=size(Uh,2);g0=double(cfg.c0);gh=double(gh);X=Xall(in,:);z=zall(in);zp=pred(in);Xc=X-mean(X,1);Y0=Xc*U0;Yh=Xc*Uh;
M=U0'*Uh;s=svd(M);fm.subspace_coverage=sum(s.^2)/max(r0,1);fm.estimate_truth_fraction=sum(s.^2)/max(rh,1);if isempty(s),fm.min_principal_cos=NaN;else,fm.min_principal_cos=min(s);end
% Sample DGP sanity is reported, never used to exclude validation observations.
fm.sample_cluster_cert=sample_cluster_cert(sim,cfg,seed);fm.sample_rank_cert=sample_rank_cert(sim,cfg,seed,B,Lmax);
% Contingency on true inliers, excluding label 0 from fitted-cluster columns.
P=zeros(g0,gh);for k=1:g0,for h=1:gh,P(k,h)=sum(in & zall==k & pred==h);end,end
[truth_to_fit,fit_to_truth,is_core]=max_overlap_core(P); %#ok<ASGLU>
XR=table();for k=1:g0,for h=1:gh,XR=[XR;table(sid,rep,string(method),k,h,P(k,h),'VariableNames',{'scenario_id','replicate','method','true_group','fitted_cluster','count_inliers'})];end,end %#ok<AGROW>
active=unique(zp(zp>=1));fm.n_pred_labels_on_inliers=numel(active);fm.fitted_cluster_purity=sum(max(P,[],1))/max(sum(P(:)),1);
map=zeros(gh,1);for h=1:gh,[mx,k]=max(P(:,h));if mx>0,map(h)=k;end,end
zm=zeros(size(z));for i=1:numel(z),if zp(i)>=1&&zp(i)<=gh,zm(i)=map(zp(i));end,end;fm.ari_after_truth_open_merge=ari_local(z,zm);fm.merge_ari_gain=fm.ari_after_truth_open_merge-ari_local(z,zp);
fm.extra_labels_only_noninliers=(gh>g0 && fm.n_pred_labels_on_inliers<=g0);fm.pure_fragmentation=(gh>g0 && fm.n_pred_labels_on_inliers>g0 && fm.ari_after_truth_open_merge>=.95 && fm.fitted_cluster_purity>=.95);
% Cluster-level descriptions. Extra clusters are unmatched by the maximum-overlap one-to-one core assignment.
CR=table();wcss=partition_wcss(Yh,zp);Ctrue=true_centroids(sim,cfg,z,Y0);dtrue=minimum_center_distance(Ctrue);nextra=0;nextra_in=0;nextra_outonly=0;
for h=1:gh
    qin=in & pred==h;qout=out & pred==h;qret=~trim & pred==h;nin=sum(qin);nout=sum(qout);nret=sum(qret);nretin=sum(qret&in);nretout=sum(qret&out);[mx,kdom]=max(P(:,h));if mx==0,kdom=0;pur=NaN;else,pur=mx/max(nin,1);end
    core=is_core(h);extra=~core && gh>g0;if extra,nextra=nextra+1;if nin>0,nextra_in=nextra_in+1;else,nextra_outonly=nextra_outonly+1;end,end
    paired=0;du0=NaN;duh=NaN;mergerel=NaN;ss=NaN;spanratio=NaN;
    if extra && kdom>0 && truth_to_fit(kdom)>0
        paired=truth_to_fit(kdom);qa=in & zall==kdom & pred==h;qb=in & zall==kdom & pred==paired;na=sum(qa);nb=sum(qb);
        if na>0&&nb>0
            ma0=mean((Xall(qa,:)-mean(X,1))*U0,1);mb0=mean((Xall(qb,:)-mean(X,1))*U0,1);du0=norm(ma0-mb0);spanratio=du0/max(dtrue,eps);
            mah=mean((Xall(qa,:)-mean(X,1))*Uh,1);mbh=mean((Xall(qb,:)-mean(X,1))*Uh,1);duh=norm(mah-mbh);mergerel=(na*nb/(na+nb)*sum((mah-mbh).^2))/max(wcss,eps);
            qu=in & zall==kdom & ismember(pred,[h paired]);if sum(qu)>=3&&numel(unique(pred(qu)))==2,ss=silhouette_mean((Xall(qu,:)-mean(X,1))*U0,pred(qu));end
        end
    end
    sig=cluster_signature(extra,nin,nout,pur,paired);
    CR=[CR;table(sid,rep,string(method),h,logical(core),logical(extra),fit_to_truth(h),kdom,pur,nin,nout,nret,nretin,nretout,paired,du0,duh,spanratio,mergerel,ss,string(sig), ...
        'VariableNames',{'scenario_id','replicate','method','fitted_cluster','is_core_match','is_extra_cluster','core_assigned_true_group','dominant_true_group','inlier_purity','n_inliers','n_outliers','n_retained','n_retained_inliers','n_retained_outliers','paired_core_cluster','centroid_distance_u0','centroid_distance_uhat','centroid_distance_u0_over_true_minsep','merge_delta_rel_wcss','split_silhouette_u0','cluster_signature'})]; %#ok<AGROW>
end
fm.n_extra_clusters=nextra;fm.n_extra_clusters_with_inliers=nextra_in;fm.n_extra_clusters_outlier_only=nextra_outonly;
% Exact directions in span(Uhat) orthogonal to U0.
if rh>r0,[~,S,V]=svd(M);ds=diag(S);tol=max(size(M))*eps(max([ds;1]));rk=sum(ds>tol);N=V(:,rk+1:end);Ue=Uh*N;else,Ue=zeros(size(Uh,1),0);end
fm.extra_dim=size(Ue,2);DR=table();
% One random-direction reference per dataset.
rng(seed+930000,'twister');Vb=randn(size(Xall,2),B);Vb=Vb-U0*(U0'*Vb);nn=sqrt(sum(Vb.^2,1));good=nn>sqrt(eps);Vb(:,good)=Vb(:,good)./nn(good);Vb(:,~good)=0;Sc=Xc*Vb;rv=nan(B,1);for b=1:B,rv(b)=label_R2(Sc(:,b),z);end;rq95=quantile_local(rv,.95);
R=Xc-(Xc*U0)*U0';Vr=[];if norm(R,'fro')>0,[~,~,VV]=svd(R,'econ');Vr=VV(:,1:min([Lmax,size(VV,2)]));end
Xallc=Xall-mean(Xall,1);
for j=1:size(Ue,2)
    u=Ue(:,j);score=Xc*u;tr2=label_R2(score,z);fr2=label_R2(score,zp);if any(out)&&any(in),or2=label_R2(Xallc*u,double(out));else,or2=NaN;end
    varfrac=sum(score.^2)/max(sum(Yh(:).^2),eps);if isempty(Vr),pcos=0;else,pcos=max(abs(Vr'*u));end;noise_like=tr2<=rq95;
    if isfinite(or2)&&or2>rq95&&or2>=tr2,assoc="OUTLIER_STATUS";elseif tr2>rq95,assoc="TRUE_LABEL_SAMPLE_ASSOC";else,assoc="RANDOMLIKE";end
    DR=[DR;table(sid,rep,string(method),j,tr2,fr2,fr2-tr2,or2,varfrac,rq95,tr2/rq95,logical(noise_like),pcos,string(assoc), ...
        'VariableNames',{'scenario_id','replicate','method','extra_direction','true_label_R2','fitted_label_R2','fit_minus_true_R2','outlier_status_R2','variance_fraction_of_estimated_scores','random_true_R2_q95','true_R2_over_random_q95','extra_noise_like','max_cosine_top_noise_pcs','dominant_external_association'})]; %#ok<AGROW>
end
end

function sig=failure_signature(fm,r,g,r0,g0,ari,thr)
flags=strings(0,1);if ~fm.sample_cluster_cert,flags(end+1)="VALIDATION_DRAW_DGP_CLUSTER_WARNING";end;if ~fm.sample_rank_cert,flags(end+1)="VALIDATION_DRAW_DGP_RANK_WARNING";end;if isfinite(fm.subspace_coverage)&&fm.subspace_coverage<.90,flags(end+1)="FIT_SUBSPACE_DISTORTION";end
if r>r0,flags(end+1)="OVER_R";end;if g>g0&&fm.extra_labels_only_noninliers,flags(end+1)="OVER_G_NONINLIER_ONLY";elseif g>g0&&fm.pure_fragmentation,flags(end+1)="OVER_G_PURE_FRAGMENTATION";elseif g>g0,flags(end+1)="OVER_G_MIXED";end
if ari<thr&&fm.ari_after_truth_open_merge>=thr,flags(end+1)="LOW_ARI_RECOVERABLE_BY_TRUTH_OPEN_MERGE";elseif ari<thr,flags(end+1)="LOW_ARI_CROSS_TRUTH_RESIDUAL";end
if isempty(flags),flags="BENIGN_OR_EXACT";end;sig=strjoin(unique(flags,'stable'),"|");
end
function sig=cluster_signature(extra,nin,nout,purity,paired)
if ~extra,sig="CORE_MATCH";return;end;if nin==0,sig="OUTLIER_ONLY";elseif nout>nin&&isfinite(purity)&&purity>=.8,sig="CONTAMINATION_DOMINANT";elseif paired>0&&isfinite(purity)&&purity>=.9,sig="PURE_FRAGMENT";elseif paired>0&&isfinite(purity)&&purity>=.7,sig="IMPURE_FRAGMENT";else,sig="MIXED_CROSS_TRUTH";end
end

function [truth_to_fit,fit_to_truth,is_core]=max_overlap_core(P)
g0=size(P,1);gh=size(P,2);truth_to_fit=zeros(g0,1);fit_to_truth=zeros(gh,1);best=-Inf;bestvec=[];
if gh>=g0
    SS=nchoosek(1:gh,g0);for i=1:size(SS,1),PP=perms(SS(i,:));for j=1:size(PP,1),v=PP(j,:);sc=0;for k=1:g0,sc=sc+P(k,v(k));end;if sc>best||(sc==best&&lexless(v,bestvec)),best=sc;bestvec=v;end,end,end
    truth_to_fit=bestvec(:);for k=1:g0,fit_to_truth(bestvec(k))=k;end
else
    SS=nchoosek(1:g0,gh);for i=1:size(SS,1),PP=perms(SS(i,:));for j=1:size(PP,1),v=PP(j,:);sc=0;for h=1:gh,sc=sc+P(v(h),h);end;if sc>best||(sc==best&&lexless(v,bestvec)),best=sc;bestvec=v;end,end,end
    for h=1:gh,fit_to_truth(h)=bestvec(h);truth_to_fit(bestvec(h))=h;end
end
is_core=false(gh,1);is_core(truth_to_fit(truth_to_fit>0))=true;
end
function tf=lexless(a,b),if isempty(b),tf=true;return;end;tf=false;for i=1:numel(a),if a(i)<b(i),tf=true;return;elseif a(i)>b(i),return;end,end,end
function cert=sample_cluster_cert(sim,cfg,seed)
in=truth_inlier_mask(sim);X=double(sim.X(in,:));z=double(sim.labels(in));U0=true_U(sim);Xc=X-mean(X,1);Y=Xc*U0;C=true_centroids(sim,cfg,z,Y);p=accumarray(z,1,[double(cfg.c0) 1]);p=p/sum(p);C=C-(p'*C);oracle=ari_local(z,nearest_center(Y,C));Gu=kgrid(Y,z,seed+710000);R=Xc-(Xc*U0)*U0';if norm(R,'fro')>0,[~,~,V]=svd(R,'econ');zz=R*V(:,1);else,zz=zeros(size(X,1),1);end;Ga=kgrid([Y zz],z,seed+720000);g0=double(cfg.c0);[smu,cmu]=margins(Gu,g0);[sma,cma]=margins(Ga,g0);cert=oracle>=.99&&value_g(Gu,g0,"ari")>=.99&&value_g(Ga,g0,"ari")>=.99&&best_g(Gu,"sil")==g0&&best_g(Gu,"ch")==g0&&best_g(Ga,"sil")==g0&&best_g(Ga,"ch")==g0&&smu>1e-6&&cmu>1e-6&&sma>1e-6&&cma>1e-6;
end
function cert=sample_rank_cert(sim,cfg,seed,B,Lmax)
in=truth_inlier_mask(sim);X=double(sim.X(in,:));z=double(sim.labels(in));Xc=X-mean(X,1);U0=true_U(sim);r0=size(U0,2);Y=Xc*U0;g0=double(cfg.c0);[~,pi]=group_weights(z,g0);C=true_centroids(sim,cfg,z,Y);C=C-(pi'*C);W=diag(sqrt(pi))*C;sv=svd(W,'econ');tol=max(size(W))*eps(max([sv;1]));cr=sum(sv>tol);Bpop=C'*diag(pi)*C;Bpop=(Bpop+Bpop')/2;[V,D]=eig(Bpop);[~,ord]=sort(real(diag(D)),'descend');V=real(V(:,ord(1:r0)));tr=zeros(r0,1);for j=1:r0,tr(j)=label_R2(Y*V(:,j),z);end;weak=min(tr);R=Xc-(Xc*U0)*U0';npc=0;if norm(R,'fro')>0,[~,~,Vr]=svd(R,'econ');L=min([Lmax,size(Vr,2)]);q=zeros(L,1);for j=1:L,q(j)=label_R2(Xc*Vr(:,j),z);end;npc=max(q);end;rng(seed+940000,'twister');Vb=randn(size(X,2),B);Vb=Vb-U0*(U0'*Vb);nn=sqrt(sum(Vb.^2,1));gd=nn>sqrt(eps);Vb(:,gd)=Vb(:,gd)./nn(gd);Vb(:,~gd)=0;rv=zeros(B,1);Sc=Xc*Vb;for b=1:B,rv(b)=label_R2(Sc(:,b),z);end;cert=(cr==r0)&&weak>max(npc,quantile_local(rv,.95))+1e-12;
end
function G=kgrid(Y,z,seed),G=table();for g=2:min(6,size(Y,1)-1),rng(seed+g,'twister');try,[idx,~,sumd]=kmeans(Y,g,'Replicates',15,'MaxIter',300,'Start','plus','Display','off');w=sum(sumd);si=silhouette_mean(Y,idx);ch=ch_index(Y,idx,w);ar=ari_local(z,idx);catch,w=NaN;si=NaN;ch=NaN;ar=NaN;end;G=[G;table(g,ar,si,ch,w,'VariableNames',{'g','ari','sil','ch','wcss'})];end,end
function [sm,cm]=margins(G,g0),i=double(G.g)==g0;s0=double(G.sil(i));c0=double(G.ch(i));so=double(G.sil(~i));co=double(G.ch(~i));so=so(isfinite(so));co=co(isfinite(co));if isempty(so),sm=NaN;else,sm=s0-max(so);end;if isempty(co),cm=NaN;else,cm=(c0-max(co))/max(abs(c0),1);end,end
function g=best_g(G,var),x=double(G.(char(var)));x(~isfinite(x))=-Inf;mx=max(x);cand=double(G.g(abs(x-mx)<=1e-12));g=min(cand);end
function v=value_g(G,g,var),i=find(double(G.g)==g,1);if isempty(i),v=NaN;else,v=double(G.(char(var))(i));end,end
function [cnt,p]=group_weights(z,g),cnt=accumarray(z,1,[g 1]);p=cnt/sum(cnt);end
function in=truth_inlier_mask(sim),n=size(sim.X,1);if isfield(sim,'outlier'),in=~logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),in=~logical(sim.outlier_mask(:));else,in=true(n,1);end,end
function U=true_U(sim),if isfield(sim,'U'),U=double(sim.U);elseif isfield(sim,'U_true'),U=double(sim.U_true);else,error('Missing U0.');end;[U,~]=qr(U,0);end
function C=true_centroids(sim,cfg,z,Y),if isfield(sim,'centroids')&&size(sim.centroids,1)==double(cfg.c0)&&size(sim.centroids,2)==size(Y,2),C=double(sim.centroids);else,C=zeros(double(cfg.c0),size(Y,2));for k=1:double(cfg.c0),C(k,:)=mean(Y(z==k,:),1);end,end,end
function d=minimum_center_distance(C),D=center_dist(C);v=D(triu(true(size(D)),1));v=v(v>0);if isempty(v),d=NaN;else,d=min(v);end,end
function D=center_dist(C),ss=sum(C.^2,2);D=sqrt(max(0,ss+ss'-2*(C*C')));D(1:size(D,1)+1:end)=0;end
function lab=nearest_center(Y,C),D=zeros(size(Y,1),size(C,1));for k=1:size(C,1),Z=Y-C(k,:);D(:,k)=sum(Z.^2,2);end;[~,lab]=min(D,[],2);end
function r=label_R2(Y,z),if isvector(Y),Y=Y(:);end;mu=mean(Y,1);tot=sum((Y-mu).^2,'all');if tot<=eps,r=0;return;end;u=unique(z);b=0;for k=reshape(u,1,[]),q=(z==k);mk=mean(Y(q,:),1);b=b+sum(q)*sum((mk-mu).^2);end;r=b/tot;end
function w=partition_wcss(Y,lab),w=0;u=unique(lab(:));for j=1:numel(u),q=lab==u(j);if any(q),Z=Y(q,:)-mean(Y(q,:),1);w=w+sum(Z(:).^2);end,end,end
function s=silhouette_mean(Y,lab),n=size(Y,1);if n<3||numel(unique(lab))<2,s=NaN;return;end;ss=sum(Y.^2,2);D=sqrt(max(0,ss+ss'-2*(Y*Y')));vals=zeros(n,1);u=unique(lab(:));for i=1:n,qi=lab==lab(i);ni=sum(qi);if ni<=1,vals(i)=0;continue;end;a=sum(D(i,qi))/(ni-1);b=Inf;for j=1:numel(u),if u(j)==lab(i),continue;end;q=lab==u(j);if any(q),b=min(b,mean(D(i,q)));end,end;vals(i)=(b-a)/max([a b eps]);end;s=mean(vals,'omitnan');end
function ch=ch_index(Y,idx,w),n=size(Y,1);k=numel(unique(idx));if k<=1||n<=k,ch=NaN;return;end;Yc=Y-mean(Y,1);t=sum(Yc(:).^2);b=max(t-w,0);ch=(b/(k-1))/max(w/(n-k),eps);end
function q=quantile_local(x,p),x=sort(x(isfinite(x)));if isempty(x),q=NaN;return;end;pos=1+(numel(x)-1)*p;lo=floor(pos);hi=ceil(pos);if lo==hi,q=x(lo);else,q=x(lo)+(pos-lo)*(x(hi)-x(lo));end,end
function a=ari_local(x,y),x=x(:);y=y(:);n=numel(x);if n<2,a=1;return;end;[~,~,ix]=unique(x);[~,~,iy]=unique(y);T=accumarray([ix iy],1);ai=sum(T,2);bj=sum(T,1);c2=@(v)sum(v.*(v-1)/2);nij=c2(T(:));aa=c2(ai);bb=c2(bj(:));den=n*(n-1)/2;ex=aa*bb/den;mx=.5*(aa+bb);if abs(mx-ex)<eps,a=double(nij==mx);else,a=(nij-ex)/(mx-ex);end,end
