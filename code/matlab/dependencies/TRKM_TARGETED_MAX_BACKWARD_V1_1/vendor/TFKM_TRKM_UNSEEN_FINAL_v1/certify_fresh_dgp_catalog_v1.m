function out=certify_fresh_dgp_catalog_v1(output_root)
%CERTIFY_FRESH_DGP_CATALOG_V1 Independent truth-open DGP certification.
% Uses a seed namespace disjoint from the untouched validation reps.
pkg=fileparts(mfilename('fullpath')); addpath(genpath(pkg),'-begin'); output_root=char(string(output_root));
S=load(fullfile(output_root,'study.mat')); study=S.study; catalog=S.catalog;
reps=double(study.certification_reps(:)); B=300; Lmax=3;
Units=table(); RankDir=table();
fprintf('\n============================================================\n');
fprintf('FRESH DGP CERTIFICATION | independent seed namespace\n');
fprintf('NO TFKM/TRKM fits. Validation seeds remain untouched.\n');
fprintf('============================================================\n\n');
for ii=1:height(catalog)
    cfg=table2struct(catalog(ii,:)); sid=double(cfg.scenario_id);
    fprintf('sid=%d source=%d | %s | true=(%d,%d,%.3f) | delta=%.4g rho=%.4g\n',sid,cfg.source_confirm_id,string(cfg.contamination),cfg.r0,cfg.c0,cfg.alpha0,cfg.delta,cfg.rho);
    for rep=reshape(reps,1,[])
        seed=double(study.certification_base_seed)+sid*100000+rep; sim=generate_fixed_dataset_v1(cfg,seed);
        cm=cluster_metrics(sim,cfg,seed); [rm,rd]=rank_metrics(sim,cfg,sid,rep,seed,B,Lmax);
        row=table(sid,rep,seed,double(cfg.r0),double(cfg.c0),double(cfg.alpha0),double(cfg.delta),double(cfg.rho),string(cfg.contamination), ...
            cm.oracle_ari,cm.g0_ari_u,cm.g0_ari_aug,cm.best_g_u_sil,cm.best_g_u_ch,cm.best_g_aug_sil,cm.best_g_aug_ch, ...
            cm.sil_margin_u,cm.ch_margin_rel_u,cm.sil_margin_aug,cm.ch_margin_rel_aug,logical(cm.strict), ...
            rm.centroid_rank,logical(rm.population_rank_cert),rm.centroid_smin_ratio,rm.weakest_true_R2,rm.noise_pc_R2_max,rm.random_R2_q95, ...
            rm.weakest_margin_vs_noise_pc,rm.weakest_margin_vs_random_q95,logical(rm.rank_strict_sample_cert),rm.kmin_oracle_ari99, ...
            'VariableNames',{'scenario_id','cert_replicate','seed','r0','g0','alpha0','delta','rho','contamination', ...
            'oracle_ari','kmeans_g0_ari_u0','kmeans_g0_ari_aug','best_g_u0_sil','best_g_u0_ch','best_g_aug_sil','best_g_aug_ch', ...
            'sil_margin_u0','ch_margin_rel_u0','sil_margin_aug','ch_margin_rel_aug','cluster_strict_cert', ...
            'centroid_rank','population_rank_cert','centroid_smin_ratio','weakest_true_R2','noise_pc_R2_max','random_R2_q95', ...
            'weakest_margin_vs_noise_pc','weakest_margin_vs_random_q95','rank_strict_sample_cert','kmin_oracle_ari99'});
        Units=[Units;row]; RankDir=[RankDir;rd]; %#ok<AGROW>
        fprintf('  cert %d | cluster=%d oracle=%.3f g_aug=%d/%d margins=(%.3g,%.3g) | rank=%d weakR2=%.3f noise=%.3f rnd95=%.3f\n', ...
            rep,cm.strict,cm.oracle_ari,cm.best_g_aug_sil,cm.best_g_aug_ch,cm.sil_margin_aug,cm.ch_margin_rel_aug,rm.rank_strict_sample_cert,rm.weakest_true_R2,rm.noise_pc_R2_max,rm.random_R2_q95);
    end
    fprintf('\n');
end
Summary=table();
for sid=unique(double(Units.scenario_id))'
    Z=Units(double(Units.scenario_id)==sid,:); row=table(sid,height(Z),mean(double(Z.cluster_strict_cert)),mean(double(Z.rank_strict_sample_cert)), ...
        min(double(Z.oracle_ari)),min(double(Z.sil_margin_aug)),min(double(Z.ch_margin_rel_aug)),min(double(Z.weakest_margin_vs_random_q95)), ...
        'VariableNames',{'scenario_id','N','cluster_cert_rate','rank_cert_rate','min_oracle_ari','min_sil_margin_aug','min_ch_margin_rel_aug','min_rank_margin_vs_random_q95'});
    Summary=[Summary;row]; %#ok<AGROW>
end
all_pass=all(Summary.cluster_cert_rate==1 & Summary.rank_cert_rate==1);
outdir=fullfile(output_root,'dgp_certification'); writetable(Units,fullfile(outdir,'dgp_certification_units.csv')); writetable(RankDir,fullfile(outdir,'rank_signal_directions.csv')); writetable(Summary,fullfile(outdir,'dgp_certification_summary.csv'));
save(fullfile(outdir,'DGP_CERTIFICATION_PASS.mat'),'all_pass','Summary','Units','-v7.3');
fprintf('\nDGP CERTIFICATION SUMMARY\n'); disp(Summary);
if all_pass
    fprintf('DGP CERTIFICATION: PASS. Untouched validation sessions are UNLOCKED.\n');
else
    fprintf('DGP CERTIFICATION: FAIL. Do NOT launch validation sessions.\n');
end
out=struct('units',Units,'rank_directions',RankDir,'summary',Summary,'all_pass',all_pass);
end

function cm=cluster_metrics(sim,cfg,seed)
in=truth_inlier_mask(sim); X=double(sim.X(in,:)); z=double(sim.labels(in)); U0=true_U(sim); Xc=X-mean(X,1); Y=Xc*U0; C=true_centroids(sim,cfg,z,Y); p=accumarray(z,1,[double(cfg.c0) 1]); p=p/sum(p); C=C-(p'*C);
cm.oracle_ari=ari_local(z,nearest_center(Y,C)); Gu=kgrid(Y,z,seed+710000); cm.best_g_u_sil=best_g(Gu,"sil"); cm.best_g_u_ch=best_g(Gu,"ch"); cm.g0_ari_u=value_g(Gu,double(cfg.c0),"ari"); [cm.sil_margin_u,cm.ch_margin_rel_u]=margins(Gu,double(cfg.c0));
R=Xc-(Xc*U0)*U0'; if norm(R,'fro')>0,[~,~,V]=svd(R,'econ'); zz=R*V(:,1);else,zz=zeros(size(X,1),1);end; zz=zz-mean(zz);
Ga=kgrid([Y zz],z,seed+720000); cm.best_g_aug_sil=best_g(Ga,"sil"); cm.best_g_aug_ch=best_g(Ga,"ch"); cm.g0_ari_aug=value_g(Ga,double(cfg.c0),"ari"); [cm.sil_margin_aug,cm.ch_margin_rel_aug]=margins(Ga,double(cfg.c0));
cm.strict=cm.oracle_ari>=.99 && cm.g0_ari_u>=.99 && cm.g0_ari_aug>=.99 && cm.best_g_u_sil==double(cfg.c0) && cm.best_g_u_ch==double(cfg.c0) && cm.best_g_aug_sil==double(cfg.c0) && cm.best_g_aug_ch==double(cfg.c0) && cm.sil_margin_u>1e-6 && cm.ch_margin_rel_u>1e-6 && cm.sil_margin_aug>1e-6 && cm.ch_margin_rel_aug>1e-6;
end
function G=kgrid(Y,z,seed)
G=table(); for g=2:min(6,size(Y,1)-1),rng(seed+g,'twister');try,[idx,~,sumd]=kmeans(Y,g,'Replicates',15,'MaxIter',300,'Start','plus','Display','off');w=sum(sumd);si=silhouette_mean(Y,idx);ch=ch_index(Y,idx,w);ar=ari_local(z,idx);catch,w=NaN;si=NaN;ch=NaN;ar=NaN;end;G=[G;table(g,ar,si,ch,w,'VariableNames',{'g','ari','sil','ch','wcss'})];end
end
function [sm,cm]=margins(G,g0)
i=double(G.g)==g0; s0=double(G.sil(i)); c0=double(G.ch(i)); so=double(G.sil(~i)); co=double(G.ch(~i)); so=so(isfinite(so));co=co(isfinite(co));if isempty(so),sm=NaN;else,sm=s0-max(so);end;if isempty(co),cm=NaN;else,cm=(c0-max(co))/max(abs(c0),1);end
end
function g=best_g(G,var),x=double(G.(char(var)));x(~isfinite(x))=-Inf;mx=max(x);cand=double(G.g(abs(x-mx)<=1e-12));g=min(cand);end
function v=value_g(G,g,var),i=find(double(G.g)==g,1);if isempty(i),v=NaN;else,v=double(G.(char(var))(i));end,end

function [rm,DIR]=rank_metrics(sim,cfg,sid,rep,seed,B,Lmax)
in=truth_inlier_mask(sim); X=double(sim.X(in,:)); z=double(sim.labels(in)); Xc=X-mean(X,1); U0=true_U(sim); r0=size(U0,2); p=size(X,2); Y=Xc*U0; g0=double(cfg.c0);
[~,pi]=group_weights(z,g0); Cpop=true_centroids(sim,cfg,z,Y); Cpop=Cpop-(pi'*Cpop); Wc=diag(sqrt(pi))*Cpop; sv=svd(Wc,'econ'); tol=max(size(Wc))*eps(max([sv;1])); crank=sum(sv>tol); if isempty(sv)||sv(1)<=0,sratio=NaN;else,sratio=sv(min(r0,numel(sv)))/sv(1);end
Bpop=Cpop'*diag(pi)*Cpop; Bpop=(Bpop+Bpop')/2; [V,Dm]=eig(Bpop); ev=real(diag(Dm));[ev,ord]=sort(ev,'descend');V=real(V(:,ord)); ev=ev(:); V=V(:,1:r0); ev=ev(1:r0);
trueR2=nan(r0,1); DIR=table(); for j=1:r0,s=Y*V(:,j);trueR2(j)=label_R2(s,z);DIR=[DIR;table(sid,rep,j,ev(j),ev(j)/max(sum(ev),eps),trueR2(j),'VariableNames',{'scenario_id','cert_replicate','signal_axis','population_between_eigenvalue','population_between_energy_fraction','true_label_R2'})];end
weak=min(trueR2); R=Xc-(Xc*U0)*U0'; npc=[]; if norm(R,'fro')>0,[~,~,Vr]=svd(R,'econ');L=min([Lmax,size(Vr,2)]);npc=nan(L,1);for j=1:L,npc(j)=label_R2(Xc*Vr(:,j),z);end,end;if isempty(npc),noiseMax=0;else,noiseMax=max(npc);end
rng(seed+910000,'twister'); Vb=randn(p,B);Vb=Vb-U0*(U0'*Vb);nn=sqrt(sum(Vb.^2,1));good=nn>sqrt(eps);Vb(:,good)=Vb(:,good)./nn(good);Vb(:,~good)=0;Sc=Xc*Vb;rv=nan(B,1);for b=1:B,rv(b)=label_R2(Sc(:,b),z);end;rq95=quantile_local(rv,.95);
aris=nan(r0,1); for k=1:r0,Vk=V(:,1:k);aris(k)=ari_local(z,nearest_center(Y*Vk,Cpop*Vk));end;kmin=find(aris>=.99,1,'first');if isempty(kmin),kmin=NaN;end
rm=struct('centroid_rank',crank,'population_rank_cert',crank==r0,'centroid_smin_ratio',sratio,'weakest_true_R2',weak,'noise_pc_R2_max',noiseMax,'random_R2_q95',rq95,'weakest_margin_vs_noise_pc',weak-noiseMax,'weakest_margin_vs_random_q95',weak-rq95,'rank_strict_sample_cert',(crank==r0)&&weak>max(noiseMax,rq95)+1e-12,'kmin_oracle_ari99',kmin);
end
function [cnt,p]=group_weights(z,g),cnt=accumarray(z,1,[g 1]);p=cnt/sum(cnt);end
function in=truth_inlier_mask(sim),n=size(sim.X,1);if isfield(sim,'outlier'),in=~logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),in=~logical(sim.outlier_mask(:));else,in=true(n,1);end,end
function U=true_U(sim),if isfield(sim,'U'),U=double(sim.U);elseif isfield(sim,'U_true'),U=double(sim.U_true);else,error('Missing true U.');end;[U,~]=qr(U,0);end
function C=true_centroids(sim,cfg,z,Y),if isfield(sim,'centroids')&&size(sim.centroids,1)==double(cfg.c0)&&size(sim.centroids,2)==size(Y,2),C=double(sim.centroids);else,C=zeros(double(cfg.c0),size(Y,2));for k=1:double(cfg.c0),C(k,:)=mean(Y(z==k,:),1);end,end,end
function lab=nearest_center(Y,C),D=zeros(size(Y,1),size(C,1));for k=1:size(C,1),Z=Y-C(k,:);D(:,k)=sum(Z.^2,2);end;[~,lab]=min(D,[],2);end
function r=label_R2(Y,z),if isvector(Y),Y=Y(:);end;mu=mean(Y,1);tot=sum((Y-mu).^2,'all');if tot<=eps,r=0;return;end;u=unique(z);b=0;for k=reshape(u,1,[]),q=(z==k);mk=mean(Y(q,:),1);b=b+sum(q)*sum((mk-mu).^2);end;r=b/tot;end
function s=silhouette_mean(Y,lab),n=size(Y,1);if n<3||numel(unique(lab))<2,s=NaN;return;end;ss=sum(Y.^2,2);D=sqrt(max(0,ss+ss'-2*(Y*Y')));vals=zeros(n,1);u=unique(lab(:));for i=1:n,qi=(lab==lab(i));ni=sum(qi);if ni<=1,vals(i)=0;continue;end;a=sum(D(i,qi))/(ni-1);b=Inf;for j=1:numel(u),if u(j)==lab(i),continue;end;q=(lab==u(j));if any(q),b=min(b,mean(D(i,q)));end,end;vals(i)=(b-a)/max([a b eps]);end;s=mean(vals,'omitnan');end
function ch=ch_index(Y,idx,w),n=size(Y,1);k=numel(unique(idx));if k<=1||n<=k,ch=NaN;return;end;Yc=Y-mean(Y,1);t=sum(Yc(:).^2);b=max(t-w,0);ch=(b/(k-1))/max(w/(n-k),eps);end
function q=quantile_local(x,p),x=sort(x(isfinite(x)));if isempty(x),q=NaN;return;end;pos=1+(numel(x)-1)*p;lo=floor(pos);hi=ceil(pos);if lo==hi,q=x(lo);else,q=x(lo)+(pos-lo)*(x(hi)-x(lo));end,end
function a=ari_local(x,y),x=x(:);y=y(:);n=numel(x);if n<2,a=1;return;end;[~,~,ix]=unique(x);[~,~,iy]=unique(y);T=accumarray([ix iy],1);ai=sum(T,2);bj=sum(T,1);c2=@(v)sum(v.*(v-1)/2);nij=c2(T(:));aa=c2(ai);bb=c2(bj(:));den=n*(n-1)/2;ex=aa*bb/den;mx=.5*(aa+bb);if abs(mx-ex)<eps,a=double(nij==mx);else,a=(nij-ex)/(mx-ex);end,end
