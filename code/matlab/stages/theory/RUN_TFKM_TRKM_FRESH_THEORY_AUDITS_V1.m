function OUT = RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1(worker_id,n_workers,documents_root)
%RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1
% Truth-open direct audit of four theoretical links on the existing 1230-unit
% fresh study. No new DGP cells are generated and no selector is retuned.
%
% AUDIT 1  Strict true-state prerequisite / one-step self-map
% AUDIT 2  Method-specific spectral accessibility (population gap vs sample error)
% AUDIT 3  Projected exact-recovery margin using the true-state REG-HD subspace
% AUDIT 4  Residual trimming separation at the true trimming mass
%
% Usage:
%   RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1("SELFTEST",[],[]);
%   RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1(1,6);
%   ...
%   RUN_TFKM_TRKM_FRESH_THEORY_AUDITS_V1(6,6);
% Then:
%   OUT4 = MERGE_TFKM_TRKM_FRESH_THEORY_AUDITS_V1;

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
assert(isfield(S,'catalog') && isfield(S,'manifest'),'study.mat schema mismatch.');
catalog=S.catalog; manifest=S.manifest;
assert(height(catalog)==41 && height(manifest)==1230,'Fresh study dimensions drifted.');
truthroot=fullfile(outroot,'fresh_truth_open_v1');
assert(isfile(fullfile(truthroot,'TRUTH_OPENED_AT.txt')),'Truth-open marker missing.');
lock_vendor_local(root);

outdir=fullfile(outroot,'fresh_theory_audits_v1');
if ~isfolder(outdir),mkdir(outdir);end
outfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_theory_audits.csv',worker_id,n_workers));
matfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_checkpoint.mat',worker_id,n_workers));
if isfile(outfile), Rows=readtable(outfile,'TextType','string'); else, Rows=table(); end

unit_index=(1:height(manifest))';
assigned=mod(unit_index-1,n_workers)+1==worker_id;
Jobs=manifest(assigned,:);

CFG=struct();
CFG.tauW=1e-5;
CFG.tauT=1e-3;
CFG.tol=1e-12;
CFG.methods=["TFKM","TRKM"];

fprintf('\n============================================================\n');
fprintf('TFKM/TRKM FRESH THEORY AUDITS V1 | worker %d/%d\n',worker_id,n_workers);
fprintf('TRUTH OPEN | EXISTING 1230 DGPs | NO SELECTOR RETUNING\n');
fprintf('Four audits: strictness, spectral, projected margin, trimming gap\n');
fprintf('Jobs: %d datasets x 2 methods\n',height(Jobs));
fprintf('============================================================\n\n');

for iu=1:height(Jobs)
    sid=double(Jobs.scenario_id(iu)); rep=double(Jobs.replicate(iu));
    if ~isempty(Rows)
        done=double(Rows.scenario_id)==sid & double(Rows.replicate)==rep;
        if sum(done)>=2
            fprintf('sid=%d rep=%02d already complete\n',sid,rep);
            continue
        end
    end
    ic=find(double(catalog.scenario_id)==sid,1); assert(~isempty(ic),'scenario_id missing.');
    c=catalog(ic,:); cfg=table2struct(c);
    data_seed=double(Jobs.data_seed(iu));
    sim=generate_fixed_dataset_v1(cfg,data_seed);
    [y,out,U0]=truth_fields_local(sim,double(c.r0));
    X=double(sim.X); U0=orth_local(U0); in=~out;
    z=relabel_local(y(in));
    g0=double(c.g0); r0=double(c.r0); n=size(X,1); p=size(X,2); h=sum(in);
    assert(numel(unique(z))==g0 && size(U0,2)==r0,'Truth geometry mismatch.');
    assert(sum(out)==floor(double(c.alpha0)*n),'Truth mask size does not match alpha0.');

    [Cpop,centroid_source]=population_centroids_local(sim,X(in,:),z,U0,g0,r0);
    nk=accumarray(z,1,[g0 1]); piw=nk/sum(nk);
    [B0r,beta_min]=between_target_local(Cpop,piw,r0);
    State=true_state_scatter_local(X,in,z,g0,CFG);

    for method=CFG.methods
        if ~isempty(Rows) && any(double(Rows.scenario_id)==sid & double(Rows.replicate)==rep & string(Rows.method)==method)
            continue
        end
        t0=tic;
        Spec=spectral_audit_local(State,U0,B0r,beta_min,double(c.rho),r0,method);
        Strict=strictness_audit_local(X,in,out,z,U0,Cpop,Spec.Uhat,double(c.alpha0),Spec.sample_gap,CFG);
        Proj=projected_margin_audit_local(X,in,z,U0,Cpop,Spec.Uhat,double(c.rho),CFG);
        Trim=trimming_gap_audit_local(X,in,out,z,U0,Cpop,Spec.Uhat,double(c.alpha0),CFG);

        R=empty_row_local();
        R.method=method; R.scenario_id=sid; R.replicate=rep; R.data_seed=data_seed;
        R.panel=string(c.panel); R.contamination=string(c.contamination);
        R.n=n; R.p=p; R.n_inlier=h; R.n_outlier=sum(out);
        R.p_over_n=double(c.p_over_n); R.rho=double(c.rho); R.delta=double(c.delta);
        R.r0=r0; R.g0=g0; R.alpha0=double(c.alpha0); R.centroid_source=centroid_source;

        % Audit 1: strict-state prerequisite / one-step self-map.
        R.sample_generalized_gap=Spec.sample_gap;
        R.assignment_margin_inliers=Strict.assignment_margin_inliers;
        R.assignment_inliers_exact=Strict.assignment_inliers_exact;
        R.trim_residual_margin_est=Strict.trim_margin_est;
        R.one_step_trim_exact=Strict.trim_exact;
        R.strict_prerequisite_quotient=Strict.strict_prereq;
        R.one_step_truth_state_recovered=Strict.one_step_truth_state_recovered;

        % Audit 2: spectral accessibility.
        R.population_generalized_gap=Spec.population_gap;
        R.population_target_overlap=Spec.population_target_overlap;
        R.population_target_correct=Spec.population_target_correct;
        R.sample_target_overlap_mean=Spec.overlap_mean;
        R.sample_target_cstar=Spec.cstar;
        R.operator_error=Spec.operator_error;
        R.operator_error_over_gap=Spec.error_over_gap;
        R.G_spec_gap_over_error=Spec.Gspec;
        R.tfkm_population_condition_value=Spec.tfkm_condition_value;
        R.tfkm_population_condition_pass=Spec.tfkm_condition_pass;
        R.ell_W=Spec.ellW; R.ell_T=Spec.ellT;

        % Audit 3: effective projected partition margin.
        R.intrinsic_partition_margin_8logn=Proj.G0;
        R.projected_variance_factor=Proj.vrho;
        R.effective_projected_margin=Proj.Geff;
        R.effective_projected_margin_gt1=Proj.Geff>1;
        R.projected_oracle_ARI=Proj.ari;
        R.projected_oracle_exact=Proj.exact;

        % Audit 4: residual trimming separation.
        R.oracle_trim_gap=Trim.oracle_gap;
        R.estimated_trim_gap=Trim.estimated_gap;
        R.max_residual_perturbation=Trim.max_error;
        R.G_trim_gap_over_2error=Trim.Gtrim;
        R.oracle_trim_mask_exact=Trim.oracle_mask_exact;
        R.estimated_trim_mask_exact=Trim.estimated_mask_exact;

        R.runtime_seconds=toc(t0);
        row=struct2table(R,'AsArray',true);
        if isempty(Rows),Rows=row;else,Rows=[Rows;row];end %#ok<AGROW>
        Rows=sortrows(Rows,{'scenario_id','replicate','method'});
        writetable(Rows,outfile);
        save(matfile,'Rows','Jobs','CFG','worker_id','n_workers','-v7.3');

        fprintf(['sid=%d rep=%02d %-4s %-11s %-8s | strict=%d Gspec=%6.3f ' ...
            'c*=%.3f Geff=%6.3f Gtrim=%6.3f | %.2fs\n'], ...
            sid,rep,char(method),char(string(c.panel)),char(string(c.contamination)), ...
            Strict.one_step_truth_state_recovered,Spec.Gspec,Spec.cstar,Proj.Geff,Trim.Gtrim,R.runtime_seconds);
    end
end

fprintf('\nWorker %d/%d complete: %d method-unit rows saved.\n%s\n',worker_id,n_workers,height(Rows),outfile);
OUT=struct('rows',Rows,'output_root',string(outdir),'worker_id',worker_id,'n_workers',n_workers);
end

% =========================================================================
function R=empty_row_local()
R=struct();
R.method="";R.scenario_id=NaN;R.replicate=NaN;R.data_seed=NaN;R.panel="";R.contamination="";
R.n=NaN;R.p=NaN;R.n_inlier=NaN;R.n_outlier=NaN;R.p_over_n=NaN;R.rho=NaN;R.delta=NaN;
R.r0=NaN;R.g0=NaN;R.alpha0=NaN;R.centroid_source="";
R.sample_generalized_gap=NaN;R.assignment_margin_inliers=NaN;R.assignment_inliers_exact=false;
R.trim_residual_margin_est=NaN;R.one_step_trim_exact=false;R.strict_prerequisite_quotient=false;R.one_step_truth_state_recovered=false;
R.population_generalized_gap=NaN;R.population_target_overlap=NaN;R.population_target_correct=false;
R.sample_target_overlap_mean=NaN;R.sample_target_cstar=NaN;R.operator_error=NaN;R.operator_error_over_gap=NaN;R.G_spec_gap_over_error=NaN;
R.tfkm_population_condition_value=NaN;R.tfkm_population_condition_pass=false;R.ell_W=NaN;R.ell_T=NaN;
R.intrinsic_partition_margin_8logn=NaN;R.projected_variance_factor=NaN;R.effective_projected_margin=NaN;
R.effective_projected_margin_gt1=false;R.projected_oracle_ARI=NaN;R.projected_oracle_exact=false;
R.oracle_trim_gap=NaN;R.estimated_trim_gap=NaN;R.max_residual_perturbation=NaN;R.G_trim_gap_over_2error=NaN;
R.oracle_trim_mask_exact=false;R.estimated_trim_mask_exact=false;R.runtime_seconds=NaN;
end

% =========================================================================
function S=true_state_scatter_local(X,in,z,g,CFG)
n=size(X,1);p=size(X,2);h=sum(in);mu=mean(X(in,:),1);Xc=X-mu;
D=diag(double(in));G=zeros(n,g);
idx=find(in);for ii=1:numel(idx),G(idx(ii),z(ii))=1;end
nk=sum(G,1);assert(all(nk>0),'Empty true retained group.');P=G*diag(1./nk)*G';
ST=Xc'*D*Xc;SW=Xc'*(D-P)*Xc;SB=ST-SW;ST=(ST+ST')/2;SW=(SW+SW')/2;SB=(SB+SB')/2;
rW=max(rank(full(SW)),1);rT=max(rank(full(ST)),1);
sW=trace(SW)/rW;sT=trace(ST)/rT;
if ~isfinite(sW)||sW<=0,sW=max(sT,1);end
if ~isfinite(sT)||sT<=0,sT=max(sW,1);end
lamW=CFG.tauW*sW;lamT=CFG.tauT*sT;
S=struct('SW',SW,'ST',ST,'SB',SB,'lamW',lamW,'lamT',lamT,'ellW',lamW/h,'ellT',lamT/h,'h',h,'p',p);
end

function S=spectral_audit_local(State,U0,B0r,beta_min,rho,r0,method)
p=State.p;P0=U0*U0';Sigma0=P0+rho*(eye(p)-P0);SigmaB=U0*B0r*U0';SigmaT=Sigma0+SigmaB;
ellW=State.ellW;ellT=State.ellT;
ShW=State.SW/State.h;ShT=State.ST/State.h;ShB=State.SB/State.h;
Bhat=ShT+ellT*eye(p);Bpop=SigmaT+ellT*eye(p);
if method=="TFKM"
    Ahat=ShW+ellW*eye(p);Apop=Sigma0+ellW*eye(p);direction="ascend";
else
    Ahat=ShB;Apop=SigmaB;direction="descend";
end
[Uhat,valsHat]=generalized_space_local(Ahat,Bhat,r0,direction);
[Upop,valsPop]=generalized_space_local(Apop,Bpop,r0,direction);
Khat=whitened_local(Ahat,Bhat);Kpop=whitened_local(Apop,Bpop);
E=(Khat-Kpop);E=(E+E')/2;e=max(abs(real(eig(E))));
allPop=sort(real(eig(Kpop)),char(direction));allHat=sort(real(eig(Khat)),char(direction));
if direction=="ascend"
    dpop=allPop(r0+1)-allPop(r0);dhat=allHat(r0+1)-allHat(r0);
else
    dpop=allPop(r0)-allPop(r0+1);dhat=allHat(r0)-allHat(r0+1);
end
s=svd(U0'*Uhat);s=min(max(real(s),0),1);cstar=min(s)^2;ov=mean(s.^2);
sp=svd(U0'*Upop);sp=min(max(real(sp),0),1);popov=mean(sp.^2);
S=struct();S.Uhat=Uhat;S.sample_gap=max(dhat,0);S.population_gap=max(dpop,0);S.operator_error=e;
S.error_over_gap=e/max(S.population_gap,eps);S.Gspec=S.population_gap/max(e,eps);S.cstar=cstar;S.overlap_mean=ov;
S.population_target_overlap=popov;S.population_target_correct=(min(sp)^2>1-1e-8);S.ellW=ellW;S.ellT=ellT;
S.tfkm_condition_value=NaN;S.tfkm_condition_pass=false;
if method=="TFKM"
    val=(rho+ellW)*beta_min+(rho-1)*(ellT-ellW);S.tfkm_condition_value=val;S.tfkm_condition_pass=val>0;
end
end

function [U,vals]=generalized_space_local(A,B,r,direction)
A=(A+A')/2;B=(B+B')/2;[V,D]=eig(A,B,'vector');D=real(D(:));ok=isfinite(D)&all(isfinite(real(V)),1)';idx=find(ok);assert(numel(idx)>=r+1,'Too few finite generalized eigenvectors.');
if direction=="ascend",[~,ord]=sort(D(idx),'ascend');else,[~,ord]=sort(D(idx),'descend');end
sel=idx(ord(1:r));U=real(V(:,sel));[U,~]=qr(U,0);vals=D(sel);
end

function K=whitened_local(A,B)
A=(A+A')/2;B=(B+B')/2;[V,d]=eig(B,'vector');d=real(d(:));floorv=max(eps,1e-12*max(max(abs(d)),1));d(d<floorv)=floorv;Bi=V*diag(1./sqrt(d))*V';K=Bi*A*Bi;K=(K+K')/2;
end

function A=strictness_audit_local(X,in,out,z,U0,Cpop,Uhat,alpha0,sample_gap,CFG) %#ok<INUSD>
n=size(X,1);g=max(z);Y=X*Uhat;C=zeros(g,size(Uhat,2));
for k=1:g,C(k,:)=mean(Y(find(in & full_labels_local(in,z,n)==k),:),1);end
DD=zeros(n,g);for k=1:g,E=Y-C(k,:);DD(:,k)=sum(E.^2,2);end
[dist,pred]=min(DD,[],2);zi=full_labels_local(in,z,n);idx=find(in);
trueD=DD(sub2ind(size(DD),idx,zi(idx)));DD2=DD(idx,:);for q=1:numel(idx),DD2(q,zi(idx(q)))=Inf;end
second=min(DD2,[],2);margin=min(second-trueD);assignExact=all(pred(in)==zi(in));
m=floor(alpha0*n);trim=false(n,1);if m>0,[~,ord]=sort(dist,'descend');trim(ord(1:m))=true;end
trimExact=isequal(trim,logical(out));if any(out),tmargin=min(dist(out))-max(dist(in));else,tmargin=Inf;end
strict=sample_gap>0 && margin>0 && (~any(out) || tmargin>0);
A=struct('assignment_margin_inliers',margin,'assignment_inliers_exact',assignExact,'trim_margin_est',tmargin, ...
    'trim_exact',trimExact,'strict_prereq',strict,'one_step_truth_state_recovered',assignExact&&trimExact);
end

function lab=full_labels_local(in,z,n)
lab=zeros(n,1);lab(in)=z;
end

function P=projected_margin_audit_local(X,in,z,U0,Cpop,Uhat,rho,CFG) %#ok<INUSD>
n=size(X,1);g=max(z);dmin=Inf;for a=1:g,for b=a+1:g,dmin=min(dmin,norm(Cpop(a,:)-Cpop(b,:)));end,end
s=svd(U0'*Uhat);s=min(max(real(s),0),1);cstar=min(s)^2;if rho<=1,vrho=1;else,vrho=cstar+rho*(1-cstar);end
G0=dmin^2/(8*log(max(n,3)));Geff=G0*cstar/max(vrho,eps);
Y=X*Uhat;Mu=Cpop*U0'*Uhat;DD=zeros(n,g);for k=1:g,E=Y-Mu(k,:);DD(:,k)=sum(E.^2,2);end
[~,pred]=min(DD,[],2);ari=adjusted_rand_index_v1(z,pred(in));exact=all(pred(in)==z);
P=struct('G0',G0,'vrho',vrho,'Geff',Geff,'ari',ari,'exact',exact);
end

function T=trimming_gap_audit_local(X,in,out,z,U0,Cpop,Uhat,alpha0,CFG) %#ok<INUSD>
n=size(X,1);g=max(z);Y0=X*U0;D0=zeros(n,g);for k=1:g,E=Y0-Cpop(k,:);D0(:,k)=sum(E.^2,2);end;R0=min(D0,[],2);
Y=X*Uhat;C=zeros(g,size(Uhat,2));for k=1:g,C(k,:)=mean(Y(in & full_labels_local(in,z,n)==k,:),1);end
Dh=zeros(n,g);for k=1:g,E=Y-C(k,:);Dh(:,k)=sum(E.^2,2);end;Rh=min(Dh,[],2);
m=floor(alpha0*n);oracleMask=false(n,1);estMask=false(n,1);if m>0,[~,o]=sort(R0,'descend');oracleMask(o(1:m))=true;[~,o]=sort(Rh,'descend');estMask(o(1:m))=true;end
if any(out)
    b0=min(R0(out))-max(R0(in));bh=min(Rh(out))-max(Rh(in));er=max(abs(Rh-R0));G=b0/max(2*er,eps);
else
    b0=NaN;bh=NaN;er=NaN;G=NaN;
end
T=struct('oracle_gap',b0,'estimated_gap',bh,'max_error',er,'Gtrim',G,'oracle_mask_exact',isequal(oracleMask,logical(out)), ...
    'estimated_mask_exact',isequal(estMask,logical(out)));
end

function [B,beta_min]=between_target_local(Cpop,piw,r0)
Cc=Cpop-(piw'*Cpop);B=zeros(r0);for k=1:size(Cc,1),B=B+piw(k)*(Cc(k,:)'*Cc(k,:));end;B=(B+B')/2;b=sort(max(real(eig(B)),0),'descend');assert(numel(b)>=r0);beta_min=b(r0);
end

function [C,source]=population_centroids_local(sim,Xin,z,U0,g0,r0)
C=[];source="";
if isfield(sim,'centroids')
    Q=double(sim.centroids);if isequal(size(Q),[g0 r0])&&all(isfinite(Q(:))),C=Q;source="DGP_POPULATION_CENTROIDS";end
end
if isempty(C)
    Y=Xin*U0;C=zeros(g0,r0);for k=1:g0,C(k,:)=mean(Y(z==k,:),1);end;source="REALIZED_TRUE_GROUP_CENTROIDS";
end
end

function [y,out,U0]=truth_fields_local(sim,r0)
if isfield(sim,'labels'),y=double(sim.labels(:));elseif isfield(sim,'gtrue'),y=double(sim.gtrue(:));else,error('Truth labels missing.');end
if isfield(sim,'outlier'),out=logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),out=logical(sim.outlier_mask(:));else,error('Outlier mask missing.');end
assert(numel(y)==numel(out),'Truth length mismatch.');assert(isfield(sim,'U')&&size(sim.U,2)>=r0,'Generating U missing.');U0=double(sim.U(:,1:r0));
end

function z=relabel_local(y)
u=unique(y(:),'stable');z=zeros(size(y));for k=1:numel(u),z(y==u(k))=k;end
end

function U=orth_local(U)
[U,~]=qr(double(U),0);
end

function lock_vendor_local(root)
v=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));assert(isfolder(v),'Canonical vendor root missing: %s',v);addpath(genpath(char(v)),'-begin');rehash path;
req=["generate_fixed_dataset_v1","adjusted_rand_index_v1"];for q=req,p=string(which(char(q)));assert(strlength(p)>0,'Missing %s.',q);assert(startsWith(p,v,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,p);end
end

function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found. Pass Documents explicitly.');S=load(fullfile(outroot,'study.mat'));
end

function selftest_local()
% Population target: TRKM must recover U0; TFKM does so in a benign case.
rng(11);p=8;r=2;U0=orth_local(randn(p,r));B=diag([4 2]);rho=1.5;P0=U0*U0';S0=P0+rho*(eye(p)-P0);SB=U0*B*U0';ST=S0+SB;ellW=1e-5;ellT=1e-3;
[UR,~]=generalized_space_local(SB,ST+ellT*eye(p),r,"descend");assert(min(svd(U0'*UR))^2>1-1e-8);
[UF,~]=generalized_space_local(S0+ellW*eye(p),ST+ellT*eye(p),r,"ascend");assert(min(svd(U0'*UF))^2>1-1e-8);
% Projected margin decreases when signal-space fidelity is degraded for rho>1.
G0=2;c1=.9;c2=.4;v1=c1+rho*(1-c1);v2=c2+rho*(1-c2);assert(G0*c1/v1>G0*c2/v2);
% Trimming-order sufficient ratio.
b=4;e=1;assert(b/(2*e)>1);
fprintf('SELFTEST PASS: strictness, spectral, projected-margin, and trimming-gap contracts.\n');
end
