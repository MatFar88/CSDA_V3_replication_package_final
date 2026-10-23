function OUT = RUN_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1(worker_id,n_workers,documents_root)
%RUN_TFKM_TRKM_FRESH_THEORY_CLOSURE_AUDITS_V1
% Three truth-open closure audits on the frozen 1230-DGP fresh study:
%   A) empirical component-respecting g-path replay,
%   B) literal midpoint-threshold cross-fitted rank replay,
%   C) boundary-directed local perturbation-and-return.
%
% No selector retuning. No new DGP cells. No 100-start refit.

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

outdir=fullfile(outroot,'fresh_theory_closure_audits_v1');
if ~isfolder(outdir),mkdir(outdir);end
grfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_grank.csv',worker_id,n_workers));
locfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_local_return.csv',worker_id,n_workers));
matfile=fullfile(outdir,sprintf('worker_%02d_of_%02d_checkpoint.mat',worker_id,n_workers));
if isfile(grfile), GR=readtable(grfile,'TextType','string'); else, GR=table(); end
if isfile(locfile), LR=readtable(locfile,'TextType','string'); else, LR=table(); end

CFG=struct();
CFG.tauW=1e-5;
CFG.tauT=1e-3;
CFG.tol=1e-12;
CFG.methods=["TFKM","TRKM"];
CFG.perturb_counts=[1 2 4];
CFG.binary_random_starts=8;
CFG.binary_max_iter=100;

unit_index=(1:height(manifest))';
assigned=mod(unit_index-1,n_workers)+1==worker_id;
Jobs=manifest(assigned,:);

fprintf('\n============================================================\n');
fprintf('TFKM/TRKM FRESH THEORY CLOSURE AUDITS V1 | worker %d/%d\n',worker_id,n_workers);
fprintf('TRUTH OPEN | EXISTING 1230 DGPs | NO SELECTOR RETUNING\n');
fprintf('Audits: empirical g-path, literal rank replay, local perturbation/return\n');
fprintf('Jobs: %d datasets\n',height(Jobs));
fprintf('============================================================\n\n');

for iu=1:height(Jobs)
    sid=double(Jobs.scenario_id(iu)); rep=double(Jobs.replicate(iu));
    ic=find(double(catalog.scenario_id)==sid,1); assert(~isempty(ic),'scenario_id missing.');
    c=catalog(ic,:); cfg=table2struct(c);
    alpha0=double(c.alpha0);
    expected_local_per_unit=20; % 2 methods * [baseline + 3 counts * 3 perturbation types]
    grDone=~isempty(GR) && any(double(GR.scenario_id)==sid & double(GR.replicate)==rep);
    if isempty(LR)
        nloc=0;
    else
        nloc=sum(double(LR.scenario_id)==sid & double(LR.replicate)==rep);
    end
    if grDone && nloc==expected_local_per_unit
        fprintf('sid=%d rep=%02d already complete\n',sid,rep);
        continue
    end

    data_seed=double(Jobs.data_seed(iu));
    sim=generate_fixed_dataset_v1(cfg,data_seed);
    [y,out,U0]=truth_fields_local(sim,double(c.r0));
    X=double(sim.X); U0=orth_local(U0); in=~out;
    z=relabel_local(y(in));
    g0=double(c.g0); r0=double(c.r0); n=size(X,1); p=size(X,2);
    assert(g0<=5,'This closure audit is pre-specified for the frozen study with g0<=5.');
    assert(numel(unique(z))==g0 && size(U0,2)==r0,'Truth geometry mismatch.');
    assert(sum(out)==floor(alpha0*n),'Truth mask size does not match alpha0.');
    zi=full_labels_local(in,z,n);

    [Cpop,centroid_source]=population_centroids_local(sim,X(in,:),z,U0,g0,r0);
    nk=accumarray(z,1,[g0 1]); piw=nk/sum(nk);

    % ------------------------------------------------------------------
    % A) Empirical component-respecting g path.
    % ------------------------------------------------------------------
    if ~grDone
        G=gpath_audit_local(X(in,:)*U0,z,g0,data_seed,CFG);
        Rk=rank_audit_local(X,in,z,U0,Cpop,piw,r0,data_seed);
        row=table(sid,rep,data_seed,string(c.panel),string(c.contamination),n,p,double(c.p_over_n), ...
            double(c.rho),double(c.delta),r0,g0,alpha0,centroid_source, ...
            G.R1,G.R2,G.R3,G.R4,G.R5,G.R6,G.min_pre_gain,G.first_overfit_gain,G.Gg_emp, ...
            G.change_point_separated,G.binary_split_gain_max,G.binary_split_group, ...
            Rk.beta_min,Rk.cross_error,Rk.G_r_cross,Rk.midpoint_threshold,Rk.rhat_midpoint, ...
            Rk.rank_exact_midpoint,Rk.weyl_sufficient_pass,Rk.lambda_r0,Rk.lambda_r0p1, ...
            'VariableNames',{'scenario_id','replicate','data_seed','panel','contamination','n','p','p_over_n', ...
            'rho','delta','r0','g0','alpha0','centroid_source', ...
            'Rhat_g1','Rhat_g2','Rhat_g3','Rhat_g4','Rhat_g5','Rhat_g6','min_pre_gain','first_overfit_gain','Gg_emp', ...
            'g_change_point_separated','binary_split_gain_max','binary_split_group', ...
            'beta_min','cross_error','G_r_cross','midpoint_threshold','rhat_midpoint', ...
            'rank_exact_midpoint','weyl_sufficient_pass','lambda_r0','lambda_r0p1'});
        if isempty(GR),GR=row;else,GR=[GR;row];end %#ok<AGROW>
        GR=sortrows(GR,{'scenario_id','replicate'});
        writetable(GR,grfile);
    end

    % ------------------------------------------------------------------
    % C) Local perturbation-and-return, method specific.
    % ------------------------------------------------------------------
    % If this unit is partially present from an interrupted run, remove the
    % partial rows and deterministically rebuild the complete 20-row block.
    if nloc>0 && nloc~=expected_local_per_unit
        keep=~(double(LR.scenario_id)==sid & double(LR.replicate)==rep);
        LR=LR(keep,:);
    elseif nloc==expected_local_per_unit
        save(matfile,'GR','LR','Jobs','CFG','worker_id','n_workers','-v7.3');
        continue
    end

    for method=CFG.methods
        Base=one_step_local(X,in,zi,g0,r0,sum(out),method,CFG);
        baseARI=adjusted_rand_index_v1(z,Base.pred(in));
        baseAssign=abs(baseARI-1)<=1e-12;
        baseTrim=isequal(~Base.in_next,logical(out));
        baseFull=baseAssign && baseTrim;

        LR=append_local_row_local(LR,c,sid,rep,data_seed,method,"BASELINE",0,true, ...
            0,0,0,0,baseARI,baseAssign,baseTrim,baseFull,Base.sample_gap);

        for kk=CFG.perturb_counts
            % LABEL perturbation.
            [inP,labP,meta]=make_perturbation_local(X,in,zi,out,Base,kk,"LABEL",g0);
            if meta.applicable
                U=one_step_local(X,inP,labP,g0,r0,sum(out),method,CFG);
                ari=adjusted_rand_index_v1(z,U.pred(in));
                ae=abs(ari-1)<=1e-12; te=isequal(~U.in_next,logical(out)); fe=ae&&te;
            else
                ari=NaN;ae=false;te=false;fe=false;U=struct('sample_gap',NaN);
            end
            LR=append_local_row_local(LR,c,sid,rep,data_seed,method,"LABEL",kk,meta.applicable, ...
                meta.label_flips,meta.trim_swaps,meta.label_disagreement_rate,meta.trim_disagreement_rate, ...
                ari,ae,te,fe,U.sample_gap);

            % TRIM perturbation.
            [inP,labP,meta]=make_perturbation_local(X,in,zi,out,Base,kk,"TRIM",g0);
            if meta.applicable
                U=one_step_local(X,inP,labP,g0,r0,sum(out),method,CFG);
                ari=adjusted_rand_index_v1(z,U.pred(in));
                ae=abs(ari-1)<=1e-12; te=isequal(~U.in_next,logical(out)); fe=ae&&te;
            else
                ari=NaN;ae=false;te=false;fe=false;U=struct('sample_gap',NaN);
            end
            LR=append_local_row_local(LR,c,sid,rep,data_seed,method,"TRIM",kk,meta.applicable, ...
                meta.label_flips,meta.trim_swaps,meta.label_disagreement_rate,meta.trim_disagreement_rate, ...
                ari,ae,te,fe,U.sample_gap);

            % JOINT perturbation.
            [inP,labP,meta]=make_perturbation_local(X,in,zi,out,Base,kk,"JOINT",g0);
            if meta.applicable
                U=one_step_local(X,inP,labP,g0,r0,sum(out),method,CFG);
                ari=adjusted_rand_index_v1(z,U.pred(in));
                ae=abs(ari-1)<=1e-12; te=isequal(~U.in_next,logical(out)); fe=ae&&te;
            else
                ari=NaN;ae=false;te=false;fe=false;U=struct('sample_gap',NaN);
            end
            LR=append_local_row_local(LR,c,sid,rep,data_seed,method,"JOINT",kk,meta.applicable, ...
                meta.label_flips,meta.trim_swaps,meta.label_disagreement_rate,meta.trim_disagreement_rate, ...
                ari,ae,te,fe,U.sample_gap);
        end
    end

    LR=sortrows(LR,{'scenario_id','replicate','method','perturb_type','perturb_count'});
    writetable(LR,locfile);
    save(matfile,'GR','LR','Jobs','CFG','worker_id','n_workers','-v7.3');

    gu=GR(double(GR.scenario_id)==sid & double(GR.replicate)==rep,:);
    fprintf('sid=%d rep=%02d | Gg_emp=%7.3f gsep=%d | rhat=%d/%d rankOK=%d | local rows=%d\n', ...
        sid,rep,double(gu.Gg_emp(1)),logical(gu.g_change_point_separated(1)), ...
        double(gu.rhat_midpoint(1)),r0,logical(gu.rank_exact_midpoint(1)),expected_local_per_unit);
end

fprintf('\nWorker %d/%d complete: %d structural rows, %d local rows.\n',worker_id,n_workers,height(GR),height(LR));
fprintf('%s\n%s\n',grfile,locfile);
OUT=struct('grank',GR,'local_return',LR,'output_root',string(outdir),'worker_id',worker_id,'n_workers',n_workers);
end

% =========================================================================
function T=append_local_row_local(T,c,sid,rep,data_seed,method,ptype,k,applicable,label_flips,trim_swaps,label_rate,trim_rate,ari,assign_exact,trim_exact,full_exact,sample_gap)
state_distance=label_rate+trim_rate;
row=table(sid,rep,data_seed,method,string(c.panel),string(c.contamination),double(c.p_over_n),double(c.rho),double(c.delta), ...
    double(c.r0),double(c.g0),double(c.alpha0),string(ptype),double(k),logical(applicable),double(label_flips),double(trim_swaps), ...
    double(label_rate),double(trim_rate),double(state_distance),double(ari),logical(assign_exact),logical(trim_exact),logical(full_exact),double(sample_gap), ...
    'VariableNames',{'scenario_id','replicate','data_seed','method','panel','contamination','p_over_n','rho','delta', ...
    'r0','g0','alpha0','perturb_type','perturb_count','applicable','label_flips','trim_swaps', ...
    'label_disagreement_rate','trim_disagreement_rate','designed_state_distance','return_inlier_ARI','return_assignment_exact', ...
    'return_trim_exact','return_full_state_exact','sample_generalized_gap_after_perturb'});
if isempty(T),T=row;else,T=[T;row];end %#ok<AGROW>
end

% =========================================================================
function G=gpath_audit_local(Y,z,g0,data_seed,CFG)
Y=double(Y);z=double(z(:));
mu=zeros(g0,size(Y,2));nk=zeros(g0,1);W0=0;
for k=1:g0
    Ik=(z==k);nk(k)=sum(Ik);assert(nk(k)>=2,'Too few observations in true group %d.',k);
    mu(k,:)=mean(Y(Ik,:),1);E=Y(Ik,:)-mu(k,:);W0=W0+sum(E(:).^2);
end
R=nan(1,6);
for g=1:g0
    if g==g0
        R(g)=W0;
    else
        R(g)=W0+best_component_merge_extra_local(mu,nk,g);
    end
end
bestGain=-Inf;bestGroup=NaN;
for k=1:g0
    Yk=Y(z==k,:);
    gain=binary_split_gain_local(Yk,data_seed+7919*k,CFG);
    if gain>bestGain
        bestGain=gain;bestGroup=k;
    end
end
R(g0+1)=max(W0-bestGain,0);
D=nan(1,6);
for g=2:min(g0+1,6),D(g)=R(g-1)-R(g);end
pre=D(2:g0);pre=pre(isfinite(pre));assert(~isempty(pre),'No pre-g0 gains.');
minPre=min(pre);post=D(g0+1);
Gemp=minPre/max(post,eps);
G=struct('R1',R(1),'R2',R(2),'R3',R(3),'R4',R(4),'R5',R(5),'R6',R(6), ...
    'min_pre_gain',minPre,'first_overfit_gain',post,'Gg_emp',Gemp, ...
    'change_point_separated',minPre>post+1e-12,'binary_split_gain_max',bestGain,'binary_split_group',bestGroup);
end

function best=best_component_merge_extra_local(mu,nk,g)
g0=size(mu,1);assert(g>=1 && g<g0);
% Enumerate all surjective maps from the g0 true components to g blocks.
% Duplicate labelings are harmless because only the minimum is retained.
N=g^g0;best=Inf;
for code=0:N-1
    q=code;map=zeros(g0,1);
    for j=1:g0
        map(j)=mod(q,g)+1;q=floor(q/g);
    end
    if numel(unique(map))<g,continue,end
    extra=0;
    for b=1:g
        I=find(map==b);nb=sum(nk(I));
        mb=sum(mu(I,:).*nk(I),1)/nb;
        E=mu(I,:)-mb;
        extra=extra+sum(nk(I).*sum(E.^2,2));
    end
    if extra<best,best=extra;end
end
assert(isfinite(best),'Failed component-respecting partition enumeration.');
end

function gain=binary_split_gain_local(Y,seed,CFG)
Y=double(Y);n=size(Y,1);assert(n>=4,'Binary split requires at least four observations.');
mu=mean(Y,1);E=Y-mu;sse1=sum(E(:).^2);
starts=cell(0,1);
% PC-sign start.
[~,~,V]=svd(E,'econ');v=V(:,1);q=E*v;I1=q<=0;I2=~I1;
if any(I1)&&any(I2),starts{end+1}=[mean(Y(I1,:),1);mean(Y(I2,:),1)];end %#ok<AGROW>
% Farthest pair start.
D2=sqdist_local(Y,Y);D2(1:n+1:end)=-Inf;[~,ix]=max(D2(:));[i,j]=ind2sub([n n],ix);starts{end+1}=[Y(i,:);Y(j,:)]; %#ok<AGROW>
% Deterministic random-pair starts.
s=RandStream('mt19937ar','Seed',mod(floor(abs(seed)),2^31-2)+1);
for b=1:CFG.binary_random_starts
    ij=randperm(s,n,2);starts{end+1}=[Y(ij(1),:);Y(ij(2),:)]; %#ok<AGROW>
end
best=sse1;
for b=1:numel(starts)
    val=lloyd2_sse_local(Y,starts{b},CFG.binary_max_iter);
    if val<best,best=val;end
end
gain=max(sse1-best,0);
end

function sse=lloyd2_sse_local(Y,C,maxit)
n=size(Y,1);lab=zeros(n,1);
for it=1:maxit
    D=sqdist_local(Y,C);[~,newlab]=min(D,[],2);
    if any(accumarray(newlab,1,[2 1])==0)
        % Deterministic repair: move the point farthest from its assigned center.
        mind=min(D,[],2);[~,ii]=max(mind);newlab(ii)=3-newlab(ii);
    end
    Cnew=C;
    for k=1:2,Cnew(k,:)=mean(Y(newlab==k,:),1);end
    if isequal(newlab,lab) || max(abs(Cnew(:)-C(:)))<1e-12,C=Cnew;lab=newlab;break,end
    C=Cnew;lab=newlab;
end
D=sqdist_local(Y,C);sse=sum(min(D,[],2));
end

function D=sqdist_local(A,B)
D=sum(A.^2,2)+sum(B.^2,2)'-2*(A*B');D=max(D,0);
end

% =========================================================================
function R=rank_audit_local(X,in,z,U0,Cpop,piw,r0,data_seed)
p=size(X,2);g0=max(z);
[B0r,beta_min]=between_target_local(Cpop,piw,r0);B0=U0*B0r*U0';B0=(B0+B0')/2;
idx=find(in);lab=zeros(size(in));lab(in)=z;
muA=zeros(g0,p);muB=zeros(g0,p);
for k=1:g0
    Ik=idx(z==k);nk=numel(Ik);assert(nk>=4,'Need >=4 inliers per group for stratified split.');
    s=RandStream('mt19937ar','Seed',mod(floor(abs(data_seed+104729*k)),2^31-2)+1);
    ord=Ik(randperm(s,nk));
    na=floor(nk/2);A=ord(1:na);B=ord(na+1:end);
    muA(k,:)=mean(X(A,:),1);muB(k,:)=mean(X(B,:),1);
end
mA=piw'*muA;mB=piw'*muB;
Bx=zeros(p);
for k=1:g0
    a=muA(k,:)-mA;b=muB(k,:)-mB;
    Bx=Bx+.5*piw(k)*(a'*b+b'*a);
end
Bx=(Bx+Bx')/2;E=(Bx-B0);E=(E+E')/2;
epsilon=max(abs(real(eig(E))));
lam=sort(real(eig(Bx)),'descend');
t=beta_min/2;rhat=sum(lam>t+1e-12);
if r0<=numel(lam),lr=lam(r0);else,lr=NaN;end
if r0+1<=numel(lam),lr1=lam(r0+1);else,lr1=NaN;end
R=struct('beta_min',beta_min,'cross_error',epsilon,'G_r_cross',beta_min/max(2*epsilon,eps), ...
    'midpoint_threshold',t,'rhat_midpoint',rhat,'rank_exact_midpoint',rhat==r0, ...
    'weyl_sufficient_pass',epsilon<beta_min/2,'lambda_r0',lr,'lambda_r0p1',lr1);
end

% =========================================================================
function O=one_step_local(X,inState,labState,g0,r0,mtrim,method,CFG)
assert(numel(inState)==size(X,1) && numel(labState)==size(X,1),'State length mismatch.');
assert(all(labState(inState)>=1 & labState(inState)<=g0),'Invalid retained labels.');
assert(all(accumarray(labState(inState),1,[g0 1])>0),'Empty retained group.');
S=state_scatter_from_full_labels_local(X,inState,labState,g0,CFG);
if method=="TFKM"
    A=S.SW+S.lamW*eye(S.p);B=S.ST+S.lamT*eye(S.p);direction="ascend";
else
    A=S.SB;B=S.ST+S.lamT*eye(S.p);direction="descend";
end
[U,~,gap]=generalized_space_local(A,B,r0,direction);
Y=X*U;C=zeros(g0,r0);
for k=1:g0,C(k,:)=mean(Y(inState & labState==k,:),1);end
DD=sqdist_local(Y,C);[dist,pred]=min(DD,[],2);
trim=false(size(inState));if mtrim>0,[~,ord]=sort(dist,'descend');trim(ord(1:mtrim))=true;end
O=struct('U',U,'C',C,'DD',DD,'dist',dist,'pred',pred,'in_next',~trim,'sample_gap',gap);
end

function [inP,labP,M]=make_perturbation_local(X,in,zi,out,Base,k,ptype,g0) %#ok<INUSD>
n=numel(in);inP=in;labP=zi;M=struct('applicable',true,'label_flips',0,'trim_swaps',0,'label_disagreement_rate',0,'trim_disagreement_rate',0);
ptype=string(ptype);
if (ptype=="TRIM" || ptype=="JOINT") && (~any(out) || k>sum(out) || k>sum(in))
    M.applicable=false;return
end

swappedOutRetained=[];
if ptype=="TRIM" || ptype=="JOINT"
    ii=find(in);oo=find(out);
    [~,a]=sort(Base.dist(ii),'descend');takeIn=ii(a(1:k));
    [~,b]=sort(Base.dist(oo),'ascend');takeOut=oo(b(1:k));
    inP(takeIn)=false;labP(takeIn)=0;
    inP(takeOut)=true;labP(takeOut)=Base.pred(takeOut);
    swappedOutRetained=takeIn;
    M.trim_swaps=k;
    M.trim_disagreement_rate=mean(inP~=in);
end

if ptype=="LABEL" || ptype=="JOINT"
    cand=find(in & inP);
    if ~isempty(swappedOutRetained),cand=setdiff(cand,swappedOutRetained,'stable');end
    trueLab=zi(cand);DD=Base.DD(cand,:);
    DDcomp=DD;
    for q=1:numel(cand),DDcomp(q,trueLab(q))=Inf;end
    [secondD,comp]=min(DDcomp,[],2);
    firstD=Base.DD(sub2ind(size(Base.DD),cand,trueLab));
    margin=secondD-firstD;
    [~,ord]=sort(margin,'ascend');
    chosen=[];
    for jj=ord(:)'
        idx=cand(jj);old=labP(idx);new=comp(jj);
        if old==new,continue,end
        % Preserve nonempty retained groups.
        if sum(inP & labP==old)<=1,continue,end
        chosen(end+1)=jj; %#ok<AGROW>
        labP(idx)=new;
        if numel(chosen)==k,break,end
    end
    if numel(chosen)<k,M.applicable=false;return,end
    M.label_flips=k;
    M.label_disagreement_rate=k/n;
end

assert(sum(inP)==sum(in),'Perturbation must preserve trimming mass.');
assert(all(accumarray(labP(inP),1,[g0 1])>0),'Perturbation emptied a retained group.');
end

% =========================================================================
function S=state_scatter_from_full_labels_local(X,in,lab,g,CFG)
n=size(X,1);p=size(X,2);h=sum(in);mu=mean(X(in,:),1);Xc=X-mu;
D=diag(double(in));G=zeros(n,g);idx=find(in);for ii=1:numel(idx),G(idx(ii),lab(idx(ii)))=1;end
nk=sum(G,1);assert(all(nk>0),'Empty retained group.');P=G*diag(1./nk)*G';
ST=Xc'*D*Xc;SW=Xc'*(D-P)*Xc;SB=ST-SW;ST=(ST+ST')/2;SW=(SW+SW')/2;SB=(SB+SB')/2;
rW=max(rank(full(SW)),1);rT=max(rank(full(ST)),1);sW=trace(SW)/rW;sT=trace(ST)/rT;
if ~isfinite(sW)||sW<=0,sW=max(sT,1);end
if ~isfinite(sT)||sT<=0,sT=max(sW,1);end
S=struct('SW',SW,'ST',ST,'SB',SB,'lamW',CFG.tauW*sW,'lamT',CFG.tauT*sT,'h',h,'p',p);
end

function [U,vals,gap]=generalized_space_local(A,B,r,direction)
A=(A+A')/2;B=(B+B')/2;[V,D]=eig(A,B,'vector');D=real(D(:));ok=isfinite(D)&all(isfinite(real(V)),1)';idx=find(ok);assert(numel(idx)>=r+1,'Too few finite generalized eigenvectors.');
if direction=="ascend",[~,ord]=sort(D(idx),'ascend');else,[~,ord]=sort(D(idx),'descend');end
sel=idx(ord(1:r));U=real(V(:,sel));[U,~]=qr(U,0);vals=D(sel);allv=D(idx(ord));
if direction=="ascend",gap=allv(r+1)-allv(r);else,gap=allv(r)-allv(r+1);end
gap=max(real(gap),0);
end

% =========================================================================
function [B,beta_min]=between_target_local(Cpop,piw,r0)
Cc=Cpop-(piw'*Cpop);B=zeros(r0);for k=1:size(Cc,1),B=B+piw(k)*(Cc(k,:)'*Cc(k,:));end;B=(B+B')/2;
b=sort(max(real(eig(B)),0),'descend');assert(numel(b)>=r0);beta_min=b(r0);
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

function lab=full_labels_local(in,z,n)
lab=zeros(n,1);lab(in)=z;
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
% g-path helper: three well-separated empirical groups must show a clear
% pre/post gain separation.
rng(71);Y=[randn(20,2)*.10+[-4 0];randn(20,2)*.10+[0 0];randn(20,2)*.10+[4 0]];z=[ones(20,1);2*ones(20,1);3*ones(20,1)];
CFG=struct('binary_random_starts',4,'binary_max_iter',50);
G=gpath_audit_local(Y,z,3,991,CFG);assert(G.Gg_emp>1 && G.change_point_separated,'g-path selftest failed.');
% literal midpoint rank counting.
B=diag([4 2 0 0]);beta=2;t=beta/2;assert(sum(eig(B)>t)==2,'Rank midpoint selftest failed.');
% local update: simple three-cluster/no-trim state should be stable for both
% geometries at r=2.
X=[randn(20,3)*.05+[-4 0 0];randn(20,3)*.05+[0 4 0];randn(20,3)*.05+[4 0 0]];
in=true(60,1);lab=[ones(20,1);2*ones(20,1);3*ones(20,1)];CFG.tauW=1e-5;CFG.tauT=1e-3;
for method=["TFKM","TRKM"]
    O=one_step_local(X,in,lab,3,2,0,method,CFG);assert(all(O.in_next),'Local no-trim selftest failed.');
end
fprintf('SELFTEST PASS: empirical g-path, literal rank count, and local-return machinery.\n');
end
