function OUT = RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1(scenario_id,replicate,documents_root)
%RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1
% Truth-blind TFKM finalization from the CORRECTED fresh HIER4 cubes.
%
% Recovery provenance:
%   * The original fresh TFKM cubes used SIGN_ALL_STABLE and are excluded.
%   * Inputs here are ONLY the correction-rerun cubes under
%       TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output/fits/TFKM
%   * The frozen TFKM procedure is unchanged:
%       whole-path HIER4 winners -> modal set -> componentwise-maximal M+
%       -> frozen action/alpha rule
%       -> PROPAGATED_SUPPORT_NATIVE -> mutual-retained ARI >= .90.
%   * No truth labels, outlier mask, r0/g0/alpha0, contamination name, rho,
%     delta, or posthoc recovery metric is exposed to selection or output.
%   * Regenerated DGP truth fields are discarded immediately after X is read.
%
% Usage:
%   RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1(51001,1);
%   RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1("SOURCECHECK",[]);
%   RUN_TFKM_FRESH_HIER4_BLIND_UNIT_V1("SELFTEST",[]);

if nargin>=1 && strcmpi(string(scenario_id),"SELFTEST")
    selftest_local();
    OUT=struct('selftest_pass',true);
    return
end
if nargin>=1 && strcmpi(string(scenario_id),"SOURCECHECK")
    if nargin<3,documents_root="";end
    [root,~,~]=load_study_local(documents_root);
    sourcecheck_tfkm_local(root);
    OUT=struct('sourcecheck_pass',true);
    return
end
if nargin<2,error('Provide scenario_id and replicate.');end
if nargin<3,documents_root="";end
method="REG_TFKM_HD";
selftest_local();

[root,outroot,S]=load_study_local(documents_root);
CFG0=S.CFG;options=S.options;catalog=S.catalog;manifest=S.manifest;
assert(replicate>=1 && replicate<=CFG0.n_replicates && replicate==round(replicate),'replicate out of range.');
ic=find(double(catalog.scenario_id)==double(scenario_id),1);assert(~isempty(ic),'Unknown scenario_id.');
im=find(double(manifest.scenario_id)==double(scenario_id)&double(manifest.replicate)==double(replicate),1);assert(~isempty(im),'Manifest row missing.');
data_seed=double(manifest.data_seed(im));
pair_seed=double(manifest.tfkm_pair_seed(im));

% Regenerate X in an isolated helper. No truth/generative fields leave it.
X=regenerate_X_only_local(catalog(ic,:),data_seed,root,method);
[n,p]=size(X);

% IMPORTANT: outroot is the HIER4 correction-rerun root, never the legacy
% fresh TFKM output root.
fitdir=fullfile(outroot,'fits','TFKM',sprintf('s%d',scenario_id),sprintf('rep%02d',replicate));
matfile=fullfile(fitdir,'cube.mat');completefit=fullfile(fitdir,'COMPLETE.txt');
assert(isfile(matfile)&&isfile(completefit),'Corrected HIER4 cube incomplete/missing: %s',fitdir);
Q=load(matfile,'Candidates','Labels','Trim','unit_info');
C=Q.Candidates;Lstored=Q.Labels;Tstored=Q.Trim;
assert(height(C)==150&&all(logical(C.valid)),'Corrected candidate cube incomplete/invalid.');
assert(size(Lstored,1)==n&&size(Lstored,2)==150&&isequal(size(Tstored),[n 150]),'Stored label/trim cube shape mismatch.');
assert(all(double(C.scenario_id)==double(scenario_id))&&all(double(C.replicate)==double(replicate)),'Cube key drift.');
assert(all(double(C.n)==n)&&all(double(C.p)==p),'Cube dimension drift.');
assert(max(abs(sort(unique(double(C.alpha)))'-double(CFG0.alpha_grid)))<=1e-12,'Alpha-grid drift.');

% Per-unit HIER4 provenance gate. This makes it impossible for the legacy
% SIGN_ALL_STABLE cubes to enter this finalizer accidentally.
reqcols=["hier_cross_ari_score","hier_cross_j_score","hier_cross_H_score","selector_score"];
assert(all(ismember(reqcols,string(C.Properties.VariableNames))),'Corrected cube missing HIER4 columns.');
A=double(C.hier_cross_ari_score);J=double(C.hier_cross_j_score);H=double(C.hier_cross_H_score);Sc=double(C.selector_score);
assert(all(isfinite(A))&&all(isfinite(J))&&all(isfinite(H))&&all(isfinite(Sc)),'Corrected cube has nonfinite HIER4 evidence.');
assert(max(abs(H-(A+.2*J)))<=1e-10,'Corrected cube HIER4 A+0.2J identity failure.');
assert(max(abs(Sc-H))<=1e-10,'Corrected cube selector_score != HIER4 score.');
if isstruct(Q.unit_info) && isfield(Q.unit_info,'hier_wrapper_path')
    hp=string(Q.unit_info.hier_wrapper_path);
    assert(contains(hp,[filesep 'vendor' filesep 'hier4' filesep],'IgnoreCase',true),'Corrected cube wrapper provenance is not vendor/hier4.');
end

canon=lock_vendor_local(root,method);
req=["make_robust_initialization_v1","regularized_trimmed_fit_corrected_v1","adjusted_rand_index_v1"];
for q=req
    f=string(which(char(q)));assert(strlength(f)>0,'Missing %s.',q);assert(startsWith(f,canon,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,f);
end

CFG=blind_cfg_local(CFG0);
finalroot=fullfile(outroot,'blind_finalization_hier4_v1');
unitdir=fullfile(finalroot,'TFKM',sprintf('s%d',scenario_id),sprintf('rep%02d',replicate));
if ~isfolder(unitdir),mkdir(unitdir);end
outfile=fullfile(unitdir,'blind_final.mat');predfile=fullfile(unitdir,'prediction.csv');edgefile=fullfile(unitdir,'edges.csv');repfile=fullfile(unitdir,'reproduction.csv');donefile=fullfile(unitdir,'COMPLETE.txt');
if isfile(donefile)&&isfile(outfile)&&isfile(predfile)
    Z=load(outfile,'Prediction');
    fprintf('TFKM-HIER4 BLIND s%d rep%02d already COMPLETE.\n',scenario_id,replicate);
    OUT=struct('complete',true,'prediction',Z.Prediction,'unit_dir',string(unitdir));return
end

[Prediction,Edges,Repro,FinalLabels,FinalTrim,Audit]=finalize_tfkm_local(X,C,Lstored,Tstored,options,pair_seed,CFG,scenario_id,replicate);

assert_no_truth_columns_local(Prediction);
assert_no_truth_columns_local(Edges);
assert_no_truth_columns_local(Repro);
if ~isempty(Repro) && ismember('reproduction_pass',Repro.Properties.VariableNames)
    assert(all(logical(Repro.reproduction_pass)),'A TFKM reproduction check failed.');
end
save(outfile,'Prediction','Edges','Repro','FinalLabels','FinalTrim','Audit','-v7.3');
writetable(Prediction,predfile);
if ~isempty(Edges),writetable(Edges,edgefile);else,writetable(empty_edge_export_local(),edgefile);end
if ~isempty(Repro),writetable(Repro,repfile);else,writetable(empty_repro_export_local(),repfile);end
fid=fopen(donefile,'w');assert(fid>0);fprintf(fid,'COMPLETE\nscenario_id=%d\nreplicate=%d\nmethod=REG_TFKM_HD\nsource=HIER4_CORRECTION_RERUN\ntruth_opened=NO\n',scenario_id,replicate);fclose(fid);
fprintf('TFKM-HIER4 BLIND s%d rep%02d COMPLETE.\n',scenario_id,replicate);
OUT=struct('complete',true,'prediction',Prediction,'unit_dir',string(unitdir));
end

% =========================================================================
% TFKM FINALIZATION
% =========================================================================
function [Pred,Edges,Repro,FinalLabels,FinalTrim,Audit]=finalize_tfkm_local(X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep)
[W,modal,Mplus]=tfkm_structure_local(C,CFG);
assert(~isempty(Mplus),'TFKM M+ empty.');
A=table();
for k=1:size(Mplus,1)
    r=Mplus(k,1);g=Mplus(k,2);P=sortrows(C(double(C.r)==r&double(C.g)==g,:),'alpha');assert(height(P)==numel(CFG.alpha_grid),'Incomplete TFKM anchor path.');
    f=tfkm_features_local(P,CFG);
    emit=f.J_qualified && ((isfinite(f.Smed)&&f.Smed<CFG.tfkm_Smed_cutoff)||f.DS_signal);
    a=NaN;if emit,a=f.alpha_proposed;end
    A=[A;table(k,r,g,f.S0,f.Smed,f.DS,f.Jmax,f.alpha_proposed,emit,a, ...
        'VariableNames',{'anchor_index','r','g','S0','Smed','DS','Jmax','alpha_proposed','emit','alpha_hat'})]; %#ok<AGROW>
end

allEmit=all(logical(A.emit));allNo=all(~logical(A.emit));alphaAgree=false;aHat=NaN;
if allEmit
    aa=double(A.alpha_hat);alphaAgree=max(aa)-min(aa)<=CFG.tol;if alphaAgree,aHat=aa(1);end
end
if allNo
    actionStatus="NO_ACTION_CONCORDANT";emitted=false;actionAmb=false;
elseif allEmit && alphaAgree
    actionStatus="EMIT_CONCORDANT";emitted=true;actionAmb=false;
else
    actionStatus="ACTION_AMBIGUOUS";emitted=false;actionAmb=true;
end

Edges=table();Repro=table();FinalLabels=[];FinalTrim=[];fitAvailable=false;fitAmbiguous=false;nEdges=0;nAccepts=0;
if actionAmb
    terminalSet=Mplus;terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    terminalStatus="ACTION_AMBIGUOUS_NO_PATH";
elseif ~emitted
    terminalSet=Mplus;terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    terminalStatus="NO_ACTION_STRUCTURE_RETAINED";
else
    Cache=empty_tfkm_cache_local();terms=zeros(0,2);fits=cell(size(Mplus,1),1);trims=cell(size(Mplus,1),1);
    for k=1:size(Mplus,1)
        [term,E,R,lab,tr,Cache,ne,na]=tfkm_path_local(X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,Mplus(k,:),aHat,k,Cache);
        terms(end+1,:)=term;Edges=[Edges;E];Repro=[Repro;R];fits{k}=lab;trims{k}=tr;nEdges=nEdges+ne;nAccepts=nAccepts+na; %#ok<AGROW>
    end
    terminalSet=unique(terms,'rows');terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;
    if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    [fitAvailable,fitAmbiguous,FinalLabels,FinalTrim]=resolve_terminal_fit_local(fits,trims,CFG);
    if ~terminalUnique,terminalStatus="AMBIGUOUS_TERMINAL";elseif fitAmbiguous,terminalStatus="UNIQUE_STRUCTURE_FIT_AMBIGUOUS";else,terminalStatus="UNIQUE_TERMINAL";end
end

Pred=table("TFKM",double(sid),double(rep),string(actionStatus),logical(emitted),double(aHat),logical(actionAmb), ...
    rg_text_local(W(:,2:3)),rg_text_local(modal),rg_text_local(Mplus),size(Mplus,1),anchor_action_text_local(A), ...
    logical(terminalUnique),rg_text_local(terminalSet),double(terminalR),double(terminalG),string(terminalStatus), ...
    logical(fitAvailable),logical(fitAmbiguous),double(nEdges),double(nAccepts), ...
    'VariableNames',{'method','scenario_id','replicate','action_status','alpha_emitted','alpha_hat','action_ambiguous', ...
    'per_alpha_winner_path','modal_set','anchor_set','n_anchors','anchor_action_detail', ...
    'terminal_unique','terminal_set','terminal_r','terminal_g','terminal_status','final_fit_available','final_fit_ambiguous','N_edges','N_accepts'});
Audit=struct('per_alpha',W,'modal',modal,'Mplus',Mplus,'anchor_action',A,'frozen_multiple_anchor_rule',true);
end

function [term,Edges,Repro,FinalLabels,FinalTrim,Cache,nEdges,nAccepts]=tfkm_path_local(X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,anchor,a,anchorIndex,Cache)
cur=anchor;propIdx=NaN;Edges=table();Repro=table();nEdges=0;nAccepts=0;FinalLabels=[];FinalTrim=[];
while has_pred_local(cur)
    nEdges=nEdges+1;pred=predecessor_local(cur);
    [PF,HF,Cache,rrF]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,cur(1),cur(2),a,"FROM");
    [PT,HT,Cache,rrT]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,pred(1),pred(2),a,"TO");
    Repro=[Repro;rrF;rrT]; %#ok<AGROW>
    CM=cross_matrix_local(PF,PT);coh=isfinite(CM)&CM>=CFG.ari_cutoff-CFG.tol;
    EF=tfkm_rep_evidence_local(PF,HF,coh,true);ET=tfkm_rep_evidence_local(PT,HT,coh,false);
    fromProp=false;
    if isfinite(propIdx)
        ii=double(propIdx);assert(PF.valid(ii)&&ismember(ii,double(EF.reps)),'Invalid propagated TFKM FROM.');fromProp=true;
    else
        ii=tfkm_choose_rep_local(PF,EF);
    end
    jj=tfkm_choose_rep_local(PT,ET);
    edge=CM(ii,jj);accept=isfinite(edge)&&edge>=CFG.ari_cutoff-CFG.tol;
    action="STOP";if accept,action="ACCEPT";end
    Edges=[Edges;table("TFKM",double(sid),double(rep),double(anchorIndex),double(a),double(nEdges),cur(1),cur(2),pred(1),pred(2), ...
        ii,jj,string(PF.family(ii)),string(PT.family(jj)),double(PF.start_id(ii)),double(PT.start_id(jj)),fromProp, ...
        double(EF.support_count(ii)),double(ET.support_count(jj)),double(EF.support_pct(ii)),double(ET.support_pct(jj)), ...
        double(EF.native_pct(ii)),double(ET.native_pct(jj)),double(edge),string(action), ...
        'VariableNames',{'method','scenario_id','replicate','anchor_index','alpha_ref','edge_index','r_from','g_from','r_to','g_to', ...
        'from_index','to_index','from_family','to_family','from_start_id','to_start_id','from_was_propagated', ...
        'from_support_count','to_support_count','from_support_pct','to_support_pct','from_native_pct','to_native_pct','selected_edge_ARI','selected_action'})]; %#ok<AGROW>
    FinalLabels=double(PF.labels(:,ii));FinalTrim=logical(PF.trim(:,ii));
    if ~accept,break;end
    nAccepts=nAccepts+1;cur=pred;propIdx=jj;FinalLabels=double(PT.labels(:,jj));FinalTrim=logical(PT.trim(:,jj));
end
if isempty(FinalLabels)
    [P,H,Cache,rr]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,cur(1),cur(2),a,"TERMINAL_BOUNDARY");Repro=[Repro;rr]; %#ok<AGROW>
    ii=double(H.chosen_index);FinalLabels=double(P.labels(:,ii));FinalTrim=logical(P.trim(:,ii));
end
term=cur;
end

function [P,H,Cache,RR]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,r,g,a,state)
key=string(sprintf('%d|%d|%d|%d|%.6f',sid,rep,r,g,a));k=find(Cache.keys==key,1);
if isempty(k)
    P=generate_tfkm_pool_local(X,r,g,a,options,pair_seed,CFG);H=hier4_decompose_local(P,CFG);
    old=candidate_row_local(C,r,g,a,CFG);reproduce_tfkm_local(P,H,old,Lstored,Tstored,CFG,sid,rep,r,g,a);
    Cache.keys(end+1,1)=key;Cache.pool{end+1,1}=P;Cache.hier{end+1,1}=H;k=numel(Cache.keys);
else
    P=Cache.pool{k};H=Cache.hier{k};
end
old=candidate_row_local(C,r,g,a,CFG);
RR=table("TFKM",double(sid),double(rep),string(state),r,g,a,string(P.family(H.chosen_index)),double(P.start_id(H.chosen_index)), ...
    double(H.cross_H),double(old.selector_score),double(H.objective),double(old.objective),true, ...
    'VariableNames',{'method','scenario_id','replicate','state','r','g','alpha','selected_family','selected_start_id','new_score','stored_score','new_objective','stored_objective','reproduction_pass'});
end

function P=generate_tfkm_pool_local(X,r,g,a,options,pair_seed,CFG)
N=sum(CFG.tfkm_alloc);n=size(X,1);valid=false(N,1);labels=zeros(n,N,'uint16');trim=false(n,N);objective=nan(N,1);family=strings(N,1);start_id=nan(N,1);z=0;
for j=1:numel(CFG.init_names)
    for h=1:CFG.tfkm_alloc(j)
        z=z+1;seed=pair_seed+1000*j+h;family(z)=CFG.init_names(j);start_id(z)=h;
        try
            rng(seed,'twister');init=make_robust_initialization_v1(X,g,r,a,CFG.init_names(j),seed);
            q=regularized_trimmed_fit_corrected_v1(X,init,r,g,a,"REG_TFKM",options.max_iterations,options.tauW,options.tauT);
            ok=all(isfinite(q.labels))&&isfinite(q.objective)&&size(q.U,2)==r;valid(z)=ok;
            if ok,labels(:,z)=uint16(q.labels(:));trim(:,z)=logical(q.trim(:));objective(z)=double(q.objective);end
        catch
            valid(z)=false;
        end
    end
end
assert(any(valid),'No valid TFKM starts at (%d,%d,%.3f).',r,g,a);
P=struct('N',N,'valid',valid,'labels',labels,'trim',trim,'objective',objective,'family',family,'start_id',start_id);
end

function H=hier4_decompose_local(P,CFG)
valid=logical(P.valid(:));fam=string(P.family(:));obj=double(P.objective(:));N=numel(valid);fw=nan(numel(CFG.init_names),1);within=nan(N,1);
for f=1:numel(CFG.init_names)
    idx=find(valid&fam==CFG.init_names(f));if isempty(idx),continue,end
    M=pair_H_local(P,idx,CFG.trim_weight);s=mean(M,2,'omitnan');within(idx)=s;mx=max(s,[],'omitnan');cand=idx(abs(s-mx)<=CFG.tol);if numel(cand)>1,[~,m]=min(obj(cand));cand=cand(m);end;fw(f)=cand;
end
win=double(fw(isfinite(fw)));assert(~isempty(win),'No HIER4 family representatives.');final=nan(N,1);
if numel(win)==1
    chosen=win(1);final(chosen)=within(chosen);
else
    M=pair_H_local(P,win,CFG.trim_weight);s=mean(M,2,'omitnan');final(win)=s;mx=max(s,[],'omitnan');cand=win(abs(s-mx)<=CFG.tol);if numel(cand)>1,[~,m]=min(obj(cand));chosen=cand(m);else,chosen=cand;end
end
H=struct('chosen_index',double(chosen),'cross_H',double(final(chosen)),'within_H',double(within(chosen)),'objective',double(obj(chosen)), ...
    'family_winner_idx',fw,'score_vector',final,'within_vector',within);
end

function E=tfkm_rep_evidence_local(P,H,coh,isFrom)
reps=double(H.family_winner_idx(isfinite(H.family_winner_idx)));reps=reps(:);N=P.N;native=nan(N,1);fit=nan(N,1);within=nan(N,1);supportPct=nan(N,1);
v=double(H.score_vector(:));native(reps)=quality_high_local(v(reps));fit(reps)=quality_low_local(double(P.objective(reps)));w=double(H.within_vector(:));within(reps)=quality_high_local(w(reps));
if isFrom,support=sum(coh,2);else,support=sum(coh,1)';end
supportPct(reps)=quality_high_local(double(support(reps)));
E=struct('reps',reps,'support_count',support,'support_pct',supportPct,'native_pct',native,'fit_pct',fit,'within_pct',within);
end
function idx=tfkm_choose_rep_local(P,E)
reps=double(E.reps(:));s=double(E.support_pct(reps));n=double(E.native_pct(reps));f=double(E.fit_pct(reps));w=double(E.within_pct(reps));obj=double(P.objective(reps));
key=[-min(s,n),-(s+n)/2,-f,-w,obj,reps];[~,o]=sortrows(key,1:size(key,2));idx=reps(o(1));
end

function reproduce_tfkm_local(P,H,old,Lstored,Tstored,CFG,sid,rep,r,g,a)
assert(string(P.family(H.chosen_index))==string(old.selected_init_family),'TFKM family reproduction failure.');
assert(abs(double(H.cross_H)-double(old.selector_score))<=CFG.tfkm_score_tol,'TFKM score reproduction failure.');
assert(abs(double(H.objective)-double(old.objective))<=CFG.obj_rel_tol*max(1,abs(double(old.objective))),'TFKM objective reproduction failure.');
if isfinite(double(old.hier_within_family_score)),assert(abs(double(H.within_H)-double(old.hier_within_family_score))<=CFG.tfkm_score_tol,'TFKM within score reproduction failure.');end
idx=double(old.row_index);tr=logical(P.trim(:,H.chosen_index));assert(isequal(tr,logical(Tstored(:,idx))),'TFKM trim reproduction failure.');
keep=~tr;ar=adjusted_rand_index_v1(double(P.labels(keep,H.chosen_index)),double(Lstored(keep,idx)));assert(isfinite(ar)&&abs(ar-1)<=1e-12,'TFKM labels reproduction failure sid=%d rep=%d (%d,%d,%.3f).',sid,rep,r,g,a);
end

function [W,modal,Mplus]=tfkm_structure_local(C,CFG)
W=zeros(numel(CFG.alpha_grid),4);
for i=1:numel(CFG.alpha_grid)
    a=CFG.alpha_grid(i);Q=C(abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(Q)==15,'TFKM alpha slice not 15.');h=double(Q.selector_score);mx=max(h);ids=find(abs(h-mx)<=CFG.tol);
    key=[double(Q.r(ids))+double(Q.g(ids)),double(Q.g(ids)),double(Q.r(ids))];[~,o]=sortrows(key,[1 2 3]);k=ids(o(1));W(i,:)=[a,double(Q.r(k)),double(Q.g(k)),h(k)];
end
[u,~,ic]=unique(W(:,2:3),'rows');cnt=accumarray(ic,1);modal=sortrows(u(cnt==max(cnt),:),[1 2]);Mplus=componentwise_maximal_local(modal);
end
function f=tfkm_features_local(P,CFG)
a=double(P.alpha);s=double(P.hier_cross_ari_score);j=double(P.hier_cross_j_score);f.S0=s(1);f.Smed=median(s);f.DS=NaN;f.DS_defined=isfinite(s(1))&&s(1)>0&&all(isfinite(s));if f.DS_defined,f.DS=(max(s)-min(s))/s(1);end
pos=find(a>0);[f.Jmax,k]=max(j(pos));f.alpha_proposed=a(pos(k));f.J_qualified=isfinite(f.Jmax)&&f.Jmax>=CFG.tfkm_J_cutoff;f.DS_signal=f.DS_defined&&f.DS>=CFG.tfkm_DS_cutoff;
end

% =========================================================================
% TRKM FINALIZATION
% =========================================================================
function [Pred,Edges,Repro,FinalLabels,FinalTrim,Audit]=finalize_trkm_local(X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep)
[archW,archModal]=trkm_archived_structure_local(C,CFG);
Pmodal=sortrows(C(double(C.r)==archModal(1)&double(C.g)==archModal(2),:),'alpha');assert(height(Pmodal)==numel(CFG.alpha_grid),'Incomplete archived modal path.');rv=arrayfun(@(k)stored_rv_local(Pmodal(k,:)),(1:height(Pmodal))');aa=double(Pmodal.alpha);rvmax=max(rv);deltaRV=rvmax-min(rv);ii=find(rv>=rvmax-CFG.trkm_onset_tol-CFG.tol,1,'first');assert(~isempty(ii),'Archived onset undefined.');A01=aa(ii);archEmit=(A01>0)&&(deltaRV>=CFG.trkm_delta_cutoff-CFG.tol);archA=NaN;if archEmit,archA=A01;end

q=radial_q_local(X);[radEmit,radA,wHi,dBIC,sep]=radial_gmm_bic_local(q,CFG);
if radEmit
    emitted=true;aHat=radA;actionSource="RADIAL_GMM_OVERRIDE";
elseif archEmit
    emitted=true;aHat=archA;actionSource="ARCHIVED_FALLBACK";
else
    emitted=false;aHat=NaN;actionSource="NO_ACTION";
end

Edges=table();Repro=table();FinalLabels=[];FinalTrim=[];fitAvailable=false;fitAmbiguous=false;nEdges=0;nAccepts=0;rawW=zeros(0,4);anchors=zeros(0,2);
if ~emitted
    terminalSet=archModal;terminalUnique=true;terminalR=archModal(1);terminalG=archModal(2);terminalStatus="NO_ACTION_ARCHIVED_MODAL";
else
    rawW=trkm_raw_winners_local(C,CFG);
    anchors=componentwise_maximal_local(rawW(rawW(:,1)>=aHat-CFG.tol,2:3));assert(~isempty(anchors),'Emitted TRKM unit has empty RAW tail envelope.');
    Cache=empty_trkm_cache_local();terms=zeros(0,2);fits=cell(size(anchors,1),1);trims=cell(size(anchors,1),1);
    for k=1:size(anchors,1)
        [term,E,R,lab,tr,Cache,ne,na]=trkm_path_local(X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,anchors(k,:),aHat,k,Cache);
        terms(end+1,:)=term;Edges=[Edges;E];Repro=[Repro;R];fits{k}=lab;trims{k}=tr;nEdges=nEdges+ne;nAccepts=nAccepts+na; %#ok<AGROW>
    end
    terminalSet=unique(terms,'rows');terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    [fitAvailable,fitAmbiguous,FinalLabels,FinalTrim]=resolve_terminal_fit_local(fits,trims,CFG);
    if ~terminalUnique,terminalStatus="AMBIGUOUS_TERMINAL";elseif fitAmbiguous,terminalStatus="UNIQUE_STRUCTURE_FIT_AMBIGUOUS";else,terminalStatus="UNIQUE_TERMINAL";end
end

Pred=table("TRKM",double(sid),double(rep),string(actionSource),logical(emitted),double(aHat),false, ...
    rg_text_local(archW(:,2:3)),rg_text_local(archModal),rg_text_local(anchors),size(anchors,1), ...
    logical(archEmit),double(archA),double(deltaRV),logical(radEmit),double(radA),double(wHi),double(dBIC),double(sep), ...
    logical(terminalUnique),rg_text_local(terminalSet),double(terminalR),double(terminalG),string(terminalStatus), ...
    logical(fitAvailable),logical(fitAmbiguous),double(nEdges),double(nAccepts), ...
    'VariableNames',{'method','scenario_id','replicate','action_status','alpha_emitted','alpha_hat','action_ambiguous', ...
    'archived_per_alpha_winner_path','archived_modal_set','anchor_set','n_anchors', ...
    'archived_emit','archived_alpha_hat','archived_delta_RV','radial_emit','radial_alpha_hat','radial_high_weight','radial_delta_BIC','radial_separation_sd', ...
    'terminal_unique','terminal_set','terminal_r','terminal_g','terminal_status','final_fit_available','final_fit_ambiguous','N_edges','N_accepts'});
Audit=struct('archived_per_alpha',archW,'archived_modal',archModal,'raw_per_alpha',rawW,'radial_q',q,'frozen_multiple_anchor_rule',true);
end

function [term,Edges,Repro,FinalLabels,FinalTrim,Cache,nEdges,nAccepts]=trkm_path_local(X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,anchor,a,anchorIndex,Cache)
cur=anchor;propIdx=NaN;Edges=table();Repro=table();nEdges=0;nAccepts=0;FinalLabels=[];FinalTrim=[];
while has_pred_local(cur)
    nEdges=nEdges+1;pred=predecessor_local(cur);
    [PF,Cache,rrF]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,cur(1),cur(2),a,"FROM");
    [PT,Cache,rrT]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,pred(1),pred(2),a,"TO");
    Repro=[Repro;rrF;rrT]; %#ok<AGROW>
    CM=cross_matrix_local(PF,PT);force=NaN;fromProp=false;if isfinite(propIdx),force=propIdx;fromProp=true;end
    [ii,jj,status]=trkm_choose_support_pair_local(PF,PT,CM,force,CFG);edge=CM(ii,jj);accept=isfinite(edge)&&edge>=CFG.ari_cutoff-CFG.tol;action="STOP";if accept,action="ACCEPT";end
    Edges=[Edges;table("TRKM",double(sid),double(rep),double(anchorIndex),double(a),double(nEdges),cur(1),cur(2),pred(1),pred(2),ii,jj, ...
        string(PF.family(ii)),string(PT.family(jj)),double(PF.start_id(ii)),double(PT.start_id(jj)),fromProp,string(status),double(edge),string(action), ...
        'VariableNames',{'method','scenario_id','replicate','anchor_index','alpha_ref','edge_index','r_from','g_from','r_to','g_to', ...
        'from_index','to_index','from_family','to_family','from_start_id','to_start_id','from_was_propagated','selector_status','selected_edge_ARI','selected_action'})]; %#ok<AGROW>
    FinalLabels=double(PF.labels(:,ii));FinalTrim=logical(PF.trim(:,ii));
    if ~accept,break;end
    nAccepts=nAccepts+1;cur=pred;propIdx=jj;FinalLabels=double(PT.labels(:,jj));FinalTrim=logical(PT.trim(:,jj));
end
if isempty(FinalLabels)
    [P,Cache,rr]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,cur(1),cur(2),a,"TERMINAL_BOUNDARY");Repro=[Repro;rr]; %#ok<AGROW>
    ii=P.chosen;FinalLabels=double(P.labels(:,ii));FinalTrim=logical(P.trim(:,ii));
end
term=cur;
end

function [P,Cache,RR]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,pair_seed,CFG,sid,rep,r,g,a,state)
key=string(sprintf('%d|%d|%d|%d|%.6f',sid,rep,r,g,a));k=find(Cache.keys==key,1);
if isempty(k)
    old=candidate_row_local(C,r,g,a,CFG);P=generate_trkm_pool_local(X,r,g,a,options,pair_seed,CFG);reproduce_trkm_local(P,old,Lstored,Tstored,CFG,sid,rep,r,g,a);
    Cache.keys(end+1,1)=key;Cache.pool{end+1,1}=P;k=numel(Cache.keys);
else
    P=Cache.pool{k};
end
old=candidate_row_local(C,r,g,a,CFG);
RR=table("TRKM",double(sid),double(rep),string(state),r,g,a,string(P.family(P.chosen)),double(P.start_id(P.chosen)),double(P.cv(P.chosen)),stored_rv_local(old),double(P.objective(P.chosen)),double(old.objective),true, ...
    'VariableNames',{'method','scenario_id','replicate','state','r','g','alpha','selected_family','selected_start_id','new_score','stored_score','new_objective','stored_objective','reproduction_pass'});
end

function P=generate_trkm_pool_local(X,r,g,a,options,pair_seed,CFG)
N=sum(CFG.trkm_alloc);n=size(X,1);valid=false(N,1);labels=zeros(n,N,'uint16');trim=false(n,N);objective=nan(N,1);family=strings(N,1);start_id=nan(N,1);candidates=cell(N,1);z=0;
for j=1:numel(CFG.init_names)
    for h=1:CFG.trkm_alloc(j)
        z=z+1;seed=pair_seed+1000*j+h;family(z)=CFG.init_names(j);start_id(z)=h;
        try
            rng(seed,'twister');init=make_robust_initialization_v1(X,g,r,a,CFG.init_names(j),seed);q=regularized_trimmed_fit_corrected_v1(X,init,r,g,a,"REG_TRKM",options.max_iterations,options.tauW,options.tauT);
            q.init_family=CFG.init_names(j);q.start_id=h;q.start_seed=seed;q.valid=all(isfinite(q.labels))&&isfinite(q.objective)&&size(q.U,2)==r;candidates{z}=q;valid(z)=q.valid;
            if q.valid,labels(:,z)=uint16(q.labels(:));trim(:,z)=logical(q.trim(:));objective(z)=double(q.objective);end
        catch
            candidates{z}=struct('valid',false,'init_family',CFG.init_names(j),'start_id',h,'start_seed',seed);
        end
    end
end
idx=find(valid);assert(~isempty(idx),'No valid TRKM starts at (%d,%d,%.3f).',r,g,a);objrank=objective_rank_local(objective,valid);
pool=struct('valid',valid,'family',family,'start_id',start_id,'objective',objective,'objective_rank',objrank,'labels',labels,'trim',trim,'trim_stability_weight',options.trim_stability_weight);
pool.current_candidate_index=legacy_candidate_local(pool,idx,options.trim_stability_weight);sel=select_robust_view_v3(pool,X,g,r,"RV_CV");chosen=sel.chosen;score=selector_vector_local(sel,'score',N);cv=selector_vector_local(sel,'cv_score',N);if all(~isfinite(cv)),cv=score;end
P=struct('N',N,'valid',valid,'family',family,'start_id',start_id,'objective',objective,'objective_rank',objrank,'labels',labels,'trim',trim,'score',score,'cv',cv,'chosen',chosen,'r',r,'g',g,'alpha',a);
end

function [ii,jj,status]=trkm_choose_support_pair_local(F,T,C,forceFrom,CFG)
famF=string(F.family(F.chosen));famT=string(T.family(T.chosen));maskF=logical(F.valid(:))&string(F.family)==famF;maskT=logical(T.valid(:))&string(T.family)==famT;
coh=(maskF*maskT'>0)&isfinite(C)&C>=CFG.ari_cutoff-CFG.tol;qRVF=quality_high_mask_local(double(F.cv),logical(F.valid));qRVT=quality_high_mask_local(double(T.cv),logical(T.valid));
sF=sum(coh,2);sT=sum(coh,1)';qSF=zeros(size(sF));qST=zeros(size(sT));idxF=find(maskF&sF>0);idxT=find(maskT&sT>0);if ~isempty(idxF),qSF(idxF)=quality_high_local(double(sF(idxF)));end;if ~isempty(idxT),qST(idxT)=quality_high_local(double(sT(idxT)));end
if isfinite(forceFrom)
    ii=double(forceFrom);assert(ii>=1&&ii<=F.N&&F.valid(ii),'Invalid propagated TRKM FROM.');assert(string(F.family(ii))==famF,'Propagated TRKM FROM left production family.');cand=find(coh(ii,:));
    if isempty(cand),jj=T.chosen;status="NO_COHERENT_PARTNER_FALLBACK_CURRENT_TO";return,end
    edge=C(ii,cand)';key=[-qST(cand),-qRVT(cand),-edge,T.objective_rank(cand),cand(:)];[~,o]=sortrows(key,1:size(key,2));jj=cand(o(1));status="PROPAGATED_SUPPORT_PARTNER";return
end
[a,b]=find(coh);if isempty(a),ii=F.chosen;jj=T.chosen;status="NO_COHERENT_CURRENT_FAMILY_PAIR_FALLBACK_CURRENT";return,end
supportMin=min(qSF(a),qST(b));rvMin=min(qRVF(a),qRVT(b));rvMean=(qRVF(a)+qRVT(b))/2;edge=C(sub2ind(size(C),a,b));key=[-supportMin,-rvMin,-rvMean,-edge,a,b];[~,o]=sortrows(key,1:size(key,2));ii=a(o(1));jj=b(o(1));status="SUPPORT_RV_MIN";
end

function reproduce_trkm_local(P,old,Lstored,Tstored,CFG,sid,rep,r,g,a)
assert(string(P.family(P.chosen))==string(old.selected_init_family),'TRKM family reproduction failure.');assert(abs(double(P.cv(P.chosen))-stored_rv_local(old))<=CFG.trkm_score_tol,'TRKM RV reproduction failure.');assert(abs(double(P.objective(P.chosen))-double(old.objective))<=CFG.obj_rel_tol*max(1,abs(double(old.objective))),'TRKM objective reproduction failure.');
idx=double(old.row_index);tr=logical(P.trim(:,P.chosen));assert(isequal(tr,logical(Tstored(:,idx))),'TRKM trim reproduction failure.');keep=~tr;ar=adjusted_rand_index_v1(double(P.labels(keep,P.chosen)),double(Lstored(keep,idx)));assert(isfinite(ar)&&abs(ar-1)<=1e-12,'TRKM labels reproduction failure sid=%d rep=%d (%d,%d,%.3f).',sid,rep,r,g,a);
end

function [W,modal]=trkm_archived_structure_local(C,CFG)
W=zeros(numel(CFG.alpha_grid),4);
for i=1:numel(CFG.alpha_grid)
    a=CFG.alpha_grid(i);P=C(abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(P)==15,'TRKM archived alpha slice not 15.');rv=arrayfun(@(k)stored_rv_local(P(k,:)),(1:height(P))');mx=max(rv);ids=find(abs(rv-mx)<=CFG.tol);
    r2=double(P.R2(ids));m2=max(r2);ids=ids(abs(r2-m2)<=CFG.tol);key=[double(P.r(ids))+double(P.g(ids)),double(P.g(ids)),double(P.r(ids))];[~,o]=sortrows(key,[1 2 3]);k=ids(o(1));W(i,:)=[a,double(P.r(k)),double(P.g(k)),rv(k)];
end
[u,~,ic]=unique(W(:,2:3),'rows');cnt=accumarray(ic,1);cand=u(cnt==max(cnt),:);key=[sum(cand,2),cand(:,2),cand(:,1)];[~,o]=sortrows(key,[-1 -2 -3]);modal=cand(o(1),:);
end
function W=trkm_raw_winners_local(C,CFG)
W=zeros(numel(CFG.alpha_grid),4);
for i=1:numel(CFG.alpha_grid)
    a=CFG.alpha_grid(i);P=C(abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(P)==15,'TRKM RAW alpha slice not 15.');sel=select_rg_frozen_v1(P,"REG_TRKM_HD");q=P(double(P.r)==sel.r&double(P.g)==sel.g,:);assert(height(q)==1,'RAW winner not unique.');W(i,:)=[a,sel.r,sel.g,stored_rv_local(q)];
end
end

function rv=stored_rv_local(q)
if ismember('selector_cv_score',q.Properties.VariableNames)&&isfinite(double(q.selector_cv_score(1))),rv=double(q.selector_cv_score(1));else,rv=double(q.selector_score(1));end
end

function q=radial_q_local(X)
med=median(X,1);sc=median(abs(X-med),1);pos=sc(isfinite(sc)&sc>0);if isempty(pos),floor_sc=1;else,floor_sc=max(eps,1e-12*median(pos));end;sc(~isfinite(sc)|sc<floor_sc)=floor_sc;Z=(X-med)./sc;q=mean(Z.^2,2);
end
function [emit,aGrid,wHi,dBIC,sep]=radial_gmm_bic_local(q,CFG)
x=log(max(double(q(:)),eps));x=x(isfinite(x));n=numel(x);assert(n>=20,'Too few radial scores.');mu=mean(x);v=max(mean((x-mu).^2),1e-8);ll1=sum(lognormpdf_local(x,mu,v));bic1=-2*ll1+2*log(n);
starts=[.30 .70 .50;.20 .80 .25;.40 .90 .20];bestLL=-Inf;best=[];
for s=1:size(starts,1)
    m1=quantile_simple_local(x,starts(s,1));m2=quantile_simple_local(x,starts(s,2));vv=max(v,1e-6);v1=vv;v2=vv;p2=starts(s,3);p1=1-p2;last=-Inf;
    for it=1:250
        l1=log(max(p1,eps))+lognormpdf_local(x,m1,v1);l2=log(max(p2,eps))+lognormpdf_local(x,m2,v2);mx=max([l1 l2],[],2);den=mx+log(exp(l1-mx)+exp(l2-mx));r2=exp(l2-den);r1=1-r2;n1=sum(r1);n2=sum(r2);if n1<2||n2<2,break,end;p1=n1/n;p2=n2/n;m1=sum(r1.*x)/n1;m2=sum(r2.*x)/n2;v1=max(sum(r1.*(x-m1).^2)/n1,1e-6);v2=max(sum(r2.*(x-m2).^2)/n2,1e-6);ll=sum(den);if abs(ll-last)<=1e-10*(1+abs(ll)),break,end;last=ll;
    end
    l1=log(max(p1,eps))+lognormpdf_local(x,m1,v1);l2=log(max(p2,eps))+lognormpdf_local(x,m2,v2);mx=max([l1 l2],[],2);ll=sum(mx+log(exp(l1-mx)+exp(l2-mx)));if ll>bestLL,bestLL=ll;best=[p1 p2 m1 m2 v1 v2];end
end
assert(~isempty(best),'GMM fit failed.');p1=best(1);p2=best(2);m1=best(3);m2=best(4);v1=best(5);v2=best(6);if m1>m2,tmp=p1;p1=p2;p2=tmp;tmp=m1;m1=m2;m2=tmp;tmp=v1;v1=v2;v2=tmp;end
wHi=p2;bic2=-2*bestLL+5*log(n);dBIC=bic1-bic2;sep=(m2-m1)/sqrt(max(.5*(v1+v2),eps));emit=isfinite(dBIC)&&dBIC>0&&isfinite(wHi)&&wHi>0&&wHi<=max(CFG.pos_grid)+CFG.tol;if emit,aGrid=project_grid_local(wHi,CFG.pos_grid);else,aGrid=NaN;end
end

% =========================================================================
% COMMON HELPERS
% =========================================================================
function CFG=blind_cfg_local(CFG0)
CFG=struct();CFG.alpha_grid=double(CFG0.alpha_grid);CFG.pos_grid=CFG.alpha_grid(CFG.alpha_grid>0);CFG.tol=1e-12;CFG.ari_cutoff=.90;CFG.trim_weight=.2;CFG.tfkm_Smed_cutoff=.44;CFG.tfkm_DS_cutoff=.12;CFG.tfkm_J_cutoff=.50;CFG.trkm_onset_tol=.01;CFG.trkm_delta_cutoff=.025;CFG.init_names=["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"];CFG.tfkm_alloc=[50 35 10 5];CFG.trkm_alloc=[40 30 15 15];CFG.tfkm_score_tol=5e-8;CFG.trkm_score_tol=1e-8;CFG.obj_rel_tol=5e-8;
end
function row=candidate_row_local(C,r,g,a,CFG)
row=C(double(C.r)==r&double(C.g)==g&abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(row)==1,'Candidate row missing/duplicated (%d,%d,%.3f).',r,g,a);
end
function Cmat=cross_matrix_local(P,Q)
Cmat=nan(P.N,Q.N);for i=1:P.N,if ~P.valid(i),continue,end;for j=1:Q.N,if ~Q.valid(j),continue,end;keep=~logical(P.trim(:,i))&~logical(Q.trim(:,j));if sum(keep)>1,Cmat(i,j)=adjusted_rand_index_v1(double(P.labels(keep,i)),double(Q.labels(keep,j)));end,end,end
end
function M=pair_H_local(P,idx,w)
idx=double(idx(:));m=numel(idx);A=zeros(m);J=zeros(m);for a=1:m,ia=idx(a);for b=a+1:m,ib=idx(b);ta=logical(P.trim(:,ia));tb=logical(P.trim(:,ib));keep=~ta&~tb;ar=0;if sum(keep)>1,ar=adjusted_rand_index_v1(double(P.labels(keep,ia)),double(P.labels(keep,ib)));if ~isfinite(ar),ar=0;end,end;un=ta|tb;if any(un),jac=sum(ta&tb)/sum(un);else,jac=1;end;A(a,b)=ar;A(b,a)=ar;J(a,b)=jac;J(b,a)=jac;end,end;M=A+w*J;
end
function q=quality_high_local(x)
x=double(x(:));q=nan(size(x));ok=isfinite(x);vals=x(ok);if isempty(vals),return,end;for i=find(ok)',q(i)=(sum(vals<x(i))+.5*sum(vals==x(i)))/numel(vals);end
end
function q=quality_low_local(x),q=quality_high_local(-double(x(:)));end
function q=quality_high_mask_local(x,valid)
x=double(x(:));valid=logical(valid(:));q=nan(size(x));idx=find(valid&isfinite(x));assert(~isempty(idx),'No finite quality values.');q(idx)=quality_high_local(x(idx));
end
function r=objective_rank_local(obj,valid),r=nan(numel(obj),1);idx=find(valid);[~,o]=sort(obj(idx),'ascend');r(idx(o))=(1:numel(idx))';end
function v=selector_vector_local(sel,name,N),v=nan(N,1);if isfield(sel,name),x=double(sel.(name)(:));v(1:min(N,numel(x)))=x(1:min(N,numel(x)));end,end
function chosen=legacy_candidate_local(pool,valid_idx,w)
obj=pool.objective(valid_idx);[~,ord]=sort(obj,'ascend');ordered=valid_idx(ord);top_n=min(numel(ordered),max(3,ceil(numel(ordered)/4)));top=ordered(1:top_n);if numel(top)==1,chosen=top(1);return,end;S=zeros(numel(top));for a=1:numel(top),for b=a+1:numel(top),ta=logical(pool.trim(:,top(a)));tb=logical(pool.trim(:,top(b)));in=~ta&~tb;if sum(in)>1,ca=adjusted_rand_index_v1(double(pool.labels(in,top(a))),double(pool.labels(in,top(b))));if ~isfinite(ca),ca=0;end;else,ca=0;end;un=ta|tb;if any(un),tj=sum(ta&tb)/sum(un);else,tj=1;end;S(a,b)=ca+w*tj;S(b,a)=S(a,b);end,end;[~,m]=max(mean(S,2));chosen=top(m);
end
function M=componentwise_maximal_local(P)
if isempty(P),M=zeros(0,2);return,end;P=unique(P,'rows');keep=true(size(P,1),1);for k=1:size(P,1),keep(k)=~any(P(:,1)>=P(k,1)&P(:,2)>=P(k,2)&(P(:,1)>P(k,1)|P(:,2)>P(k,2)));end;M=sortrows(P(keep,:),[1 2]);
end
function p=predecessor_local(q)
r=q(1);g=q(2);if g>r+1,p=[r,g-1];elseif g==r+1&&r>1,p=[r-1,g-1];else,p=[];end
end
function tf=has_pred_local(q),tf=~isempty(predecessor_local(q));end
function s=rg_text_local(R)
if isempty(R),s="{}";return,end;R=unique(double(R),'rows','stable');parts=strings(size(R,1),1);for i=1:size(R,1),parts(i)=sprintf('(%d,%d)',R(i,1),R(i,2));end;s="{"+strjoin(parts,';')+"}";
end
function s=anchor_action_text_local(A)
parts=strings(height(A),1);for i=1:height(A),if logical(A.emit(i)),aa=sprintf('%.3f',double(A.alpha_hat(i)));else,aa='NA';end;parts(i)=sprintf('(%d,%d):emit=%d:a=%s:Smed=%.6g:DS=%.6g:Jmax=%.6g',double(A.r(i)),double(A.g(i)),logical(A.emit(i)),aa,double(A.Smed(i)),double(A.DS(i)),double(A.Jmax(i)));end;s=strjoin(parts,' | ');
end
function [available,amb,lab,tr]=resolve_terminal_fit_local(fits,trims,CFG)
available=false;amb=false;lab=[];tr=[];if isempty(fits),return,end;good=~cellfun(@isempty,fits);fits=fits(good);trims=trims(good);if isempty(fits),return,end;if numel(fits)==1,available=true;lab=fits{1};tr=trims{1};return,end
refL=fits{1};refT=trims{1};equiv=true;for k=2:numel(fits),if ~isequal(refT,trims{k}),equiv=false;break,end;keep=~refT;if sum(keep)>1,ar=adjusted_rand_index_v1(double(refL(keep)),double(fits{k}(keep)));else,ar=1;end;if ~isfinite(ar)||abs(ar-1)>CFG.tol,equiv=false;break,end;end
if equiv,available=true;lab=refL;tr=refT;else,amb=true;end
end
function y=lognormpdf_local(x,mu,v),y=-.5*(log(2*pi*v)+(x-mu).^2./v);end
function q=quantile_simple_local(x,p),x=sort(double(x(:)));n=numel(x);h=1+(n-1)*p;lo=floor(h);hi=ceil(h);if lo==hi,q=x(lo);else,q=x(lo)+(h-lo)*(x(hi)-x(lo));end,end
function a=project_grid_local(w,grid),[~,i]=min(abs(double(grid)-double(w)));a=grid(i);end
function Cache=empty_tfkm_cache_local(),Cache=struct('keys',strings(0,1),'pool',{cell(0,1)},'hier',{cell(0,1)});end
function Cache=empty_trkm_cache_local(),Cache=struct('keys',strings(0,1),'pool',{cell(0,1)});end

function X=regenerate_X_only_local(catrow,data_seed,root,method)
% This is the only helper allowed to read DGP-generative catalog fields.
canon=lock_vendor_local(root,method);f=string(which('generate_fixed_dataset_v1'));assert(startsWith(f,canon,'IgnoreCase',true),'Generator path lock failed.');cfg=table2struct(catrow);sim=generate_fixed_dataset_v1(cfg,data_seed);X=sim.X;clear sim cfg
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";
for c=cands(:)'
    orig=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
    hier=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
    if isfile(fullfile(orig,'study.mat')) && isfolder(hier),root=c;outroot=string(hier);S=load(fullfile(orig,'study.mat'));break,end
end
assert(strlength(root)>0,'Original fresh study + HIER4 correction-rerun output not found.');
end
function sourcecheck_tfkm_local(root)
canon=lock_vendor_local(root,"REG_TFKM_HD");
req=["generate_fixed_dataset_v1","regularized_trimmed_fit_corrected_v1","make_robust_initialization_v1","adjusted_rand_index_v1","run_hd_multistart_top_hier_stable4_unseen_v1","select_tfkm_hier_stable4_unseen_v1"];
paths=strings(numel(req),1);
for k=1:numel(req),f=string(which(char(req(k))));assert(strlength(f)>0,'Missing %s.',req(k));assert(startsWith(f,canon,'IgnoreCase',true),'SOURCE LOCK FAIL %s -> %s',req(k),f);paths(k)=f;end
h=paths(req=="run_hd_multistart_top_hier_stable4_unseen_v1");assert(contains(h,[filesep 'hier4' filesep],'IgnoreCase',true),'HIER4 wrapper is not vendor/hier4.');
fprintf('\n============================================================\n');fprintf('TFKM HIER4 BLIND FINALIZER SOURCE CHECK - TRUTH SEALED\n');fprintf('============================================================\n');disp(table(req(:),paths,'VariableNames',{'function_name','resolved_path'}));fprintf('STATUS: HIER4_BLIND_SOURCE_LOCK_CONFIRMED\n');
end
function canon=lock_vendor_local(root,method)
tf=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));tr=string(fullfile(root,'TRKM_TARGETED_MAX_BACKWARD_V1_1','vendor','TFKM_TRKM_UNSEEN_FINAL_v1'));
if method=="REG_TFKM_HD"
    assert(isfolder(tf),'TFKM vendor missing.');if isfolder(tr),safe_rmpath_tree_local(tr);end;addpath(genpath(char(tf)),'-begin');canon=tf;
else
    assert(isfolder(tr),'TRKM vendor missing.');if isfolder(tf),safe_rmpath_tree_local(tf);end;addpath(genpath(char(tr)),'-begin');canon=tr;
end;rehash path;
end
function safe_rmpath_tree_local(tree)
tree=char(tree);dirs=string(strsplit(genpath(tree),pathsep));dirs=dirs(strlength(dirs)>0);cur=string(strsplit(path,pathsep));for i=1:numel(dirs),if any(strcmpi(cur,dirs(i))),rmpath(char(dirs(i)));end,end
end
function assert_no_truth_columns_local(T)
if isempty(T),return,end;n=lower(string(T.Properties.VariableNames));bad=contains(n,"truth")|contains(n,"posthoc")|contains(n,"outlier")|ismember(n,["r0","g0","c0","alpha0","contamination","rho","delta","p_over_n"]);assert(~any(bad),'TRUTH FIREWALL FAIL: forbidden output column(s): %s',strjoin(n(bad),', '));
end
function T=empty_edge_export_local(),T=table('Size',[0 3],'VariableTypes',{'string','double','double'},'VariableNames',{'method','scenario_id','replicate'});end
function T=empty_repro_export_local(),T=table('Size',[0 3],'VariableTypes',{'string','double','double'},'VariableNames',{'method','scenario_id','replicate'});end

function selftest_local()
assert(isequal(predecessor_local([5 6]),[4 5]));assert(isequal(predecessor_local([4 6]),[4 5]));assert(isequal(predecessor_local([4 5]),[3 4]));assert(isempty(predecessor_local([1 2])));
M=componentwise_maximal_local([3 4;3 5;4 5]);assert(isequal(M,[4 5]));M=componentwise_maximal_local([4 5;3 6]);assert(size(M,1)==2);
q=quality_high_local([1;2;3]);assert(all(abs(q-[1/6;3/6;5/6])<1e-12));
CFG=struct('pos_grid',[.025 .05 .075 .10 .125 .15 .20 .25 .30],'tol',1e-12);rng(1);qq=[exp(.2*randn(90,1));exp(4+.2*randn(10,1))];[e,a,w,d,s]=radial_gmm_bic_local(qq,CFG);assert(e&&isfinite(a)&&w>0&&d>0&&s>0);
fprintf('SELFTEST PASS: HIER4 blind finalization geometry + frozen contracts.\n');
end
