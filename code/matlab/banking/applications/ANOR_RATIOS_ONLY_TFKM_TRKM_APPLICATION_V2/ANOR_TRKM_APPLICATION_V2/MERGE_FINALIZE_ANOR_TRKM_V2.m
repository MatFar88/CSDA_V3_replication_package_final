function OUT = MERGE_FINALIZE_ANOR_TRKM_V2()
%MERGE_FINALIZE_ANOR_TRKM_V2 Merge 150 ANOR ratios-only fits and apply frozen TRKM finalizer.
%
% Exact operational architecture:
% archived RV path -> archived A_.01 branch + radial GMM override ->
% RAW EXPL_NATIVE_010 upper-tail anchors -> propagated SUPPORT+RV_MIN ->
% mutual-retained ARI >= .90 gate.
%
% No bank names enter selection. If no positive trimming action is selected,
% the structural pair is retained but no primary terminal partition is created.

root=fileparts(mfilename('fullpath'));
CFG=anor_trkm_config_v2(root);
anor_trkm_path_lock_v2(root);
assert(isfile(CFG.input_file),'Build the input first.');
S=load(CFG.input_file,'X','original_row_index','protocol_version','pair_seed','options');
X=double(S.X);options=S.options;original_row_index=double(S.original_row_index(:));
assert(isequal(size(X),[CFG.n CFG.p]),'Input size drift.');
assert(string(S.protocol_version)==CFG.version,'Protocol version drift.');
assert(double(S.pair_seed)==CFG.pair_seed,'Pair-seed drift.');

[C,Labels,Trim]=load_cube_local(CFG);
assert(height(C)==150,'Expected 150 candidate rows.');
assert(all(logical(C.valid)),'At least one candidate is invalid. Investigate before finalization.');
assert(size(Labels,2)==150 && size(Trim,2)==150,'Stored cube shape mismatch.');

if ~isfolder(CFG.final_root),mkdir(CFG.final_root);end
[Prediction,Edges,Repro,FinalLabels,FinalTrim,FinalU,Audit]=finalize_trkm_local(X,C,Labels,Trim,options,CFG);

save(fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_FINAL.mat'), ...
    'Prediction','Edges','Repro','FinalLabels','FinalTrim','FinalU','Audit','-v7.3');
writetable(Prediction,fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_SELECTION.csv'));
if isempty(Edges),writetable(empty_edge_export_local(),fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_EDGES.csv'));else,writetable(Edges,fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_EDGES.csv'));end
if isempty(Repro),writetable(empty_repro_export_local(),fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_REPRODUCTION.csv'));else,writetable(Repro,fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_REPRODUCTION.csv'));end

if ~isempty(FinalLabels)
    P=table((1:CFG.n)',original_row_index,double(FinalLabels(:)),logical(FinalTrim(:)), ...
        'VariableNames',{'bank_index','original_365_row','cluster','trimmed'});
    writetable(P,fullfile(CFG.final_root,'ANOR_TRKM_PRIMARY_PARTITION_ANONYMOUS.csv'));
end

% Anonymous selector diagnostics useful for later interpretation.
if ~isempty(Audit.archived_per_alpha)
    A=array2table(Audit.archived_per_alpha,'VariableNames',{'alpha','r','g','RV_CV'});
    writetable(A,fullfile(CFG.final_root,'ANOR_TRKM_ARCHIVED_PER_ALPHA.csv'));
end
if ~isempty(Audit.raw_per_alpha)
    R=array2table(Audit.raw_per_alpha,'VariableNames',{'alpha','r','g','RV_CV'});
    writetable(R,fullfile(CFG.final_root,'ANOR_TRKM_RAW_PER_ALPHA.csv'));
end
Rad=table((1:CFG.n)',original_row_index,double(Audit.radial_q(:)), ...
    'VariableNames',{'bank_index','original_365_row','radial_q'});
writetable(Rad,fullfile(CFG.final_root,'ANOR_TRKM_RADIAL_SCORES_ANONYMOUS.csv'));

fprintf('\n============================================================\n');
fprintf('ANOR TRKM PRIMARY FINALIZATION COMPLETE\n');
fprintf('============================================================\n');
disp(Prediction);
if ~logical(Prediction.final_fit_available)
    if logical(Prediction.terminal_unique) && ~logical(Prediction.action_ambiguous)
        fprintf('\nNo primary terminal partition is defined by the selector.\n');
        fprintf('A pre-specified secondary descriptive refit is admissible. Run:\n');
        fprintf('  RUN_ANOR_TRKM_DESCRIPTIVE_REFIT_V2\n');
    else
        fprintf('\nNo unique structural state is available for a descriptive refit.\n');
    end
end
OUT=struct('prediction',Prediction,'final_root',CFG.final_root,'primary_fit_available',logical(Prediction.final_fit_available));
end

% =========================================================================
% Merge candidate cube
% =========================================================================
function [C,Labels,Trim]=load_cube_local(CFG)
C=table();Labels=zeros(CFG.n,150,'uint16');Trim=false(CFG.n,150);
for k=1:150
    f=fullfile(CFG.candidate_root,sprintf('candidate_%03d.mat',k));
    d=fullfile(CFG.candidate_root,sprintf('candidate_%03d.COMPLETE.txt',k));
    assert(isfile(f)&&isfile(d),'Candidate %d is incomplete.',k);
    Q=load(f,'Row','Labels','Trim','CandidateMeta');
    assert(isfield(Q,'CandidateMeta'),'Candidate %d lacks application-contract metadata.',k);
    assert(string(Q.CandidateMeta.protocol_version)==CFG.version && double(Q.CandidateMeta.pair_seed)==CFG.pair_seed && ...
        double(Q.CandidateMeta.n)==CFG.n && double(Q.CandidateMeta.p)==CFG.p && isequal(double(Q.CandidateMeta.drop_raw_rows(:))',CFG.drop_raw_rows), ...
        'Candidate %d application-contract mismatch.',k);
    assert(double(Q.Row.candidate_index)==k,'Candidate index mismatch in file %d.',k);
    assert(double(Q.Row.seed)==CFG.pair_seed && double(Q.Row.n)==CFG.n && double(Q.Row.p)==CFG.p,'Candidate %d seed/dimension mismatch.',k);
    assert(string(Q.Row.method)=="REG_TRKM_HD" && string(Q.Row.dataset)==CFG.dataset,'Candidate %d method/dataset mismatch.',k);
    if isempty(C),C=Q.Row;else,C=[C;Q.Row];end %#ok<AGROW>
    Labels(:,k)=uint16(Q.Labels(:));Trim(:,k)=logical(Q.Trim(:));
end
C=sortrows(C,'candidate_index');assert(isequal(double(C.candidate_index(:))',(1:150)),'Candidate ordering mismatch.');
end

% =========================================================================
% Exact TRKM finalization adapted from the 7-Oct-2026 blind finalizer.
% =========================================================================
function [Pred,Edges,Repro,FinalLabels,FinalTrim,FinalU,Audit]=finalize_trkm_local(X,C,Lstored,Tstored,options,CFG)
[archW,archModal]=trkm_archived_structure_local(C,CFG);
Pmodal=sortrows(C(double(C.r)==archModal(1)&double(C.g)==archModal(2),:),'alpha');
assert(height(Pmodal)==numel(CFG.alpha_grid),'Incomplete archived modal path.');
rv=arrayfun(@(k)stored_rv_local(Pmodal(k,:)),(1:height(Pmodal))');aa=double(Pmodal.alpha);
rvmax=max(rv);deltaRV=rvmax-min(rv);ii=find(rv>=rvmax-CFG.trkm_onset_tol-CFG.tol,1,'first');assert(~isempty(ii),'Archived onset undefined.');
A01=aa(ii);archEmit=(A01>0)&&(deltaRV>=CFG.trkm_delta_cutoff-CFG.tol);archA=NaN;if archEmit,archA=A01;end

q=radial_q_local(X);[radEmit,radA,wHi,dBIC,sep]=radial_gmm_bic_local(q,CFG);
if radEmit
    emitted=true;aHat=radA;actionSource="RADIAL_GMM_OVERRIDE";
elseif archEmit
    emitted=true;aHat=archA;actionSource="ARCHIVED_FALLBACK";
else
    emitted=false;aHat=NaN;actionSource="NO_ACTION";
end

Edges=table();Repro=table();FinalLabels=[];FinalTrim=[];FinalU=[];
fitAvailable=false;fitAmbiguous=false;subspaceAvailable=false;subspaceAmbiguous=false;nEdges=0;nAccepts=0;rawW=zeros(0,4);anchors=zeros(0,2);
if ~emitted
    terminalSet=archModal;terminalUnique=true;terminalR=archModal(1);terminalG=archModal(2);terminalStatus="NO_ACTION_ARCHIVED_MODAL";
else
    rawW=trkm_raw_winners_local(C,CFG);
    anchors=componentwise_maximal_local(rawW(rawW(:,1)>=aHat-CFG.tol,2:3));
    assert(~isempty(anchors),'Emitted TRKM unit has empty RAW tail envelope.');
    Cache=empty_trkm_cache_local();terms=zeros(0,2);fits=cell(size(anchors,1),1);trims=cell(size(anchors,1),1);Us=cell(size(anchors,1),1);
    for k=1:size(anchors,1)
        [term,E,R,lab,tr,U,Cache,ne,na]=trkm_path_local(X,C,Lstored,Tstored,options,CFG,anchors(k,:),aHat,k,Cache);
        terms(end+1,:)=term;Edges=[Edges;E];Repro=[Repro;R];fits{k}=lab;trims{k}=tr;Us{k}=U;nEdges=nEdges+ne;nAccepts=nAccepts+na; %#ok<AGROW>
    end
    terminalSet=unique(terms,'rows');terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    [fitAvailable,fitAmbiguous,FinalLabels,FinalTrim,subspaceAvailable,subspaceAmbiguous,FinalU]=resolve_terminal_fit_local(fits,trims,Us,CFG);
    if ~terminalUnique,terminalStatus="AMBIGUOUS_TERMINAL";elseif fitAmbiguous,terminalStatus="UNIQUE_STRUCTURE_FIT_AMBIGUOUS";else,terminalStatus="UNIQUE_TERMINAL";end
end

Pred=table("TRKM",CFG.dataset,string(actionSource),logical(emitted),double(aHat),false, ...
    rg_text_local(archW(:,2:3)),rg_text_local(archModal),rg_text_local(anchors),size(anchors,1), ...
    logical(archEmit),double(archA),double(deltaRV),logical(radEmit),double(radA),double(wHi),double(dBIC),double(sep), ...
    logical(terminalUnique),rg_text_local(terminalSet),double(terminalR),double(terminalG),string(terminalStatus), ...
    logical(fitAvailable),logical(fitAmbiguous),logical(subspaceAvailable),logical(subspaceAmbiguous),double(nEdges),double(nAccepts), ...
    'VariableNames',{'method','dataset','action_status','alpha_emitted','alpha_hat','action_ambiguous', ...
    'archived_per_alpha_winner_path','archived_modal_set','anchor_set','n_anchors', ...
    'archived_emit','archived_alpha_hat','archived_delta_RV','radial_emit','radial_alpha_hat','radial_high_weight','radial_delta_BIC','radial_separation_sd', ...
    'terminal_unique','terminal_set','terminal_r','terminal_g','terminal_status','final_fit_available','final_fit_ambiguous', ...
    'final_subspace_available','final_subspace_ambiguous','N_edges','N_accepts'});
Audit=struct('archived_per_alpha',archW,'archived_modal',archModal,'raw_per_alpha',rawW,'anchors',anchors,'radial_q',q, ...
    'pair_seed',CFG.pair_seed,'protocol_version',CFG.version,'frozen_multiple_anchor_rule',true);
end

function [term,Edges,Repro,FinalLabels,FinalTrim,FinalU,Cache,nEdges,nAccepts]=trkm_path_local(X,C,Lstored,Tstored,options,CFG,anchor,a,anchorIndex,Cache)
cur=anchor;propIdx=NaN;Edges=table();Repro=table();nEdges=0;nAccepts=0;FinalLabels=[];FinalTrim=[];FinalU=[];
while has_pred_local(cur)
    nEdges=nEdges+1;pred=predecessor_local(cur);
    [PF,Cache,rrF]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,cur(1),cur(2),a,"FROM");
    [PT,Cache,rrT]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,pred(1),pred(2),a,"TO");
    Repro=[Repro;rrF;rrT]; %#ok<AGROW>
    CM=cross_matrix_local(PF,PT);force=NaN;fromProp=false;if isfinite(propIdx),force=propIdx;fromProp=true;end
    [ii,jj,status]=trkm_choose_support_pair_local(PF,PT,CM,force,CFG);edge=CM(ii,jj);accept=isfinite(edge)&&edge>=CFG.ari_cutoff-CFG.tol;action="STOP";if accept,action="ACCEPT";end
    Edges=[Edges;table("TRKM",CFG.dataset,double(anchorIndex),double(a),double(nEdges),cur(1),cur(2),pred(1),pred(2),ii,jj, ...
        string(PF.family(ii)),string(PT.family(jj)),double(PF.start_id(ii)),double(PT.start_id(jj)),fromProp,string(status),double(edge),string(action), ...
        'VariableNames',{'method','dataset','anchor_index','alpha_ref','edge_index','r_from','g_from','r_to','g_to', ...
        'from_index','to_index','from_family','to_family','from_start_id','to_start_id','from_was_propagated','selector_status','selected_edge_ARI','selected_action'})]; %#ok<AGROW>
    FinalLabels=double(PF.labels(:,ii));FinalTrim=logical(PF.trim(:,ii));FinalU=double(PF.U{ii});
    if ~accept,break;end
    nAccepts=nAccepts+1;cur=pred;propIdx=jj;FinalLabels=double(PT.labels(:,jj));FinalTrim=logical(PT.trim(:,jj));FinalU=double(PT.U{jj});
end
if isempty(FinalLabels)
    [P,Cache,rr]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,cur(1),cur(2),a,"TERMINAL_BOUNDARY");Repro=[Repro;rr]; %#ok<AGROW>
    ii=P.chosen;FinalLabels=double(P.labels(:,ii));FinalTrim=logical(P.trim(:,ii));FinalU=double(P.U{ii});
end
term=cur;
end

function [P,Cache,RR]=get_trkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,r,g,a,state)
key=string(sprintf('%d|%d|%.6f',r,g,a));k=find(Cache.keys==key,1);
if isempty(k)
    old=candidate_row_local(C,r,g,a,CFG);P=generate_trkm_pool_local(X,r,g,a,options,CFG);reproduce_trkm_local(P,old,Lstored,Tstored,CFG,r,g,a);
    Cache.keys(end+1,1)=key;Cache.pool{end+1,1}=P;k=numel(Cache.keys);
else
    P=Cache.pool{k};
end
old=candidate_row_local(C,r,g,a,CFG);
RR=table("TRKM",CFG.dataset,string(state),r,g,a,string(P.family(P.chosen)),double(P.start_id(P.chosen)),double(P.cv(P.chosen)),stored_rv_local(old),double(P.objective(P.chosen)),double(old.objective),true, ...
    'VariableNames',{'method','dataset','state','r','g','alpha','selected_family','selected_start_id','new_score','stored_score','new_objective','stored_objective','reproduction_pass'});
end

function P=generate_trkm_pool_local(X,r,g,a,options,CFG)
N=sum(CFG.options.hd_trkm_allocation);n=size(X,1);valid=false(N,1);labels=zeros(n,N,'uint16');trim=false(n,N);objective=nan(N,1);family=strings(N,1);start_id=nan(N,1);U=cell(N,1);z=0;
for j=1:numel(CFG.options.hd_init_names)
    for h=1:CFG.options.hd_trkm_allocation(j)
        z=z+1;seed=CFG.pair_seed+1000*j+h;family(z)=CFG.options.hd_init_names(j);start_id(z)=h;
        try
            rng(seed,'twister');init=make_robust_initialization_v1(X,g,r,a,CFG.options.hd_init_names(j),seed);
            q=regularized_trimmed_fit_corrected_v1(X,init,r,g,a,"REG_TRKM",options.max_iterations,options.tauW,options.tauT);
            ok=all(isfinite(q.labels))&&isfinite(q.objective)&&size(q.U,2)==r;valid(z)=ok;
            if ok,labels(:,z)=uint16(q.labels(:));trim(:,z)=logical(q.trim(:));objective(z)=double(q.objective);U{z}=double(q.U);end
        catch
            valid(z)=false;
        end
    end
end
idx=find(valid);assert(~isempty(idx),'No valid TRKM starts at (%d,%d,%.3f).',r,g,a);objrank=objective_rank_local(objective,valid);
pool=struct('valid',valid,'family',family,'start_id',start_id,'objective',objective,'objective_rank',objrank,'labels',labels,'trim',trim,'trim_stability_weight',options.trim_stability_weight);
pool.current_candidate_index=legacy_candidate_local(pool,idx,options.trim_stability_weight);
sel=select_robust_view_v3(pool,X,g,r,"RV_CV");chosen=sel.chosen;score=selector_vector_local(sel,'score',N);cv=selector_vector_local(sel,'cv_score',N);if all(~isfinite(cv)),cv=score;end
P=struct('N',N,'valid',valid,'family',family,'start_id',start_id,'objective',objective,'objective_rank',objrank,'labels',labels,'trim',trim,'U',{U},'score',score,'cv',cv,'chosen',chosen,'r',r,'g',g,'alpha',a);
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

function reproduce_trkm_local(P,old,Lstored,Tstored,CFG,r,g,a)
assert(string(P.family(P.chosen))==string(old.selected_init_family),'TRKM family reproduction failure.');
assert(abs(double(P.cv(P.chosen))-stored_rv_local(old))<=CFG.trkm_score_tol,'TRKM RV reproduction failure.');
assert(abs(double(P.objective(P.chosen))-double(old.objective))<=CFG.obj_rel_tol*max(1,abs(double(old.objective))),'TRKM objective reproduction failure.');
idx=double(old.candidate_index);tr=logical(P.trim(:,P.chosen));assert(isequal(tr,logical(Tstored(:,idx))),'TRKM trim reproduction failure.');
keep=~tr;ar=adjusted_rand_index_v1(double(P.labels(keep,P.chosen)),double(Lstored(keep,idx)));assert(isfinite(ar)&&abs(ar-1)<=1e-12,'TRKM labels reproduction failure (%d,%d,%.3f).',r,g,a);
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
% Common helpers
% =========================================================================
function row=candidate_row_local(C,r,g,a,CFG)
row=C(double(C.r)==r&double(C.g)==g&abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(row)==1,'Candidate row missing/duplicated (%d,%d,%.3f).',r,g,a);
end
function Cmat=cross_matrix_local(P,Q)
Cmat=nan(P.N,Q.N);for i=1:P.N,if ~P.valid(i),continue,end;for j=1:Q.N,if ~Q.valid(j),continue,end;keep=~logical(P.trim(:,i))&~logical(Q.trim(:,j));if sum(keep)>1,Cmat(i,j)=adjusted_rand_index_v1(double(P.labels(keep,i)),double(Q.labels(keep,j)));end,end,end
end
function q=quality_high_local(x)
x=double(x(:));q=nan(size(x));ok=isfinite(x);vals=x(ok);if isempty(vals),return,end;for i=find(ok)',q(i)=(sum(vals<x(i))+.5*sum(vals==x(i)))/numel(vals);end
end
function q=quality_high_mask_local(x,valid)
x=double(x(:));valid=logical(valid(:));q=nan(size(x));idx=find(valid&isfinite(x));assert(~isempty(idx),'No finite quality values.');q(idx)=quality_high_local(x(idx));
end
function r=objective_rank_local(obj,valid)
r=nan(numel(obj),1);idx=find(valid);[~,o]=sort(obj(idx),'ascend');r(idx(o))=(1:numel(idx))';
end
function v=selector_vector_local(sel,name,N)
v=nan(N,1);if isfield(sel,name),x=double(sel.(name)(:));v(1:min(N,numel(x)))=x(1:min(N,numel(x)));end
end
function chosen=legacy_candidate_local(pool,valid_idx,w)
obj=pool.objective(valid_idx);[~,ord]=sort(obj,'ascend');ordered=valid_idx(ord);top_n=min(numel(ordered),max(3,ceil(numel(ordered)/4)));top=ordered(1:top_n);if numel(top)==1,chosen=top(1);return,end
S=zeros(numel(top));for a=1:numel(top),for b=a+1:numel(top),ta=logical(pool.trim(:,top(a)));tb=logical(pool.trim(:,top(b)));in=~ta&~tb;if sum(in)>1,ca=adjusted_rand_index_v1(double(pool.labels(in,top(a))),double(pool.labels(in,top(b))));if ~isfinite(ca),ca=0;end;else,ca=0;end;un=ta|tb;if any(un),tj=sum(ta&tb)/sum(un);else,tj=1;end;S(a,b)=ca+w*tj;S(b,a)=S(a,b);end,end;[~,m]=max(mean(S,2));chosen=top(m);
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
function [available,amb,lab,tr,Uavailable,Uamb,U]=resolve_terminal_fit_local(fits,trims,Us,CFG)
available=false;amb=false;lab=[];tr=[];Uavailable=false;Uamb=false;U=[];if isempty(fits),return,end;good=~cellfun(@isempty,fits);fits=fits(good);trims=trims(good);Us=Us(good);if isempty(fits),return,end
if numel(fits)==1,available=true;lab=fits{1};tr=trims{1};U=Us{1};Uavailable=~isempty(U);return,end
refL=fits{1};refT=trims{1};equiv=true;for k=2:numel(fits),if ~isequal(refT,trims{k}),equiv=false;break,end;keep=~refT;if sum(keep)>1,ar=adjusted_rand_index_v1(double(refL(keep)),double(fits{k}(keep)));else,ar=1;end;if ~isfinite(ar)||abs(ar-1)>CFG.tol,equiv=false;break,end;end
if equiv
    available=true;lab=refL;tr=refT;refU=Us{1};sameU=~isempty(refU);if sameU,P0=refU*refU';for k=2:numel(Us),if isempty(Us{k})||size(Us{k},2)~=size(refU,2)||norm(Us{k}*Us{k}'-P0,2)>1e-8,sameU=false;break,end,end,end
    if sameU,Uavailable=true;U=refU;else,Uamb=true;end
else
    amb=true;
end
end
function y=lognormpdf_local(x,mu,v),y=-.5*(log(2*pi*v)+(x-mu).^2./v);end
function q=quantile_simple_local(x,p),x=sort(double(x(:)));n=numel(x);h=1+(n-1)*p;lo=floor(h);hi=ceil(h);if lo==hi,q=x(lo);else,q=x(lo)+(h-lo)*(x(hi)-x(lo));end,end
function a=project_grid_local(w,grid),[~,i]=min(abs(double(grid)-double(w)));a=grid(i);end
function Cache=empty_trkm_cache_local(),Cache=struct('keys',strings(0,1),'pool',{cell(0,1)});end
function T=empty_edge_export_local(),T=table('Size',[0 2],'VariableTypes',{'string','string'},'VariableNames',{'method','dataset'});end
function T=empty_repro_export_local(),T=table('Size',[0 2],'VariableTypes',{'string','string'},'VariableNames',{'method','dataset'});end
