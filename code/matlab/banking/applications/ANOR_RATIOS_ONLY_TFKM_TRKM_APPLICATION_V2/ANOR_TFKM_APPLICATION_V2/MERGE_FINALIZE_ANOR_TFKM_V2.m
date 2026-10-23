function OUT = MERGE_FINALIZE_ANOR_TFKM_V2()
%MERGE_FINALIZE_ANOR_TFKM_V2 Merge all 150 fits and apply the locked TFKM finalizer.
%
% The selection logic is the paper's HIER4 -> modal/M+ -> action/alpha ->
% propagated SUPPORT+NATIVE -> mutual-retained ARI >= .90 procedure.
% No bank names enter selection. If no positive trimming action is selected,
% the structural pair is retained but no terminal partition is manufactured.

root=fileparts(mfilename('fullpath'));
CFG=anor_tfkm_config_v2(root);
anor_tfkm_path_lock_v2(root);
assert(isfile(CFG.input_file),'Build the input first.');
S=load(CFG.input_file,'X','original_row_index','protocol_version','pair_seed','options');
X=double(S.X); options=S.options; original_row_index=double(S.original_row_index(:));
assert(isequal(size(X),[CFG.n CFG.p]),'Input size drift.');
assert(string(S.protocol_version)==CFG.version,'Protocol version drift.');
assert(double(S.pair_seed)==CFG.pair_seed,'Pair seed drift.');

[C,Labels,Trim]=load_cube_local(CFG);
assert(height(C)==150,'Expected 150 candidate rows.');
assert(all(logical(C.valid)),'At least one candidate is invalid. Do not finalize before investigating invalid states.');
assert(size(Labels,2)==150 && size(Trim,2)==150,'Stored cube shape mismatch.');

if ~isfolder(CFG.final_root),mkdir(CFG.final_root);end
[Prediction,Edges,Repro,FinalLabels,FinalTrim,FinalU,Audit]= ...
    finalize_tfkm_local(X,C,Labels,Trim,options,CFG);

save(fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_FINAL.mat'), ...
    'Prediction','Edges','Repro','FinalLabels','FinalTrim','FinalU','Audit','-v7.3');
writetable(Prediction,fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_SELECTION.csv'));
if isempty(Edges)
    writetable(empty_edge_export_local(),fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_EDGES.csv'));
else
    writetable(Edges,fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_EDGES.csv'));
end
if isempty(Repro)
    writetable(empty_repro_export_local(),fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_REPRODUCTION.csv'));
else
    writetable(Repro,fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_REPRODUCTION.csv'));
end

if ~isempty(FinalLabels)
    P=table((1:CFG.n)',original_row_index,double(FinalLabels(:)),logical(FinalTrim(:)), ...
        'VariableNames',{'bank_index','original_365_row','cluster','trimmed'});
    writetable(P,fullfile(CFG.final_root,'ANOR_TFKM_PRIMARY_PARTITION_ANONYMOUS.csv'));
end

fprintf('\n============================================================\n');
fprintf('ANOR TFKM PRIMARY FINALIZATION COMPLETE\n');
fprintf('============================================================\n');
disp(Prediction);
if ~logical(Prediction.final_fit_available)
    if logical(Prediction.terminal_unique) && ~logical(Prediction.action_ambiguous)
        fprintf('\nNo primary terminal partition is defined by the selector.\n');
        fprintf('A pre-specified secondary descriptive refit is admissible. Run:\n');
        fprintf('  RUN_ANOR_TFKM_DESCRIPTIVE_REFIT_V2\n');
    else
        fprintf('\nNo unique structural state is available for a descriptive refit.\n');
    end
end

OUT=struct('prediction',Prediction,'final_root',CFG.final_root, ...
    'primary_fit_available',logical(Prediction.final_fit_available));
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
    assert(string(Q.CandidateMeta.protocol_version)==CFG.version && ...
        double(Q.CandidateMeta.pair_seed)==CFG.pair_seed && ...
        double(Q.CandidateMeta.n)==CFG.n && double(Q.CandidateMeta.p)==CFG.p && ...
        isequal(double(Q.CandidateMeta.drop_raw_rows(:))',CFG.drop_raw_rows), ...
        'Candidate %d application-contract mismatch.',k);
    assert(double(Q.Row.candidate_index)==k,'Candidate index mismatch in file %d.',k);
    assert(double(Q.Row.seed)==CFG.pair_seed && double(Q.Row.n)==CFG.n && double(Q.Row.p)==CFG.p, ...
        'Candidate %d seed/dimension mismatch.',k);
    assert(string(Q.Row.method)=="REG_TFKM_HD" && string(Q.Row.dataset)==CFG.dataset, ...
        'Candidate %d method/dataset mismatch.',k);
    if isempty(C),C=Q.Row;else,C=[C;Q.Row];end %#ok<AGROW>
    Labels(:,k)=uint16(Q.Labels(:)); Trim(:,k)=logical(Q.Trim(:));
end
C=sortrows(C,'candidate_index');
assert(isequal(double(C.candidate_index(:))',(1:150)),'Candidate ordering mismatch.');
end

% =========================================================================
% Exact TFKM finalization adapted from the blind fresh-study finalizer.
% =========================================================================
function [Pred,Edges,Repro,FinalLabels,FinalTrim,FinalU,Audit]=finalize_tfkm_local(X,C,Lstored,Tstored,options,CFG)
[W,modal,Mplus]=tfkm_structure_local(C,CFG);
assert(~isempty(Mplus),'TFKM M+ empty.');
A=table();
for k=1:size(Mplus,1)
    r=Mplus(k,1);g=Mplus(k,2);
    P=sortrows(C(double(C.r)==r&double(C.g)==g,:),'alpha');
    assert(height(P)==numel(CFG.alpha_grid),'Incomplete TFKM anchor path.');
    f=tfkm_features_local(P,CFG);
    emit=f.J_qualified && ((isfinite(f.Smed)&&f.Smed<CFG.tfkm_Smed_cutoff)||f.DS_signal);
    a=NaN;if emit,a=f.alpha_proposed;end
    A=[A;table(k,r,g,f.S0,f.Smed,f.DS,f.Jmax,f.alpha_proposed,emit,a, ...
        'VariableNames',{'anchor_index','r','g','S0','Smed','DS','Jmax','alpha_proposed','emit','alpha_hat'})]; %#ok<AGROW>
end

allEmit=all(logical(A.emit)); allNo=all(~logical(A.emit)); alphaAgree=false; aHat=NaN;
if allEmit
    aa=double(A.alpha_hat); alphaAgree=max(aa)-min(aa)<=CFG.tol; if alphaAgree,aHat=aa(1);end
end
if allNo
    actionStatus="NO_ACTION_CONCORDANT"; emitted=false; actionAmb=false;
elseif allEmit && alphaAgree
    actionStatus="EMIT_CONCORDANT"; emitted=true; actionAmb=false;
else
    actionStatus="ACTION_AMBIGUOUS"; emitted=false; actionAmb=true;
end

Edges=table();Repro=table();FinalLabels=[];FinalTrim=[];FinalU=[];
fitAvailable=false;fitAmbiguous=false;subspaceAvailable=false;subspaceAmbiguous=false;nEdges=0;nAccepts=0;
if actionAmb
    terminalSet=Mplus;terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;
    if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    terminalStatus="ACTION_AMBIGUOUS_NO_PATH";
elseif ~emitted
    terminalSet=Mplus;terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;
    if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    terminalStatus="NO_ACTION_STRUCTURE_RETAINED";
else
    Cache=empty_tfkm_cache_local();terms=zeros(0,2);fits=cell(size(Mplus,1),1);trims=cell(size(Mplus,1),1);Us=cell(size(Mplus,1),1);
    for k=1:size(Mplus,1)
        [term,E,R,lab,tr,U,Cache,ne,na]=tfkm_path_local(X,C,Lstored,Tstored,options,CFG,Mplus(k,:),aHat,k,Cache);
        terms(end+1,:)=term;Edges=[Edges;E];Repro=[Repro;R];fits{k}=lab;trims{k}=tr;Us{k}=U;nEdges=nEdges+ne;nAccepts=nAccepts+na; %#ok<AGROW>
    end
    terminalSet=unique(terms,'rows');terminalUnique=size(terminalSet,1)==1;terminalR=NaN;terminalG=NaN;
    if terminalUnique,terminalR=terminalSet(1,1);terminalG=terminalSet(1,2);end
    [fitAvailable,fitAmbiguous,FinalLabels,FinalTrim,subspaceAvailable,subspaceAmbiguous,FinalU]=resolve_terminal_fit_local(fits,trims,Us,CFG);
    if ~terminalUnique,terminalStatus="AMBIGUOUS_TERMINAL";elseif fitAmbiguous,terminalStatus="UNIQUE_STRUCTURE_FIT_AMBIGUOUS";else,terminalStatus="UNIQUE_TERMINAL";end
end

Pred=table("TFKM",CFG.dataset,string(actionStatus),logical(emitted),double(aHat),logical(actionAmb), ...
    rg_text_local(W(:,2:3)),rg_text_local(modal),rg_text_local(Mplus),size(Mplus,1),anchor_action_text_local(A), ...
    logical(terminalUnique),rg_text_local(terminalSet),double(terminalR),double(terminalG),string(terminalStatus), ...
    logical(fitAvailable),logical(fitAmbiguous),logical(subspaceAvailable),logical(subspaceAmbiguous),double(nEdges),double(nAccepts), ...
    'VariableNames',{'method','dataset','action_status','alpha_emitted','alpha_hat','action_ambiguous', ...
    'per_alpha_winner_path','modal_set','anchor_set','n_anchors','anchor_action_detail', ...
    'terminal_unique','terminal_set','terminal_r','terminal_g','terminal_status','final_fit_available','final_fit_ambiguous', ...
    'final_subspace_available','final_subspace_ambiguous','N_edges','N_accepts'});
Audit=struct('per_alpha',W,'modal',modal,'Mplus',Mplus,'anchor_action',A,'multiple_anchor_rule_fixed',true, ...
    'pair_seed',CFG.pair_seed,'protocol_version',CFG.version);
end

function [term,Edges,Repro,FinalLabels,FinalTrim,FinalU,Cache,nEdges,nAccepts]=tfkm_path_local(X,C,Lstored,Tstored,options,CFG,anchor,a,anchorIndex,Cache)
cur=anchor;propIdx=NaN;Edges=table();Repro=table();nEdges=0;nAccepts=0;FinalLabels=[];FinalTrim=[];FinalU=[];
while has_pred_local(cur)
    nEdges=nEdges+1;pred=predecessor_local(cur);
    [PF,HF,Cache,rrF]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,cur(1),cur(2),a,"FROM");
    [PT,HT,Cache,rrT]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,pred(1),pred(2),a,"TO");
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
    Edges=[Edges;table("TFKM",CFG.dataset,double(anchorIndex),double(a),double(nEdges),cur(1),cur(2),pred(1),pred(2), ...
        ii,jj,string(PF.family(ii)),string(PT.family(jj)),double(PF.start_id(ii)),double(PT.start_id(jj)),fromProp, ...
        double(EF.support_count(ii)),double(ET.support_count(jj)),double(EF.support_pct(ii)),double(ET.support_pct(jj)), ...
        double(EF.native_pct(ii)),double(ET.native_pct(jj)),double(edge),string(action), ...
        'VariableNames',{'method','dataset','anchor_index','alpha_ref','edge_index','r_from','g_from','r_to','g_to', ...
        'from_index','to_index','from_family','to_family','from_start_id','to_start_id','from_was_propagated', ...
        'from_support_count','to_support_count','from_support_pct','to_support_pct','from_native_pct','to_native_pct','selected_edge_ARI','selected_action'})]; %#ok<AGROW>
    FinalLabels=double(PF.labels(:,ii));FinalTrim=logical(PF.trim(:,ii));FinalU=double(PF.U{ii});
    if ~accept,break;end
    nAccepts=nAccepts+1;cur=pred;propIdx=jj;FinalLabels=double(PT.labels(:,jj));FinalTrim=logical(PT.trim(:,jj));FinalU=double(PT.U{jj});
end
if isempty(FinalLabels)
    [P,H,Cache,rr]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,cur(1),cur(2),a,"TERMINAL_BOUNDARY");Repro=[Repro;rr]; %#ok<AGROW>
    ii=double(H.chosen_index);FinalLabels=double(P.labels(:,ii));FinalTrim=logical(P.trim(:,ii));FinalU=double(P.U{ii});
end
term=cur;
end

function [P,H,Cache,RR]=get_tfkm_pool_local(Cache,X,C,Lstored,Tstored,options,CFG,r,g,a,state)
key=string(sprintf('%d|%d|%.6f',r,g,a));k=find(Cache.keys==key,1);
if isempty(k)
    P=generate_tfkm_pool_local(X,r,g,a,options,CFG);H=hier4_decompose_local(P,CFG);
    old=candidate_row_local(C,r,g,a,CFG);reproduce_tfkm_local(P,H,old,Lstored,Tstored,CFG,r,g,a);
    Cache.keys(end+1,1)=key;Cache.pool{end+1,1}=P;Cache.hier{end+1,1}=H;k=numel(Cache.keys);
else
    P=Cache.pool{k};H=Cache.hier{k};
end
old=candidate_row_local(C,r,g,a,CFG);
RR=table("TFKM",CFG.dataset,string(state),r,g,a,string(P.family(H.chosen_index)),double(P.start_id(H.chosen_index)), ...
    double(H.cross_H),double(old.selector_score),double(H.objective),double(old.objective),true, ...
    'VariableNames',{'method','dataset','state','r','g','alpha','selected_family','selected_start_id','new_score','stored_score','new_objective','stored_objective','reproduction_pass'});
end

function P=generate_tfkm_pool_local(X,r,g,a,options,CFG)
N=sum(CFG.options.hd_tfkm_allocation);n=size(X,1);valid=false(N,1);labels=zeros(n,N,'uint16');trim=false(n,N);objective=nan(N,1);family=strings(N,1);start_id=nan(N,1);U=cell(N,1);z=0;
for j=1:numel(CFG.options.hd_init_names)
    for h=1:CFG.options.hd_tfkm_allocation(j)
        z=z+1;seed=CFG.pair_seed+1000*j+h;family(z)=CFG.options.hd_init_names(j);start_id(z)=h;
        try
            rng(seed,'twister');init=make_robust_initialization_v1(X,g,r,a,CFG.options.hd_init_names(j),seed);
            q=regularized_trimmed_fit_corrected_v1(X,init,r,g,a,"REG_TFKM",options.max_iterations,options.tauW,options.tauT);
            ok=all(isfinite(q.labels))&&isfinite(q.objective)&&size(q.U,2)==r;valid(z)=ok;
            if ok,labels(:,z)=uint16(q.labels(:));trim(:,z)=logical(q.trim(:));objective(z)=double(q.objective);U{z}=double(q.U);end
        catch
            valid(z)=false;
        end
    end
end
assert(any(valid),'No valid TFKM starts at (%d,%d,%.3f).',r,g,a);
P=struct('N',N,'valid',valid,'labels',labels,'trim',trim,'objective',objective,'family',family,'start_id',start_id,'U',{U});
end

function H=hier4_decompose_local(P,CFG)
valid=logical(P.valid(:));fam=string(P.family(:));obj=double(P.objective(:));N=numel(valid);fw=nan(numel(CFG.options.hd_init_names),1);within=nan(N,1);
for f=1:numel(CFG.options.hd_init_names)
    idx=find(valid&fam==CFG.options.hd_init_names(f));if isempty(idx),continue,end
    M=pair_H_local(P,idx,CFG.trim_weight);s=mean(M,2,'omitnan');within(idx)=s;mx=max(s,[],'omitnan');cand=idx(abs(s-mx)<=CFG.tol);if numel(cand)>1,[~,m]=min(obj(cand));cand=cand(m);end;fw(f)=cand;
end
win=double(fw(isfinite(fw)));assert(~isempty(win),'No HIER4 family representatives.');final=nan(N,1);crossA=nan(N,1);crossJ=nan(N,1);
if numel(win)==1
    chosen=win(1);final(chosen)=within(chosen);crossA(chosen)=0;crossJ(chosen)=0;
else
    [A,J,M]=pair_components_local(P,win,CFG.trim_weight);s=mean(M,2,'omitnan');final(win)=s;crossA(win)=mean(A,2,'omitnan');crossJ(win)=mean(J,2,'omitnan');mx=max(s,[],'omitnan');cand=win(abs(s-mx)<=CFG.tol);if numel(cand)>1,[~,m]=min(obj(cand));chosen=cand(m);else,chosen=cand;end
end
H=struct('chosen_index',double(chosen),'cross_H',double(final(chosen)),'cross_A',double(crossA(chosen)),'cross_J',double(crossJ(chosen)), ...
    'within_H',double(within(chosen)),'objective',double(obj(chosen)),'family_winner_idx',fw,'score_vector',final,'within_vector',within);
end

function reproduce_tfkm_local(P,H,old,Lstored,Tstored,CFG,r,g,a)
assert(string(P.family(H.chosen_index))==string(old.selected_init_family),'TFKM family reproduction failure.');
assert(abs(double(H.cross_H)-double(old.selector_score))<=CFG.tfkm_score_tol,'TFKM score reproduction failure.');
assert(abs(double(H.cross_A)-double(old.hier_cross_ari_score))<=CFG.tfkm_score_tol,'TFKM ARI-component reproduction failure.');
assert(abs(double(H.cross_J)-double(old.hier_cross_j_score))<=CFG.tfkm_score_tol,'TFKM J-component reproduction failure.');
assert(abs(double(H.objective)-double(old.objective))<=CFG.obj_rel_tol*max(1,abs(double(old.objective))),'TFKM objective reproduction failure.');
idx=double(old.candidate_index);tr=logical(P.trim(:,H.chosen_index));assert(isequal(tr,logical(Tstored(:,idx))),'TFKM trim reproduction failure.');
keep=~tr;ar=adjusted_rand_index_v1(double(P.labels(keep,H.chosen_index)),double(Lstored(keep,idx)));assert(isfinite(ar)&&abs(ar-1)<=1e-12, ...
    'TFKM labels reproduction failure (%d,%d,%.3f).',r,g,a);
end

function [W,modal,Mplus]=tfkm_structure_local(C,CFG)
W=zeros(numel(CFG.alpha_grid),4);
for i=1:numel(CFG.alpha_grid)
    a=CFG.alpha_grid(i);Q=C(abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(Q)==15,'TFKM alpha slice not 15.');
    h=double(Q.selector_score);mx=max(h);ids=find(abs(h-mx)<=CFG.tol);key=[double(Q.r(ids))+double(Q.g(ids)),double(Q.g(ids)),double(Q.r(ids))];[~,o]=sortrows(key,[1 2 3]);k=ids(o(1));W(i,:)=[a,double(Q.r(k)),double(Q.g(k)),h(k)];
end
[u,~,ic]=unique(W(:,2:3),'rows');cnt=accumarray(ic,1);modal=sortrows(u(cnt==max(cnt),:),[1 2]);Mplus=componentwise_maximal_local(modal);
end

function f=tfkm_features_local(P,CFG)
a=double(P.alpha);s=double(P.hier_cross_ari_score);j=double(P.hier_cross_j_score);f.S0=s(1);f.Smed=median(s);f.DS=NaN;f.DS_defined=isfinite(s(1))&&s(1)>0&&all(isfinite(s));if f.DS_defined,f.DS=(max(s)-min(s))/s(1);end
pos=find(a>0);[f.Jmax,k]=max(j(pos));f.alpha_proposed=a(pos(k));f.J_qualified=isfinite(f.Jmax)&&f.Jmax>=CFG.tfkm_J_cutoff;f.DS_signal=f.DS_defined&&f.DS>=CFG.tfkm_DS_cutoff;
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

% =========================================================================
% Common helpers
% =========================================================================
function row=candidate_row_local(C,r,g,a,CFG)
row=C(double(C.r)==r&double(C.g)==g&abs(double(C.alpha)-a)<=CFG.tol,:);assert(height(row)==1,'Candidate row missing/duplicated (%d,%d,%.3f).',r,g,a);
end
function Cmat=cross_matrix_local(P,Q)
Cmat=nan(P.N,Q.N);for i=1:P.N,if ~P.valid(i),continue,end;for j=1:Q.N,if ~Q.valid(j),continue,end;keep=~logical(P.trim(:,i))&~logical(Q.trim(:,j));if sum(keep)>1,Cmat(i,j)=adjusted_rand_index_v1(double(P.labels(keep,i)),double(Q.labels(keep,j)));end,end,end
end
function M=pair_H_local(P,idx,w),[~,~,M]=pair_components_local(P,idx,w);end
function [A,J,H]=pair_components_local(P,idx,w)
idx=double(idx(:));m=numel(idx);A=zeros(m);J=zeros(m);for a=1:m,ia=idx(a);for b=a+1:m,ib=idx(b);ta=logical(P.trim(:,ia));tb=logical(P.trim(:,ib));keep=~ta&~tb;ar=0;if sum(keep)>1,ar=adjusted_rand_index_v1(double(P.labels(keep,ia)),double(P.labels(keep,ib)));if ~isfinite(ar),ar=0;end,end;un=ta|tb;if any(un),jac=sum(ta&tb)/sum(un);else,jac=1;end;A(a,b)=ar;A(b,a)=ar;J(a,b)=jac;J(b,a)=jac;end,end;H=A+w*J;
end
function q=quality_high_local(x)
x=double(x(:));q=nan(size(x));ok=isfinite(x);vals=x(ok);if isempty(vals),return,end;for i=find(ok)',q(i)=(sum(vals<x(i))+.5*sum(vals==x(i)))/numel(vals);end
end
function q=quality_low_local(x),q=quality_high_local(-double(x(:)));end
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
function [available,amb,lab,tr,Uavailable,Uamb,U]=resolve_terminal_fit_local(fits,trims,Us,CFG)
available=false;amb=false;lab=[];tr=[];Uavailable=false;Uamb=false;U=[];if isempty(fits),return,end;good=~cellfun(@isempty,fits);fits=fits(good);trims=trims(good);Us=Us(good);if isempty(fits),return,end
if numel(fits)==1,available=true;lab=fits{1};tr=trims{1};U=Us{1};Uavailable=~isempty(U);return,end
refL=fits{1};refT=trims{1};equiv=true;for k=2:numel(fits),if ~isequal(refT,trims{k}),equiv=false;break,end;keep=~refT;if sum(keep)>1,ar=adjusted_rand_index_v1(double(refL(keep)),double(fits{k}(keep)));else,ar=1;end;if ~isfinite(ar)||abs(ar-1)>CFG.tol,equiv=false;break,end;end
if equiv
    available=true;lab=refL;tr=refT;
    refU=Us{1};sameU=~isempty(refU);if sameU,P0=refU*refU';for k=2:numel(Us),if isempty(Us{k})||size(Us{k},2)~=size(refU,2)||norm(Us{k}*Us{k}'-P0,2)>1e-8,sameU=false;break,end,end,end
    if sameU,Uavailable=true;U=refU;else,Uamb=true;end
else
    amb=true;
end
end
function Cache=empty_tfkm_cache_local(),Cache=struct('keys',strings(0,1),'pool',{cell(0,1)},'hier',{cell(0,1)});end
function T=empty_edge_export_local(),T=table('Size',[0 2],'VariableTypes',{'string','string'},'VariableNames',{'method','dataset'});end
function T=empty_repro_export_local(),T=table('Size',[0 2],'VariableTypes',{'string','string'},'VariableNames',{'method','dataset'});end
