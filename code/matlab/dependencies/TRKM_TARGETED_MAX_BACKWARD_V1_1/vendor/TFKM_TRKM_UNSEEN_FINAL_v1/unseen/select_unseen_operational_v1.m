function R=select_unseen_operational_v1(output_root)
%SELECT_UNSEEN_OPERATIONAL_V1 Apply ALL frozen metaparameter rules truth-blind.
%
% TFKM:
%   HIER4 per-alpha global winner -> modal pair -> SUPPORTED_WEAK one-step.
%   tau* < .01 leverage; .01<=tau*<=.05 clean; tau*>.05 bridge ambiguous.
%   leverage alpha = unique argmax J(alpha>0); clean alpha=0; bridge abstain.
%
% TRKM:
%   per alpha max RV_CV -> exact-tie R2 -> parsimonious residual tie;
%   across alpha modal, frequency tie chooses larger pair.
%   no RV break = clean; break<.01 bridge ambiguous; break>=.01 leverage.
%   leverage alpha = FIRST_RV_LOSS high alpha; clean alpha=0; bridge abstain.
%
% THIS FUNCTION NEVER READS study.catalog OR SEALED TRUTH FILES.

tol=1e-12; root=char(string(output_root));
allfile=fullfile(root,'operational_selection','operational_candidates_all.csv');
if ~isfile(allfile), M=merge_unseen_operational_v1(root); else, M=readtable(allfile,'TextType','string'); end

for forbidden=["r0","g0","alpha0","contamination","ari_inlier","trim_f1"]
    assert(~ismember(forbidden,M.Properties.VariableNames), ...
        'Truth leakage: operational cube contains %s.',forbidden);
end

units=unique(M(:,{'scenario_id','replicate'}),'rows','stable');
TF=table(); TR=table();

% ======================================================================
% TFKM
% ======================================================================
for iu=1:height(units)
    sid=double(units.scenario_id(iu)); rep=double(units.replicate(iu));
    U=M(double(M.scenario_id)==sid & double(M.replicate)==rep & string(M.method)=="REG_TFKM_HD",:);
    assert(height(U)==150,'TFKM cube incomplete sid=%d rep=%d.',sid,rep);
    av=sort(unique(double(U.alpha)),'ascend'); W=zeros(numel(av),2); winScore=zeros(numel(av),1);
    for ia=1:numel(av)
        Q=U(abs(double(U.alpha)-av(ia))<=tol,:);
        s=double(Q.selector_score); mx=max(s); cand=find(abs(s-mx)<=tol);
        key=[double(Q.r(cand))+double(Q.g(cand)),double(Q.g(cand)),double(Q.r(cand))];
        [~,o]=sortrows(key,[1 2 3]); k=cand(o(1));
        W(ia,:)=[double(Q.r(k)),double(Q.g(k))]; winScore(ia)=s(k);
    end
    [rm,gm,mvotes,runnerVotes]=modal_small(W);
    weak=mvotes<=numel(av)/2;
    gmax=max(double(U.g)); isModal=W(:,1)==rm & W(:,2)==gm; ds=W(~isModal,:);
    if isempty(ds)
        boundary=false;
    else
        simpler=ds(:,1)<=rm & ds(:,2)<=gm & (ds(:,1)<rm | ds(:,2)<gm);
        boundary=(gm==gmax) && all(simpler);
    end
    [pr,pg,canBack]=canonical_pred(rm,gm,U);
    predSupport=sum(W(:,1)==pr & W(:,2)==pg);
    supportedWeak=weak && canBack && predSupport>0;
    trigger=boundary || supportedWeak;
    rh=rm; gh=gm;
    if trigger && canBack, rh=pr; gh=pg; end

    P=U(double(U.r)==rh & double(U.g)==gh,:); P=sortrows(P,'alpha','ascend');
    R2=double(P.R2); R2max=max(R2); loss=1-R2./R2max; tauStar=max(loss);
    if tauStar<.01-tol
        regime="LEVERAGE";
    elseif tauStar<=.05+tol
        regime="CLEAN";
    else
        regime="BRIDGE_AMBIGUOUS";
    end

    PJ=P(double(P.alpha)>tol,:); J=double(PJ.hier_cross_j_score); Jmax=max(J); im=find(abs(J-Jmax)<=tol);
    Junique=numel(im)==1; alphaJ=NaN; if Junique, alphaJ=double(PJ.alpha(im)); end
    if regime=="CLEAN"
        ah=0; identified=true; reason="CLEAN_ALPHA_ZERO";
    elseif regime=="LEVERAGE" && Junique
        ah=alphaJ; identified=true; reason="UNIQUE_J_MAX";
    elseif regime=="LEVERAGE"
        ah=NaN; identified=false; reason="TIED_J_MAX";
    else
        ah=NaN; identified=false; reason="BRIDGE_ALPHA_AMBIGUOUS";
    end

    TF=[TF;table(sid,rep,"REG_TFKM_HD",rm,gm,mvotes,runnerVotes,weak,boundary, ...
        pr,pg,predSupport,supportedWeak,trigger,rh,gh,tauStar,regime, ...
        Jmax,Junique,alphaJ,ah,identified,reason, ...
        'VariableNames',{'scenario_id','replicate','method','modal_r','modal_g', ...
        'modal_votes','runner_votes','weak_mode','boundary_trigger','pred_r','pred_g', ...
        'predecessor_support_count','supported_weak_trigger','structural_backstep', ...
        'selected_r','selected_g','tau_star','predicted_regime','J_max','J_unique_global_max', ...
        'J_max_alpha','alpha_hat','alpha_identified','alpha_reason'})]; %#ok<AGROW>
end

% ======================================================================
% TRKM
% ======================================================================
for iu=1:height(units)
    sid=double(units.scenario_id(iu)); rep=double(units.replicate(iu));
    U=M(double(M.scenario_id)==sid & double(M.replicate)==rep & string(M.method)=="REG_TRKM_HD",:);
    assert(height(U)==150,'TRKM cube incomplete sid=%d rep=%d.',sid,rep);
    av=sort(unique(double(U.alpha)),'ascend'); W=zeros(numel(av),2);
    for ia=1:numel(av)
        Q=U(abs(double(U.alpha)-av(ia))<=tol,:);
        k=pick_rv_r2_lex(Q,tol); W(ia,:)=[double(Q.r(k)),double(Q.g(k))];
    end
    [rh,gh,mvotes,runnerVotes]=modal_large(W);
    P=U(double(U.r)==rh & double(U.g)==gh,:); P=sortrows(P,'alpha','descend');
    a=double(P.alpha); rv=double(P.native_score); d=rv(1:end-1)-rv(2:end);
    k=find(d>tol,1,'first');
    if isempty(k)
        broke=false; breakLoss=0; breakHigh=NaN; breakLow=NaN;
        regime="CLEAN"; ah=min(a); identified=true; reason="NO_RV_BREAK";
    else
        broke=true; breakLoss=d(k); breakHigh=a(k); breakLow=a(k+1);
        if breakLoss>=.01-tol
            regime="LEVERAGE"; ah=breakHigh; identified=true; reason="RV_BREAK_AT_LEAST_001";
        else
            regime="BRIDGE_AMBIGUOUS"; ah=NaN; identified=false; reason="RV_BREAK_BELOW_001";
        end
    end
    TR=[TR;table(sid,rep,"REG_TRKM_HD",rh,gh,mvotes,runnerVotes,broke, ...
        breakHigh,breakLow,breakLoss,regime,ah,identified,reason, ...
        'VariableNames',{'scenario_id','replicate','method','selected_r','selected_g', ...
        'modal_votes','runner_votes','rv_break_found','break_alpha_high','break_alpha_low', ...
        'break_loss_RV','predicted_regime','alpha_hat','alpha_identified','alpha_reason'})]; %#ok<AGROW>
end

% Normalize common columns for combined table.
TF.rv_break_found=false(height(TF),1); TF.break_alpha_high=nan(height(TF),1); TF.break_alpha_low=nan(height(TF),1); TF.break_loss_RV=nan(height(TF),1);
TR.tau_star=nan(height(TR),1); TR.J_max=nan(height(TR),1); TR.J_unique_global_max=false(height(TR),1); TR.J_max_alpha=nan(height(TR),1);
common={'scenario_id','replicate','method','selected_r','selected_g','predicted_regime','alpha_hat','alpha_identified','alpha_reason', ...
        'tau_star','J_max','J_unique_global_max','J_max_alpha','rv_break_found','break_alpha_high','break_alpha_low','break_loss_RV'};
Dec=[TF(:,common);TR(:,common)]; Dec=sortrows(Dec,{'method','scenario_id','replicate'});

outdir=fullfile(root,'operational_selection');
writetable(TF,fullfile(outdir,'tfkm_operational_detail.csv'));
writetable(TR,fullfile(outdir,'trkm_operational_detail.csv'));
writetable(Dec,fullfile(outdir,'operational_decisions.csv'));
selection_frozen_at=string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')); %#ok<NASGU>
rule_version="TFKM_TRKM_UNSEEN_FINAL_v1"; %#ok<NASGU>
save(fullfile(outdir,'SELECTION_FROZEN.mat'),'TF','TR','Dec','selection_frozen_at','rule_version','-v7.3');

fprintf('\n============================================================\n');
fprintf('UNSEEN OPERATIONAL SELECTION FROZEN -- TRUTH STILL SEALED\n');
fprintf('TFKM: SUPPORTED_WEAK + tau* [.01,.05] + J alpha\n');
fprintf('TRKM: RV_R2_LEX_CONSERVATIVE + .01 break rule\n');
fprintf('Decisions: %d TFKM + %d TRKM\n',height(TF),height(TR));
fprintf('NOW AND ONLY NOW may evaluate_unseen_truth_v1 be called.\n');
fprintf('============================================================\n');

R=struct('tfkm',TF,'trkm',TR,'decisions',Dec,'truth_opened',false, ...
    'rule_version',rule_version);
end

function k=pick_rv_r2_lex(Q,tol)
rv=double(Q.native_score); mx=max(rv); cand=find(abs(rv-mx)<=tol);
r2=double(Q.R2); mx2=max(r2(cand)); cand2=cand(abs(r2(cand)-mx2)<=tol);
key=[double(Q.r(cand2))+double(Q.g(cand2)),double(Q.g(cand2)),double(Q.r(cand2))];
[~,o]=sortrows(key,[1 2 3]); k=cand2(o(1));
end
function [r,g,count,runner]=modal_small(P)
[U,~,ic]=unique(P,'rows'); cnt=accumarray(ic,1); mx=max(cnt); cand=find(cnt==mx);
key=[U(cand,1)+U(cand,2),U(cand,2),U(cand,1)]; [~,o]=sortrows(key,[1 2 3]); w=cand(o(1));
r=U(w,1); g=U(w,2); count=cnt(w); other=cnt; other(w)=[]; if isempty(other),runner=0;else,runner=max(other);end
end
function [r,g,count,runner]=modal_large(P)
[U,~,ic]=unique(P,'rows'); cnt=accumarray(ic,1); mx=max(cnt); cand=find(cnt==mx);
key=[U(cand,1)+U(cand,2),U(cand,2),U(cand,1)]; [~,o]=sortrows(key,[-1 -2 -3]); w=cand(o(1));
r=U(w,1); g=U(w,2); count=cnt(w); other=cnt; other(w)=[]; if isempty(other),runner=0;else,runner=max(other);end
end
function [pr,pg,ok]=canonical_pred(r,g,V)
pr=r; pg=g; ok=false;
if g>r+1, pg=g-1; ok=true; elseif g==r+1 && r>1, pr=r-1; pg=g-1; ok=true; end
if ok, ok=any(double(V.r)==pr & double(V.g)==pg); end
end
