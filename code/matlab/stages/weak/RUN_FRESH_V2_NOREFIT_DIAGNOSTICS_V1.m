function OUT = RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1(documents_root)
%RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1
% Truth-open V2 development diagnostics on frozen fresh-study outputs.
% NO estimator fits/refits. NO V1 rule changes.
%
% Audits:
%   A1. TFKM action-feature distributions (Jmax, Smed, DS).
%   A2. TRKM radial-GMM vs archived fallback decomposition.
%   A3. (r,g) partition-equivalence sets at truth alpha0.
%
% Usage:
%   RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1;
%   RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1("C:\Users\name\Documents");
%   RUN_FRESH_V2_NOREFIT_DIAGNOSTICS_V1("SELFTEST");

if nargin>=1 && strcmpi(string(documents_root),"SELFTEST")
    assert(abs(rank_auc_local([0;1;2;3],[0;0;1;1])-1)<=1e-12);
    assert(abs(rank_auc_local([3;2;1;0],[0;0;1;1])-0)<=1e-12);
    assert(abs(rank_auc_local([1;1;1;1],[0;0;1;1])-.5)<=1e-12);
    assert(set_text_local([2 3;3 4])=="(2,3);(3,4)");
    fprintf('SELFTEST PASS: V2 no-refit diagnostic contracts.\n');
    OUT=struct('selftest_pass',true);return
end
if nargin<1,documents_root="";end
[root,outroot,S]=load_study_local(documents_root);
catalog=S.catalog;manifest=S.manifest;
assert(height(catalog)==41 && height(manifest)==1230,'Fresh study dimensions drifted.');

truthroot=fullfile(outroot,'fresh_truth_open_v1');
assert(isfile(fullfile(truthroot,'TRUTH_OPENED_AT.txt')),'Truth-open marker missing.');
postroot=fullfile(outroot,'fresh_posthoc_diagnostic_v1');
postfile=fullfile(postroot,'POSTHOC_UNIT_DIAGNOSTIC.csv');
assert(isfile(postfile),'Run broad posthoc diagnostic first: %s',postfile);
D=readtable(postfile,'TextType','string');
assert(height(D)==2460,'Expected 2460 posthoc unit rows.');

seal=fullfile(outroot,'FROZEN_PRETRUTH_20261010');
assert(isfolder(seal),'Pretruth snapshot missing.');
PT=readtable(fullfile(seal,'TFKM_HIER4_BLIND_FINAL_PREDICTIONS_FROZEN.csv'),'TextType','string');
PR=readtable(fullfile(seal,'TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv'),'TextType','string');
assert(height(PT)==1230 && height(PR)==1230,'Frozen prediction counts drifted.');

lock_vendor_local(root);
hierRoot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output');
tfBlind=fullfile(hierRoot,'blind_finalization_hier4_v1','TFKM');
trBlind=fullfile(outroot,'blind_finalization_v1','TRKM');
tfCube=fullfile(hierRoot,'fits','TFKM');
trCube=fullfile(outroot,'fits','TRKM');
for p=[string(tfBlind),string(trBlind),string(tfCube),string(trCube)]
    assert(isfolder(p),'Missing frozen source: %s',p);
end

outdiag=fullfile(outroot,'fresh_v2_diagnostic_suite_v1');
if ~isfolder(outdiag),mkdir(outdiag);end

fprintf('\n============================================================\n');
fprintf('FRESH V2 DEVELOPMENT DIAGNOSTICS V1 - NO REFITS\n');
fprintf('TRUTH OPEN | V1 FROZEN | DIAGNOSTICS ONLY\n');
fprintf('============================================================\n');

% -------------------------------------------------------------------------
% A1. TFKM action evidence, anchor-level and unit-level.
% -------------------------------------------------------------------------
Aanchor=table();Aunit=table();
for i=1:height(PT)
    sid=double(PT.scenario_id(i));rep=double(PT.replicate(i));
    ic=find(double(catalog.scenario_id)==sid,1);assert(~isempty(ic));c=catalog(ic,:);
    f=fullfile(tfBlind,sprintf('s%d',sid),sprintf('rep%02d',rep),'blind_final.mat');assert(isfile(f));
    Z=load(f,'Prediction','Audit');A=Z.Audit.anchor_action;P=Z.Prediction;assert(height(P)==1);
    lev=string(c.contamination)=="leverage";
    for k=1:height(A)
        jmargin=double(A.Jmax(k))-.50;
        smargin=.44-double(A.Smed(k));
        dmargin=double(A.DS(k))-.12;
        branch=max(smargin,dmargin);
        Aanchor=[Aanchor;table(sid,rep,string(c.panel),double(c.p_over_n),double(c.rho),double(c.delta),string(c.contamination),lev, ...
            double(A.anchor_index(k)),double(A.r(k)),double(A.g(k)),double(A.S0(k)),double(A.Smed(k)),double(A.DS(k)),double(A.Jmax(k)), ...
            double(A.alpha_proposed(k)),logical(A.emit(k)),jmargin,smargin,dmargin,branch, ...
            'VariableNames',{'scenario_id','replicate','panel','p_over_n','rho','delta','contamination_POSTHOC','leverage_POSTHOC', ...
            'anchor_index','anchor_r','anchor_g','S0','Smed','DS','Jmax','alpha_proposed','anchor_emit', ...
            'J_margin_050','Smed_margin_044','DS_margin_012','secondary_branch_margin'})]; %#ok<AGROW>
    end
    jj=double(A.Jmax);ss=double(A.Smed);dd=double(A.DS);ee=logical(A.emit);
    Aunit=[Aunit;table(sid,rep,string(c.panel),double(c.p_over_n),double(c.rho),double(c.delta),string(c.contamination),lev, ...
        height(A),min(jj),max(jj),median(jj),min(ss),max(ss),median(ss),safe_min_local(dd),safe_max_local(dd),safe_median_local(dd), ...
        all(ee),all(~ee),logical(P.action_ambiguous),logical(P.alpha_emitted),double(P.alpha_hat), ...
        'VariableNames',{'scenario_id','replicate','panel','p_over_n','rho','delta','contamination_POSTHOC','leverage_POSTHOC','n_anchors', ...
        'Jmax_min','Jmax_max','Jmax_median','Smed_min','Smed_max','Smed_median','DS_min','DS_max','DS_median', ...
        'all_anchors_emit','all_anchors_noemit','action_ambiguous','final_emit','alpha_hat'})]; %#ok<AGROW>
end
assert(height(Aunit)==1230,'TFKM action unit table size mismatch.');
Aq=quantile_summary_local(Aunit,{'panel','contamination_POSTHOC'}, ...
    {'Jmax_min','Jmax_max','Smed_min','Smed_max','DS_min','DS_max'});
Aauc=table();
features=["Jmax_min","Jmax_max","Smed_min","Smed_max","DS_min","DS_max"];
orientation=[1 1 -1 -1 1 1];
for j=1:numel(features)
    x=double(Aunit.(features(j)));if orientation(j)<0,x=-x;end
    auc=rank_auc_local(x,double(Aunit.leverage_POSTHOC));
    Aauc=[Aauc;table(features(j),orientation(j),auc,'VariableNames',{'feature','orientation_high_means_leverage','AUC_leverage_vs_nonleverage'})]; %#ok<AGROW>
end

% -------------------------------------------------------------------------
% A2. TRKM radial vs archived.
% -------------------------------------------------------------------------
Trad=table();
for i=1:height(PR)
    sid=double(PR.scenario_id(i));rep=double(PR.replicate(i));ic=find(double(catalog.scenario_id)==sid,1);assert(~isempty(ic));c=catalog(ic,:);
    f=fullfile(trBlind,sprintf('s%d',sid),sprintf('rep%02d',rep),'blind_final.mat');assert(isfile(f));Z=load(f,'Prediction');P=Z.Prediction;assert(height(P)==1);
    lev=string(c.contamination)=="leverage";a0=double(c.alpha0);
    re=logical(P.radial_emit);ae=logical(P.archived_emit);fe=logical(P.alpha_emitted);
    ra=double(P.radial_alpha_hat);aa=double(P.archived_alpha_hat);fa=double(P.alpha_hat);
    Trad=[Trad;table(sid,rep,string(c.panel),double(c.p_over_n),double(c.rho),double(c.delta),string(c.contamination),lev,a0, ...
        re,ae,fe,ra,aa,fa,double(P.radial_high_weight),double(P.radial_delta_BIC),double(P.radial_separation_sd),double(P.archived_delta_RV), ...
        re&&lev,ae&&lev,fe&&lev,re&&~lev,ae&&~lev,fe&&~lev, ...
        lev&&re&&isfinite(ra)&&abs(ra-a0)<=1e-12,lev&&ae&&isfinite(aa)&&abs(aa-a0)<=1e-12,lev&&fe&&isfinite(fa)&&abs(fa-a0)<=1e-12, ...
        'VariableNames',{'scenario_id','replicate','panel','p_over_n','rho','delta','contamination_POSTHOC','leverage_POSTHOC','alpha0_POSTHOC', ...
        'radial_emit','archived_emit','final_emit','radial_alpha_hat','archived_alpha_hat','final_alpha_hat','radial_high_weight','radial_delta_BIC','radial_separation_sd','archived_delta_RV', ...
        'radial_TP','archived_TP','final_TP','radial_FP','archived_FP','final_FP','radial_alpha_exact','archived_alpha_exact','final_alpha_exact'})]; %#ok<AGROW>
end
assert(height(Trad)==1230,'TRKM radial table size mismatch.');
TrSummary=trkm_rule_summary_local(Trad);
TrPattern=count_pattern_local(Trad);
TrQuant=quantile_summary_local(Trad,{'panel','contamination_POSTHOC'}, ...
    {'radial_delta_BIC','radial_high_weight','radial_separation_sd','archived_delta_RV'});
TrAuc=table();
for f=["radial_delta_BIC","radial_separation_sd","radial_high_weight","archived_delta_RV"]
    auc=rank_auc_local(double(Trad.(f)),double(Trad.leverage_POSTHOC));
    TrAuc=[TrAuc;table(f,auc,'VariableNames',{'feature','AUC_leverage_vs_nonleverage'})]; %#ok<AGROW>
end

% -------------------------------------------------------------------------
% A3. Partition-equivalence sets at alpha0 from already frozen state fits.
% -------------------------------------------------------------------------
Eq=table();EqPairs=table();
for im=1:2
    if im==1,method="TFKM";cubeBase=tfCube;else,method="TRKM";cubeBase=trCube;end
    for iu=1:height(manifest)
        sid=double(manifest.scenario_id(iu));rep=double(manifest.replicate(iu));ic=find(double(catalog.scenario_id)==sid,1);c=catalog(ic,:);cfg=table2struct(c);
        sim=generate_fixed_dataset_v1(cfg,double(manifest.data_seed(iu)));[y,out]=truth_local(sim);in=~out;
        cf=fullfile(cubeBase,sprintf('s%d',sid),sprintf('rep%02d',rep),'cube.mat');assert(isfile(cf));Q=load(cf,'Candidates','Labels');C=Q.Candidates;L=Q.Labels;
        a0=double(c.alpha0);ids=find(abs(double(C.alpha)-a0)<=1e-12);assert(numel(ids)==15,'Expected 15 alpha0 states.');
        truthidx=find(double(C.r)==double(c.r0)&double(C.g)==double(c.g0)&abs(double(C.alpha)-a0)<=1e-12,1);assert(~isempty(truthidx));
        li=truthidx;if ismember('row_index',C.Properties.VariableNames),li=double(C.row_index(truthidx));end
        labTrueState=double(L(:,li));
        rr=zeros(numel(ids),1);gg=rr;ariTruth=nan(size(rr));ariTS=nan(size(rr));
        for k=1:numel(ids)
            z=ids(k);ix=z;if ismember('row_index',C.Properties.VariableNames),ix=double(C.row_index(z));end
            lab=double(L(:,ix));rr(k)=double(C.r(z));gg(k)=double(C.g(z));
            ariTruth(k)=adjusted_rand_index_v1(y(in),lab(in));
            ariTS(k)=adjusted_rand_index_v1(labTrueState(in),lab(in));
        end
        best=max(ariTruth);E99=[rr(ariTruth>=.99-1e-12),gg(ariTruth>=.99-1e-12)];
        ETS=[rr(ariTS>=.99-1e-12),gg(ariTS>=.99-1e-12)];
        ENB=[rr(ariTruth>=best-.01-1e-12),gg(ariTruth>=best-.01-1e-12)];
        Eq=[Eq;table(method,sid,rep,string(c.panel),string(c.contamination),double(c.p_over_n),double(c.rho),double(c.delta),double(c.r0),double(c.g0),a0, ...
            best,numel(unique(E99(:,1))),numel(unique(E99(:,2))),size(E99,1),set_text_local(E99), ...
            numel(unique(ETS(:,1))),numel(unique(ETS(:,2))),size(ETS,1),set_text_local(ETS), ...
            numel(unique(ENB(:,1))),numel(unique(ENB(:,2))),size(ENB,1),set_text_local(ENB), ...
            contains_pair_local(E99,[double(c.r0) double(c.g0)]),all_true_g_local(E99,double(c.g0)),all_true_r_local(E99,double(c.r0)), ...
            'VariableNames',{'method','scenario_id','replicate','panel','contamination_POSTHOC','p_over_n','rho','delta','r0_POSTHOC','g0_POSTHOC','alpha0_POSTHOC', ...
            'best_alpha0_ARI_POSTHOC','E99_n_unique_r','E99_n_unique_g','E99_n_pairs','E99_pairs', ...
            'TRUESTATE_EQ99_n_unique_r','TRUESTATE_EQ99_n_unique_g','TRUESTATE_EQ99_n_pairs','TRUESTATE_EQ99_pairs', ...
            'NEARBEST01_n_unique_r','NEARBEST01_n_unique_g','NEARBEST01_n_pairs','NEARBEST01_pairs', ...
            'E99_contains_truth_pair','E99_all_pairs_have_true_g','E99_all_pairs_have_true_r'})]; %#ok<AGROW>
        for k=1:numel(ids)
            EqPairs=[EqPairs;table(method,sid,rep,string(c.panel),string(c.contamination),double(c.r0),double(c.g0),a0,rr(k),gg(k),ariTruth(k),ariTS(k), ...
                ariTruth(k)>=.99-1e-12,ariTS(k)>=.99-1e-12,ariTruth(k)>=best-.01-1e-12, ...
                'VariableNames',{'method','scenario_id','replicate','panel','contamination_POSTHOC','r0_POSTHOC','g0_POSTHOC','alpha0_POSTHOC','r','g','ARI_to_truth_POSTHOC','ARI_to_true_state_fit_POSTHOC','in_E99','in_TRUESTATE_EQ99','in_NEARBEST01'})]; %#ok<AGROW>
        end
        if mod((im-1)*height(manifest)+iu,200)==0,fprintf('Equivalence audit progress %d/2460\n',(im-1)*height(manifest)+iu);end
    end
end
EqSummary=equivalence_summary_local(Eq);
EqFreq=equivalence_pair_frequency_local(EqPairs);

% Save outputs.
writetable(Aanchor,fullfile(outdiag,'A1_TFKM_ACTION_ANCHOR_FEATURES.csv'));
writetable(Aunit,fullfile(outdiag,'A1_TFKM_ACTION_UNIT_FEATURES.csv'));
writetable(Aq,fullfile(outdiag,'A1_TFKM_ACTION_FEATURE_QUANTILES.csv'));
writetable(Aauc,fullfile(outdiag,'A1_TFKM_ACTION_FEATURE_AUC.csv'));
writetable(Trad,fullfile(outdiag,'A2_TRKM_RADIAL_ARCHIVED_UNIT.csv'));
writetable(TrSummary,fullfile(outdiag,'A2_TRKM_RULE_SUMMARY.csv'));
writetable(TrPattern,fullfile(outdiag,'A2_TRKM_EMISSION_PATTERNS.csv'));
writetable(TrQuant,fullfile(outdiag,'A2_TRKM_FEATURE_QUANTILES.csv'));
writetable(TrAuc,fullfile(outdiag,'A2_TRKM_FEATURE_AUC.csv'));
writetable(Eq,fullfile(outdiag,'A3_EQUIVALENCE_UNIT.csv'));
writetable(EqPairs,fullfile(outdiag,'A3_EQUIVALENCE_PAIR_LEVEL.csv'));
writetable(EqSummary,fullfile(outdiag,'A3_EQUIVALENCE_SUMMARY.csv'));
writetable(EqFreq,fullfile(outdiag,'A3_EQUIVALENCE_PAIR_FREQUENCY.csv'));
write_report_local(fullfile(outdiag,'V2_NOREFIT_DIAGNOSTIC_REPORT.txt'),Aauc,TrSummary,TrPattern,TrAuc,EqSummary,outdiag);

fprintf('\n========== A1 TFKM ACTION FEATURE AUC ==========\n');disp(Aauc);
fprintf('\n========== A2 TRKM RULE SUMMARY ================\n');disp(TrSummary);
fprintf('\n========== A2 TRKM EMISSION PATTERNS ============\n');disp(TrPattern);
fprintf('\n========== A3 EQUIVALENCE SUMMARY ===============\n');disp(EqSummary);
fprintf('\nNO-REFIT V2 DIAGNOSTICS COMPLETE.\nOutput: %s\n',outdiag);
OUT=struct('tfkm_anchor',Aanchor,'tfkm_unit',Aunit,'tfkm_auc',Aauc,'trkm_unit',Trad,'trkm_summary',TrSummary, ...
    'trkm_patterns',TrPattern,'equivalence_unit',Eq,'equivalence_pairs',EqPairs,'equivalence_summary',EqSummary,'output_root',string(outdiag));
end

% =========================================================================
function S=trkm_rule_summary_local(T)
rules=["RADIAL_ONLY";"ARCHIVED_ONLY";"FINAL_OR"];
S=table();
for i=1:numel(rules)
    r=rules(i);
    if r=="RADIAL_ONLY",emit=logical(T.radial_emit);ahat=double(T.radial_alpha_hat);
    elseif r=="ARCHIVED_ONLY",emit=logical(T.archived_emit);ahat=double(T.archived_alpha_hat);
    else,emit=logical(T.final_emit);ahat=double(T.final_alpha_hat);end
    lev=logical(T.leverage_POSTHOC);non=~lev;
    rec=sum(emit&lev)/sum(lev);fp=sum(emit&non)/sum(non);acc=(sum(emit&lev)+sum(~emit&non))/height(T);
    ae=lev&emit&isfinite(ahat);alphaExact=mean(abs(ahat(ae)-double(T.alpha0_POSTHOC(ae)))<=1e-12);
    addTP=NaN;addFP=NaN;
    if r=="FINAL_OR"
        addTP=sum(logical(T.archived_emit)&~logical(T.radial_emit)&lev);
        addFP=sum(logical(T.archived_emit)&~logical(T.radial_emit)&non);
    end
    S=[S;table(r,height(T),rec,fp,acc,alphaExact,sum(emit&lev),sum(emit&non),addTP,addFP, ...
        'VariableNames',{'rule','N','leverage_recall','nonleverage_false_positive_rate','action_accuracy','alpha_exact_given_emitted_leverage','N_true_positive','N_false_positive','archived_added_TP_beyond_radial','archived_added_FP_beyond_radial'})]; %#ok<AGROW>
end
end

function S=count_pattern_local(T)
pat=strings(height(T),1);
for i=1:height(T)
    r=logical(T.radial_emit(i));a=logical(T.archived_emit(i));
    if r&&a,pat(i)="BOTH";elseif r,pat(i)="RADIAL_ONLY";elseif a,pat(i)="ARCHIVED_ONLY";else,pat(i)="NEITHER";end
end
U=table(T.panel,T.contamination_POSTHOC,pat,'VariableNames',{'panel','contamination_POSTHOC','pattern'});
[G,p,c,q]=findgroups(U.panel,U.contamination_POSTHOC,U.pattern);N=splitapply(@numel,U.pattern,G);
S=table(p,c,q,N,'VariableNames',{'panel','contamination_POSTHOC','pattern','N'});
[G2,p2,c2]=findgroups(U.panel,U.contamination_POSTHOC);Den=splitapply(@numel,U.pattern,G2);K=table(p2,c2,Den,'VariableNames',{'panel','contamination_POSTHOC','Den'});
S=outerjoin(S,K,'Keys',{'panel','contamination_POSTHOC'},'MergeKeys',true,'Type','left');S.rate_within_group=S.N./S.Den;
S=sortrows(S,{'panel','contamination_POSTHOC','pattern'});
end

function S=equivalence_summary_local(E)
[G,m,p,c]=findgroups(E.method,E.panel,E.contamination_POSTHOC);
N=splitapply(@numel,E.scenario_id,G);
meanPairs=splitapply(@(x)mean(double(x),'omitnan'),E.E99_n_pairs,G);
meanR=splitapply(@(x)mean(double(x),'omitnan'),E.E99_n_unique_r,G);
meanG=splitapply(@(x)mean(double(x),'omitnan'),E.E99_n_unique_g,G);
hasTruth=splitapply(@(x)mean(double(x)),E.E99_contains_truth_pair,G);
allG=splitapply(@(x)mean(double(x)),E.E99_all_pairs_have_true_g,G);
allR=splitapply(@(x)mean(double(x)),E.E99_all_pairs_have_true_r,G);
tsPairs=splitapply(@(x)mean(double(x),'omitnan'),E.TRUESTATE_EQ99_n_pairs,G);
nbPairs=splitapply(@(x)mean(double(x),'omitnan'),E.NEARBEST01_n_pairs,G);
S=table(m,p,c,N,meanPairs,meanR,meanG,hasTruth,allG,allR,tsPairs,nbPairs, ...
    'VariableNames',{'method','panel','contamination_POSTHOC','N','mean_E99_n_pairs','mean_E99_unique_r','mean_E99_unique_g','E99_contains_truth_pair_rate','E99_all_true_g_rate','E99_all_true_r_rate','mean_true_state_equivalent_pairs','mean_nearbest01_pairs'});
S=sortrows(S,{'method','panel','contamination_POSTHOC'});
end

function S=equivalence_pair_frequency_local(P)
Q=P(P.in_E99,:);if isempty(Q),S=table();return,end
Q.dr=Q.r-Q.r0_POSTHOC;Q.dg=Q.g-Q.g0_POSTHOC;
[G,m,p,c,dr,dg]=findgroups(Q.method,Q.panel,Q.contamination_POSTHOC,Q.dr,Q.dg);N=splitapply(@numel,Q.scenario_id,G);
S=table(m,p,c,dr,dg,N,'VariableNames',{'method','panel','contamination_POSTHOC','dr','dg','N_E99_occurrences'});
S=sortrows(S,{'method','panel','contamination_POSTHOC','dr','dg'});
end

function S=quantile_summary_local(T,groupvars,features)
[G,ID]=findgroups(T(:,groupvars));
S=table();ng=height(ID);
for g=1:ng
    idx=G==g;row=ID(g,:);
    row.N=sum(idx);
    for f=1:numel(features)
        x=double(T.(features{f})(idx));x=x(isfinite(x));
        q=[NaN NaN NaN NaN NaN];if ~isempty(x),q=[qlinear_local(x,.05) qlinear_local(x,.25) qlinear_local(x,.5) qlinear_local(x,.75) qlinear_local(x,.95)];end
        nm=features{f};row.([nm '_q05'])=q(1);row.([nm '_q25'])=q(2);row.([nm '_median'])=q(3);row.([nm '_q75'])=q(4);row.([nm '_q95'])=q(5);
    end
    S=[S;row]; %#ok<AGROW>
end
end

function auc=rank_auc_local(x,y)
x=double(x(:));y=double(y(:));ok=isfinite(x)&isfinite(y);x=x(ok);y=y(ok);p=y==1;n=y==0;
if ~any(p)||~any(n),auc=NaN;return,end
[sv,ord]=sort(x);r=zeros(size(x));i=1;
while i<=numel(x),j=i;while j<numel(x)&&sv(j+1)==sv(i),j=j+1;end;rr=(i+j)/2;r(ord(i:j))=rr;i=j+1;end
auc=(sum(r(p))-sum(p)*(sum(p)+1)/2)/(sum(p)*sum(n));
end

function q=qlinear_local(x,p)
x=sort(double(x(:)));x=x(isfinite(x));if isempty(x),q=NaN;return,end;if numel(x)==1,q=x(1);return,end
t=1+(numel(x)-1)*p;lo=floor(t);hi=ceil(t);if lo==hi,q=x(lo);else,q=x(lo)+(t-lo)*(x(hi)-x(lo));end
end
function v=safe_min_local(x),x=double(x(:));x=x(isfinite(x));if isempty(x),v=NaN;else,v=min(x);end,end
function v=safe_max_local(x),x=double(x(:));x=x(isfinite(x));if isempty(x),v=NaN;else,v=max(x);end,end
function v=safe_median_local(x),x=double(x(:));x=x(isfinite(x));if isempty(x),v=NaN;else,v=median(x);end,end
function tf=contains_pair_local(M,q),tf=~isempty(M)&&any(all(M==q,2));end
function tf=all_true_g_local(M,g),tf=~isempty(M)&&all(M(:,2)==g);end
function tf=all_true_r_local(M,r),tf=~isempty(M)&&all(M(:,1)==r);end
function s=set_text_local(M)
if isempty(M),s="";return,end;M=unique(double(M),'rows');z=strings(size(M,1),1);for i=1:size(M,1),z(i)=sprintf('(%d,%d)',M(i,1),M(i,2));end;s=strjoin(z,';');
end
function [y,out]=truth_local(sim)
if isfield(sim,'labels'),y=double(sim.labels(:));elseif isfield(sim,'gtrue'),y=double(sim.gtrue(:));else,error('Truth labels missing.');end
if isfield(sim,'outlier'),out=logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),out=logical(sim.outlier_mask(:));else,error('Outlier mask missing.');end
end
function lock_vendor_local(root)
v=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));assert(isfolder(v),'Canonical vendor missing.');addpath(genpath(char(v)),'-begin');rehash path;
for q=["generate_fixed_dataset_v1","adjusted_rand_index_v1"],p=string(which(char(q)));assert(strlength(p)>0,'Missing %s.',q);assert(startsWith(p,v,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,p);end
end
function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;root="";outroot="";
for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break,end,end
assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
function write_report_local(f,Aauc,Ts,Tp,Ta,Eq,root)
fid=fopen(f,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'FRESH V2 NO-REFIT DIAGNOSTICS V1\n');fprintf(fid,'=================================\n');
fprintf(fid,'Truth open. V1 frozen. No fits/refits. Development diagnostics only.\nOutput root: %s\n\n',root);
fprintf(fid,'A1 TFKM ACTION FEATURE AUC\n%s\n',evalc('disp(Aauc)'));
fprintf(fid,'A2 TRKM RULE SUMMARY\n%s\n',evalc('disp(Ts)'));
fprintf(fid,'A2 TRKM EMISSION PATTERNS\n%s\n',evalc('disp(Tp)'));
fprintf(fid,'A2 TRKM FEATURE AUC\n%s\n',evalc('disp(Ta)'));
fprintf(fid,'A3 EQUIVALENCE SUMMARY\n%s\n',evalc('disp(Eq)'));
end
