function OUT = RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1(documents_root)
%RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1
% First truth-open evaluation of the fully frozen fresh generalization study.
%
% IMPORTANT:
%   - NO selector fitting or refitting occurs here.
%   - NO threshold/rule is estimated here.
%   - The pre-truth snapshot is verified against the active V2 frozen seal
%     BEFORE truth is regenerated.
%   - Truth is used only for evaluation after all 2460 predictions were frozen.
%
% Outputs are written to:
%   TFKM_TRKM_FRESH_GENERALIZATION_V1_output/fresh_truth_open_v1/
%
% Usage:
%   RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1;
%   RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1("C:\Users\name\Documents");
%   RUN_TFKM_TRKM_FRESH_TRUTH_OPEN_V1("SELFTEST");

if nargin>=1 && strcmpi(string(documents_root),"SELFTEST")
    assert(action_truth_local("clean")==false);
    assert(action_truth_local("bridge")==false);
    assert(action_truth_local("leverage")==true);
    assert(abs(f1_local(1,1,1)-1)<=eps);
    assert(abs(f1_local(0,1,1)-0)<=eps);
    assert(abs(f1_local(1,2,1)-(2/3))<=10*eps);
    assert(isnan(f1_local(0,0,0)));
    fprintf('SELFTEST PASS: fresh truth-open scoring contracts.\n');
    OUT=struct('selftest_pass',true);return
end
if nargin<1, documents_root=""; end
[root,outroot,S]=load_study_local(documents_root);
catalog=S.catalog;manifest=S.manifest;CFG=S.CFG;
assert(height(catalog)==41 && height(manifest)==1230,'Fresh study dimensions drifted.');

% -------------------------------------------------------------------------
% PRE-TRUTH IMMUTABILITY GATE
% -------------------------------------------------------------------------
activeSeal=fullfile(outroot,'blind_finalization_v2','FROZEN');
snapSeal=fullfile(outroot,'FROZEN_PRETRUTH_20261010');
assert(isfolder(activeSeal),'Active V2 FROZEN seal missing.');
assert(isfolder(snapSeal),'Pre-truth snapshot missing: %s',snapSeal);
files=["TFKM_HIER4_BLIND_FINAL_PREDICTIONS_FROZEN.csv", ...
       "TRKM_BLIND_FINAL_PREDICTIONS_FROZEN.csv"];
for nm=files
    a=fullfile(activeSeal,char(nm));b=fullfile(snapSeal,char(nm));
    assert(isfile(a)&&isfile(b),'Frozen prediction file missing: %s',nm);
    A=readtable(a,'TextType','string');B=readtable(b,'TextType','string');
    A=sortrows(A,{'scenario_id','replicate'});B=sortrows(B,{'scenario_id','replicate'});
    assert(isequaln(A,B),'PRETRUTH SNAPSHOT MISMATCH: %s',nm);
end
PT=readtable(fullfile(snapSeal,char(files(1))),'TextType','string');
PR=readtable(fullfile(snapSeal,char(files(2))),'TextType','string');
PT=sortrows(PT,{'scenario_id','replicate'});PR=sortrows(PR,{'scenario_id','replicate'});
assert(height(PT)==1230&&height(PR)==1230,'Frozen prediction count mismatch.');
assert(isequal(double(PT.scenario_id),double(PR.scenario_id)) && ...
       isequal(double(PT.replicate),double(PR.replicate)), ...
       'TFKM/TRKM frozen keys do not align.');

% Truth-generation functions only; no estimator routine is called below.
lock_truth_vendor_local(root);
assert(exist('generate_fixed_dataset_v1','file')==2,'generate_fixed_dataset_v1 missing.');
assert(exist('adjusted_rand_index_v1','file')==2,'adjusted_rand_index_v1 missing.');

truthroot=fullfile(outroot,'fresh_truth_open_v1');
if ~isfolder(truthroot),mkdir(truthroot);end
marker=fullfile(truthroot,'TRUTH_OPENED_AT.txt');
if ~isfile(marker)
    fid=fopen(marker,'w');assert(fid>0);
    fprintf(fid,'TFKM_TRKM_FRESH_GENERALIZATION_V1 - TRUTH OPENED\n');
    fprintf(fid,'opened_at=%s\n',char(datetime('now')));
    fprintf(fid,'pretruth_snapshot=%s\n',snapSeal);
    fprintf(fid,'TFKM_predictions=1230\nTRKM_predictions=1230\n');
    fprintf(fid,'RULE FREEZE: no selector/action/alpha/anchor/backward threshold may be changed after this point.\n');
    fclose(fid);
end
fprintf('\n============================================================\n');
fprintf('FRESH GENERALIZATION V1 - TRUTH OPENS NOW\n');
fprintf('Frozen predictions: TFKM 1230 | TRKM 1230\n');
fprintf('NO REFITS / NO RETUNING / EVALUATION ONLY\n');
fprintf('============================================================\n');

R=table();
for imethod=1:2
    if imethod==1
        P=PT;method="TFKM";
        blindBase=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_HIER4_RERUN_V1_output', ...
            'blind_finalization_hier4_v1','TFKM');
    else
        P=PR;method="TRKM";
        blindBase=fullfile(outroot,'blind_finalization_v1','TRKM');
    end
    assert(isfolder(blindBase),'Blind per-unit directory missing for %s: %s',method,blindBase);
    for i=1:height(P)
        sid=double(P.scenario_id(i));rep=double(P.replicate(i));
        ic=find(double(catalog.scenario_id)==sid,1);assert(~isempty(ic),'Catalog sid missing.');
        im=find(double(manifest.scenario_id)==sid & double(manifest.replicate)==rep,1);assert(~isempty(im),'Manifest key missing.');
        c=catalog(ic,:);cfg=table2struct(c);
        sim=generate_fixed_dataset_v1(cfg,double(manifest.data_seed(im)));
        [ytrue,out]=truth_local(sim);in=~out;
        n=numel(ytrue);assert(n==double(c.n),'Truth dataset n drift.');

        cont=string(c.contamination);isLev=action_truth_local(cont);
        emitted=logical(P.alpha_emitted(i));
        actionCorrect=(emitted==isLev);
        alphaHat=double(P.alpha_hat(i));alpha0=double(c.alpha0);
        alphaExact=false;alphaAbsErr=NaN;
        if isLev && emitted && isfinite(alphaHat)
            alphaAbsErr=abs(alphaHat-alpha0);
            alphaExact=alphaAbsErr<=1e-12;
        end
        terminalUnique=logical(P.terminal_unique(i));
        tr=double(P.terminal_r(i));tg=double(P.terminal_g(i));
        structureExact=terminalUnique && isfinite(tr) && isfinite(tg) && ...
            tr==double(c.r0) && tg==double(c.g0);
        actionAmb=false;if ismember('action_ambiguous',P.Properties.VariableNames),actionAmb=logical(P.action_ambiguous(i));end
        fitAmb=false;if ismember('final_fit_ambiguous',P.Properties.VariableNames),fitAmb=logical(P.final_fit_ambiguous(i));end
        ambiguous=(~terminalUnique)||actionAmb||fitAmb;
        tripleExact=isLev && structureExact && actionCorrect && alphaExact;
        jointCorrect=structureExact && actionCorrect && (~isLev || alphaExact);

        ari=NaN;trimPrec=NaN;trimRec=NaN;trimF1=NaN;trimFPR=NaN;trimRate=NaN;fitMaterialized=false;
        bf=fullfile(blindBase,sprintf('s%d',sid),sprintf('rep%02d',rep),'blind_final.mat');
        assert(isfile(bf),'Missing frozen blind_final.mat: %s',bf);
        Z=load(bf,'FinalLabels','FinalTrim','Prediction');
        % Firewall reproduction: prediction saved in MAT must agree on key.
        assert(height(Z.Prediction)==1 && double(Z.Prediction.scenario_id)==sid && double(Z.Prediction.replicate)==rep, ...
            'Per-unit prediction key mismatch.');
        if ~isempty(Z.FinalLabels)
            lab=double(Z.FinalLabels(:));assert(numel(lab)==n,'FinalLabels length mismatch.');
            fitMaterialized=true;
            if sum(in)>1,ari=adjusted_rand_index_v1(ytrue(in),lab(in));end
            if ~isempty(Z.FinalTrim)
                trim=logical(Z.FinalTrim(:));assert(numel(trim)==n,'FinalTrim length mismatch.');
                trimRate=mean(trim);
                trimFPR=sum(trim&in)/max(1,sum(in));
                if any(out)
                    tp=sum(trim&out);pred=sum(trim);truth=sum(out);
                    if pred>0,trimPrec=tp/pred;end
                    trimRec=tp/truth;
                    trimF1=f1_local(tp,pred,truth);
                end
            end
        end

        row=table(method,sid,rep,string(c.panel),double(c.n),double(c.p),double(c.p_over_n), ...
            double(c.r0),double(c.g0),double(c.delta),double(c.rho),cont,alpha0, ...
            isLev,emitted,actionCorrect,alphaHat,alphaExact,alphaAbsErr, ...
            terminalUnique,tr,tg,structureExact,ambiguous,tripleExact,jointCorrect, ...
            fitMaterialized,ari,trimRate,trimPrec,trimRec,trimF1,trimFPR, ...
            'VariableNames',{'method','scenario_id','replicate','panel','n','p','p_over_n','r0_POSTHOC','g0_POSTHOC', ...
            'delta_POSTHOC','rho_POSTHOC','contamination_POSTHOC','alpha0_POSTHOC','leverage_action_truth_POSTHOC', ...
            'alpha_emitted','action_correct_POSTHOC','alpha_hat','alpha_exact_POSTHOC','alpha_abs_error_POSTHOC', ...
            'terminal_unique','terminal_r','terminal_g','structure_exact_POSTHOC','ambiguous_output_POSTHOC', ...
            'triple_exact_leverage_POSTHOC','joint_correct_POSTHOC','final_fit_materialized', ...
            'final_inlier_ARI_POSTHOC','final_trim_rate_POSTHOC','final_trim_precision_POSTHOC', ...
            'final_trim_recall_POSTHOC','final_trim_F1_POSTHOC','final_trim_inlier_FPR_POSTHOC'});
        R=[R;row]; %#ok<AGROW>
    end
end
R=sortrows(R,{'method','scenario_id','replicate'});
assert(height(R)==2460,'Unit truth table must have 2460 rows.');
writetable(R,fullfile(truthroot,'FRESH_TRUTH_UNIT_RESULTS.csv'));

Overall=summarize_groups_local(R,{'method'});
ByCont=summarize_groups_local(R,{'method','contamination_POSTHOC'});
ByPanel=summarize_groups_local(R,{'method','panel'});
ByCell=summarize_groups_local(R,{'method','scenario_id','panel','p_over_n','rho_POSTHOC','delta_POSTHOC','r0_POSTHOC','g0_POSTHOC','contamination_POSTHOC'});
writetable(Overall,fullfile(truthroot,'FRESH_TRUTH_SUMMARY_OVERALL.csv'));
writetable(ByCont,fullfile(truthroot,'FRESH_TRUTH_SUMMARY_BY_CONTAMINATION.csv'));
writetable(ByPanel,fullfile(truthroot,'FRESH_TRUTH_SUMMARY_BY_PANEL.csv'));
writetable(ByCell,fullfile(truthroot,'FRESH_TRUTH_SUMMARY_BY_CELL.csv'));
write_text_report_local(fullfile(truthroot,'FRESH_TRUTH_REPORT.txt'),Overall,ByCont,ByPanel,truthroot);

fprintf('\n================ OVERALL ================\n');disp(Overall);
fprintf('\n============= BY CONTAMINATION ===========\n');disp(ByCont);
fprintf('\n================ BY PANEL ================\n');disp(ByPanel);
fprintf('\nTRUTH-OPEN EVALUATION COMPLETE.\nOutput: %s\n',truthroot);
fprintf('The frozen procedures must not be retuned after this output.\n');
OUT=struct('unit_results',R,'overall',Overall,'by_contamination',ByCont,'by_panel',ByPanel,'by_cell',ByCell,'output_root',string(truthroot));
end

function S=summarize_groups_local(R,vars)
G=findgroups(R(:,vars));ng=max(G);S=table();
for g=1:ng
    U=R(G==g,:);head=U(1,vars);
    lev=logical(U.leverage_action_truth_POSTHOC);
    emit=logical(U.alpha_emitted);fit=logical(U.final_fit_materialized);
    nLev=sum(lev);nNon=sum(~lev);nEmitLev=sum(lev&emit);
    vals={height(U),mean(double(U.structure_exact_POSTHOC)),mean(double(U.action_correct_POSTHOC)), ...
        safe_mean_local(double(~emit(~lev))),safe_mean_local(double(emit(lev))), ...
        safe_mean_local(double(U.alpha_exact_POSTHOC(lev&emit))),safe_mean_local(double(U.alpha_exact_POSTHOC(lev))), ...
        safe_mean_local(double(U.triple_exact_leverage_POSTHOC(lev))),mean(double(U.joint_correct_POSTHOC)), ...
        mean(double(U.ambiguous_output_POSTHOC)),mean(double(fit)), ...
        safe_mean_local(double(U.final_inlier_ARI_POSTHOC(fit))), ...
        safe_mean_local(double(U.final_trim_F1_POSTHOC(fit&lev))), ...
        nLev,nNon,nEmitLev};
    tail=cell2table(vals,'VariableNames',{'N','structure_exact_rate','action_accuracy', ...
        'nonleverage_abstention_rate','leverage_action_recall','alpha_exact_given_leverage_emit', ...
        'alpha_exact_unconditional_leverage','triple_exact_rate_leverage','joint_correct_rate', ...
        'ambiguous_output_rate','final_fit_materialized_rate','mean_final_inlier_ARI_when_materialized', ...
        'mean_trim_F1_leverage_when_materialized','N_leverage','N_nonleverage','N_emitted_leverage'});
    S=[S;[head,tail]]; %#ok<AGROW>
end
end

function m=safe_mean_local(x)
x=x(isfinite(x));if isempty(x),m=NaN;else,m=mean(x);end
end

function f=f1_local(tp,pred,truth)
% Standard binary F1 = 2TP / (predicted positives + true positives).
den=pred+truth;
if den<=0,f=NaN;else,f=2*tp/den;end
end

function tf=action_truth_local(cont)
tf=string(cont)=="leverage";
end

function [y,out]=truth_local(sim)
if isfield(sim,'labels'),y=double(sim.labels(:));elseif isfield(sim,'gtrue'),y=double(sim.gtrue(:));else,error('Truth labels missing.');end
if isfield(sim,'outlier'),out=logical(sim.outlier(:));elseif isfield(sim,'outlier_mask'),out=logical(sim.outlier_mask(:));else,error('Outlier mask missing.');end
assert(numel(y)==numel(out),'Truth label/outlier length mismatch.');
end

function write_text_report_local(f,O,C,P,root)
fid=fopen(f,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'TFKM/TRKM FRESH GENERALIZATION V1 - FIRST TRUTH-OPEN REPORT\n');
fprintf(fid,'============================================================\n');
fprintf(fid,'Evaluation only. No refits and no truth-driven retuning.\n');
fprintf(fid,'Output root: %s\n\n',root);
fprintf(fid,'OVERALL\n');fprintf(fid,'%s\n',evalc('disp(O)'));
fprintf(fid,'BY CONTAMINATION\n');fprintf(fid,'%s\n',evalc('disp(C)'));
fprintf(fid,'BY PANEL\n');fprintf(fid,'%s\n',evalc('disp(P)'));
fprintf(fid,'\nMetric note: correct action means EMIT for leverage and NO ACTION for clean/bridge.\n');
fprintf(fid,'Alpha exactness is evaluated only for leverage, because bridge is a correct-abstention regime in the frozen procedure.\n');
fprintf(fid,'Final ARI/trimming are reported only where the blind finalizer materialized a final fit. No post-truth fit is created.\n');
end

function lock_truth_vendor_local(root)
tf=string(fullfile(root,'RSS_V2_CANONICAL_REAPPLY_v1_1','vendor'));
assert(isfolder(tf),'Canonical truth-generation vendor missing: %s',tf);
addpath(genpath(char(tf)),'-begin');rehash path;
for q=["generate_fixed_dataset_v1","adjusted_rand_index_v1"]
    p=string(which(char(q)));assert(strlength(p)>0,'Missing %s.',q);
    assert(startsWith(p,tf,'IgnoreCase',true),'PATH LOCK FAIL %s -> %s',q,p);
end
end

function [root,outroot,S]=load_study_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";outroot="";
for c=cands(:)'
    o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
    if isfile(fullfile(o,'study.mat')),root=c;outroot=string(o);break;end
end
assert(strlength(root)>0,'Fresh study not found.');S=load(fullfile(outroot,'study.mat'));
end
