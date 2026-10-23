function E=evaluate_unseen_truth_v1(output_root)
%EVALUATE_UNSEEN_TRUTH_V1 Open sealed truth ONLY after selection freeze exists.
root=char(string(output_root)); tol=1e-12;
freeze=fullfile(root,'operational_selection','SELECTION_FROZEN.mat');
assert(isfile(freeze),'Operational selection is not frozen. Run select_unseen_operational_v1 first.');
F=load(freeze,'Dec','selection_frozen_at','rule_version'); Dec=F.Dec;
S=load(fullfile(root,'study.mat'),'catalog','manifest'); catalog=S.catalog; manifest=S.manifest;

Z=table();
for i=1:height(manifest)
    ses=string(manifest.session(i));
    f=fullfile(root,'_SEALED_TRUTH_DO_NOT_OPEN_BEFORE_SELECTION',sprintf('sealed_truth_%s.csv',char(ses)));
    assert(isfile(f),'Missing sealed truth file: %s',f);
    Q=readtable(f,'TextType','string');
    assert(height(Q)==double(manifest.expected_total(i)),'Sealed session %s incomplete.',ses);
    Z=[Z;Q]; %#ok<AGROW>
end
Z=sortrows(Z,{'method','scenario_id','replicate','r','g','alpha'});
assert(height(Z)==9600,'Expected 9600 sealed truth candidate rows.');

D=Dec;
D.r0=nan(height(D),1); D.g0=nan(height(D),1); D.alpha0=nan(height(D),1);
D.contamination=strings(height(D),1); D.truth_regime=strings(height(D),1);
D.exact_rg=false(height(D),1); D.regime_correct=false(height(D),1);
D.alpha_exact=nan(height(D),1); D.joint_success=false(height(D),1);
D.ari_inlier=nan(height(D),1); D.trim_f1=nan(height(D),1); D.trim_precision=nan(height(D),1); D.trim_recall=nan(height(D),1);

for i=1:height(D)
    sid=double(D.scenario_id(i)); rep=double(D.replicate(i)); method=string(D.method(i));
    c=catalog(double(catalog.scenario_id)==sid,:); assert(height(c)==1,'Catalog row missing.');
    r0=double(c.r0); g0=double(c.c0); a0=double(c.alpha0); cont=string(c.contamination);
    if a0<=tol, truthReg="CLEAN"; elseif contains(lower(cont),"bridge"), truthReg="BRIDGE_AMBIGUOUS"; else, truthReg="LEVERAGE"; end
    D.r0(i)=r0; D.g0(i)=g0; D.alpha0(i)=a0; D.contamination(i)=cont; D.truth_regime(i)=truthReg;
    D.exact_rg(i)=double(D.selected_r(i))==r0 && double(D.selected_g(i))==g0;
    D.regime_correct(i)=string(D.predicted_regime(i))==truthReg;
    if logical(D.alpha_identified(i))
        D.alpha_exact(i)=double(abs(double(D.alpha_hat(i))-a0)<=tol);
        q=Z(double(Z.scenario_id)==sid & double(Z.replicate)==rep & string(Z.method)==method & ...
            double(Z.r)==double(D.selected_r(i)) & double(Z.g)==double(D.selected_g(i)) & ...
            abs(double(Z.alpha)-double(D.alpha_hat(i)))<=tol,:);
        assert(height(q)==1,'Final candidate truth row not unique sid=%d rep=%d %s.',sid,rep,method);
        D.ari_inlier(i)=double(q.ari_inlier); D.trim_f1(i)=double(q.trim_f1);
        D.trim_precision(i)=double(q.trim_precision); D.trim_recall(i)=double(q.trim_recall);
    end
    if D.regime_correct(i)
        if truthReg=="BRIDGE_AMBIGUOUS"
            D.joint_success(i)=~logical(D.alpha_identified(i));
        else
            D.joint_success(i)=logical(D.alpha_identified(i)) && D.alpha_exact(i)==1;
        end
    end
end

Summary=table(); methods=unique(string(D.method),'stable');
for m=methods'
    Qm=D(string(D.method)==m,:);
    sids=unique(double(Qm.scenario_id),'stable');
    for sid=sids'
        Q=Qm(double(Qm.scenario_id)==sid,:);
        Summary=[Summary;summary_row(m,sid,string(Q.contamination(1)),Q)]; %#ok<AGROW>
    end
    Summary=[Summary;summary_row(m,0,"ALL",Qm)]; %#ok<AGROW>
end

outdir=fullfile(root,'truth_evaluation');
writetable(D,fullfile(outdir,'unseen_truth_evaluation_detail.csv'));
writetable(Summary,fullfile(outdir,'unseen_truth_evaluation_summary.csv'));
truth_opened_at=string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')); %#ok<NASGU>
save(fullfile(outdir,'UNSEEN_TRUTH_EVALUATION.mat'),'D','Summary','truth_opened_at','-v7.3');

fprintf('\n============================================================\n');
fprintf('UNSEEN TRUTH OPENED AFTER FROZEN SELECTION\n');
fprintf('Selection frozen at: %s\n',F.selection_frozen_at);
fprintf('Rule version: %s\n',F.rule_version);
fprintf('============================================================\n\n');
disp(Summary)

E=struct('detail',D,'summary',Summary,'selection_frozen_at',F.selection_frozen_at, ...
    'truth_opened_at',truth_opened_at,'rule_version',F.rule_version);
end

function row=summary_row(method,sid,label,Q)
id=logical(Q.alpha_identified);
if any(id), aexact=mean(Q.alpha_exact(id),'omitnan'); else, aexact=NaN; end
row=table(method,sid,label,height(Q),mean(Q.exact_rg),mean(Q.regime_correct), ...
    mean(Q.alpha_identified),aexact,mean(Q.joint_success), ...
    mean(Q.ari_inlier,'omitnan'),mean(Q.trim_f1,'omitnan'), ...
    'VariableNames',{'method','scenario_id','contamination','N','exact_rg_rate', ...
    'regime_accuracy','alpha_identification_rate','alpha_exact_conditional_identified', ...
    'joint_success_rate','mean_ari_when_alpha_identified','mean_trim_f1_when_alpha_identified'});
end
