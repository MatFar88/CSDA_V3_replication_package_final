function sel=select_alpha_performance_frozen_v1(P,method)
%SELECT_ALPHA_PERFORMANCE_FROZEN_V1 Frozen pre-fresh general-alpha selectors.
% Truth-blind: uses only the candidate path and method identity.
method=string(method);
if method=="REG_TFKM_HD"
    R=table("GAIN_OBJ_L0800","GAIN_OBJ",0.8,NaN, ...
        'VariableNames',{'rule','family','p1','p2'});
elseif method=="REG_TRKM_HD"
    R=table("STOP_OWN_T010","STOP_OWN",0.10,NaN, ...
        'VariableNames',{'rule','family','p1','p2'});
else
    error('generalAlphaFresh:Method','Unsupported method %s.',method);
end
sel=apply_general_alpha_rule_v2(P,R);
end
