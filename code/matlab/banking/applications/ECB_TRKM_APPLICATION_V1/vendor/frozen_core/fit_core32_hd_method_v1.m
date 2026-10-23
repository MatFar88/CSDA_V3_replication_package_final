function fit = fit_core32_hd_method_v1(X,cfg,method,options,seed)
%FIT_CORE32_HD_METHOD_V1 Dispatcher for final CORE32 fixed validation.
method=string(method);
if method=="REG_TFKM_HD"
    pair_seed=seed+110000;
elseif method=="REG_TRKM_HD"
    pair_seed=seed+220000;
else
    error('core32:Method','Unsupported method %s.',method);
end
fit=run_hd_multistart_top_v1(X,cfg.r0,cfg.c0,cfg.alpha0,method,options,pair_seed);
if fit.valid
    try
        fit=normalize_fit_contract_v1(fit,size(X,1),cfg.c0,cfg.alpha0,method);
    catch ME
        fit.valid=false;
        fit.message=string(ME.identifier)+": "+string(ME.message);
    end
end
end
