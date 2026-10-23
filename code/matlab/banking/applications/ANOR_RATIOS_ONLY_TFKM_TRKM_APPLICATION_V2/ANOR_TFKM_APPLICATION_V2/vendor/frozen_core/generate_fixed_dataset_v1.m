function sim = generate_fixed_dataset_v1(cfg,seed)
% Adapter to the validated contaminated latent generator.
c=struct('n',cfg.n,'p',cfg.p,'r',cfg.r0,'c',cfg.c0,'r_true',cfg.r0,'c_true',cfg.c0,'delta',cfg.delta,'rho',cfg.rho,'separation',cfg.delta,'spectral_ratio',cfg.rho,'contamination',cfg.contamination,'alpha',cfg.alpha0);
sim=generate_contaminated_latent_v1(c,seed);
% Normalize contract used by all modules.
if isfield(sim,'gtrue'), sim.labels=sim.gtrue; end
if isfield(sim,'outlier_mask'), sim.outlier=logical(sim.outlier_mask); end
if ~isfield(sim,'outlier'), sim.outlier=false(cfg.n,1); end
end
