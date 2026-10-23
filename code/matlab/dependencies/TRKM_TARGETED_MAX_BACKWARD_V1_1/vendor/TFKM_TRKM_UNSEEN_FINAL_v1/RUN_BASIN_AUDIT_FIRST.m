% RUN_BASIN_AUDIT_FIRST
% Exact deterministic multistart-basin audit for the frozen v11 study.
%
% 1) Set your EXISTING v11 output root (the folder containing study.mat,
%    results/, summary/, dgp_certification/).
%
% output_root = "C:\\...\\your_v11_output_root";
%
% 2) Start with the decisive clean/bridge audit of the CURRENT selected
%    fixed triplets for both methods (32 pools x 100 starts, resumable):
%
% report = run_multistart_basin_audit_v1(output_root,"PRIMARY");
%
% Faster TFKM-only first pass (16 pools):
% report = run_multistart_basin_audit_v1(output_root,"TFKM_SELECTED");
%
% Then, if desired, add TFKM true-fixed controls without rerunning completed
% selected pools:
% report = run_multistart_basin_audit_v1(output_root,"TFKM_SELECTED_TRUE");
%
% Full clean/bridge audit (TFKM selected + true-fixed; TRKM selected):
% report = run_multistart_basin_audit_v1(output_root,"FULL_CLEAN_BRIDGE");
%
% Results are written under:
%   <output_root>\multistart_basin_audit_v1\summary\
%
% The frozen estimator is never changed. Every pool is first checked for
% exact reproduction of the stored v11 current winner. Truth is opened only
% for post-fit diagnostics.

pkg=fileparts(mfilename('fullpath'));
addpath(genpath(pkg),'-begin');
fprintf('\nBasin-audit functions loaded from:\n%s\n\n',pkg);
fprintf('Recommended first command:\n');
fprintf('report = run_multistart_basin_audit_v1(output_root,"TFKM_SELECTED");\n');
