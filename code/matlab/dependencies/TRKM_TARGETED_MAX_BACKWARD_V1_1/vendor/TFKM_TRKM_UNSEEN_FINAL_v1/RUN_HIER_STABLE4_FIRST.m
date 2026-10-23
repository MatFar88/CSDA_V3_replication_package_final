pkg=fileparts(mfilename('fullpath')); addpath(genpath(pkg),'-begin');

%% TFKM HIER_STABLE4 full-pipeline development challenge
% 1) Point to the ORIGINAL frozen v11 output_root.
output_root = "C:\\PATH\\TO\\TFKM_TRKM_FINAL_FRESH_VALIDATION_v1_output";

%% 2) Recommended first gate: clean + bridge only.
% One MATLAB session:
% run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",1,1);
%
% Or four parallel MATLAB sessions, one command per session:
% run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",1,4);
% run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",2,4);
% run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",3,4);
% run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",4,4);

%% 3) Check progress.
% inspect_hier_stable4_progress_v1(output_root,"CLEAN_BRIDGE");

%% 4) When all 16 clean/bridge jobs are complete, apply the UNCHANGED full pipeline.
% report = merge_hier_stable4_full_pipeline_v1(output_root,"CLEAN_BRIDGE");

%% 5) Only if the clean/bridge development result is coherent, complete leverage too.
% run_hier_stable4_queue_v1(output_root,"ALL_TFKM",1,4);
% run_hier_stable4_queue_v1(output_root,"ALL_TFKM",2,4);
% run_hier_stable4_queue_v1(output_root,"ALL_TFKM",3,4);
% run_hier_stable4_queue_v1(output_root,"ALL_TFKM",4,4);
% report_all = merge_hier_stable4_full_pipeline_v1(output_root,"ALL_TFKM");
