pkg=fileparts(mfilename('fullpath')); addpath(genpath(pkg),'-begin');

%% ================================================================
% TFKM/TRKM FINAL UNSEEN CONFIRMATION v1
% ONE NEW ROOT FOR BOTH METHODS.
% DO NOT point this at any development output folder.
% ================================================================

output_root = "C:\Users\matteo.farne2\Documents\TFKM_TRKM_UNSEEN_FINAL_v1_output";

%% 1. Create/freeze the unseen design ONCE.
% If study.mat already exists, DO NOT rerun setup with a different root/seed.
if ~isfile(fullfile(output_root,'study.mat'))
    study = setup_unseen_validation_v1(output_root);
else
    fprintf('Existing unseen study found: %s\n',fullfile(output_root,'study.mat'));
end

%% 2. RUN A:H. Recommended: parallel MATLAB sessions.
% One session per replicate; each runs BOTH methods.
%
% Session 1:
% run_unseen_session_v1("A",output_root);
% Session 2:
% run_unseen_session_v1("B",output_root);
% ...
% Session 8:
% run_unseen_session_v1("H",output_root);
%
% Any session is resumable by rerunning the same command.

%% 3. Progress at any time:
% inspect_unseen_progress_v1(output_root)

%% 4. ONLY after A:H are all complete:
% M = merge_unseen_operational_v1(output_root);
%
% Freeze all operational metaparameter decisions, STILL WITHOUT TRUTH:
% R_UNSEEN = select_unseen_operational_v1(output_root);
%
% This creates operational_selection\SELECTION_FROZEN.mat.

%% 5. ONLY AFTER SELECTION_FROZEN.mat EXISTS:
% E_UNSEEN = evaluate_unseen_truth_v1(output_root);
% E_UNSEEN.summary
%
% Until step 4 is complete, DO NOT inspect:
%   _SEALED_TRUTH_DO_NOT_OPEN_BEFORE_SELECTION
