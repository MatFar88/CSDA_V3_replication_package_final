function T=inspect_final_fresh_progress_v1(output_root)
%INSPECT_FINAL_FRESH_PROGRESS_V1 Lightweight progress table for A:H.
output_root=char(string(output_root));S=load(fullfile(output_root,'study.mat'));M=S.manifest;T=M;T.completed_rows=zeros(height(M),1);T.complete=false(height(M),1);
for i=1:height(M)
 f=fullfile(output_root,'results',sprintf('joint_candidates_%s.csv',char(M.session(i))));cp=fullfile(output_root,'results',sprintf('joint_candidates_%s_checkpoint.mat',char(M.session(i))));
 if isfile(f),Q=readtable(f);T.completed_rows(i)=height(Q);elseif isfile(cp),Q=load(cp,'R_checkpoint');if isfield(Q,'R_checkpoint'),T.completed_rows(i)=height(Q.R_checkpoint);end,end
 T.complete(i)=T.completed_rows(i)>=double(T.expected_candidate_fits(i));
end
disp(T);fprintf('Total completed candidate rows: %d / %d\n',sum(T.completed_rows),sum(double(T.expected_candidate_fits)));
end
