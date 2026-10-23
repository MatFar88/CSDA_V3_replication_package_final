function run_hier_stable4_queue_v1(output_root,scope,worker_id,n_workers)
%RUN_HIER_STABLE4_QUEUE_V1 Run/resume a deterministic subset of challenger jobs.
%
% Examples:
%   run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",1,1)
%   % Four parallel MATLAB sessions:
%   run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",1,4)
%   run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",2,4)
%   run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",3,4)
%   run_hier_stable4_queue_v1(output_root,"CLEAN_BRIDGE",4,4)
%
% Scope ALL_TFKM runs all 4 scenarios x 8 reps after the clean/bridge gate.

if nargin<2 || strlength(string(scope))==0,scope="CLEAN_BRIDGE";end
if nargin<3 || isempty(worker_id),worker_id=1;end
if nargin<4 || isempty(n_workers),n_workers=1;end
scope=upper(string(scope)); worker_id=double(worker_id); n_workers=double(n_workers);
assert(worker_id>=1 && worker_id<=n_workers && mod(worker_id,1)==0 && mod(n_workers,1)==0,'Invalid worker assignment.');

switch scope
    case "CLEAN_BRIDGE"
        scenarios=[23001 23004];
    case "ALL_TFKM"
        scenarios=[23001 23004 23005 23007];
    otherwise
        error('hierStable4:Scope','Unknown scope %s.',scope);
end
jobs=[];
for s=scenarios
    for rep=1:8
        jobs=[jobs; s rep]; %#ok<AGROW>
    end
end
assigned=find(mod((1:size(jobs,1))-1,n_workers)==worker_id-1);
fprintf('\nHIER_STABLE4 QUEUE | scope=%s | worker=%d/%d | jobs=%d\n',scope,worker_id,n_workers,numel(assigned));
for k=reshape(assigned,1,[])
    run_hier_stable4_full_pipeline_job_v1(output_root,jobs(k,1),jobs(k,2));
end
end
