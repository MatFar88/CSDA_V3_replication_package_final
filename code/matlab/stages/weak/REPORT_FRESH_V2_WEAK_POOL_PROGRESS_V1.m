function REPORT_FRESH_V2_WEAK_POOL_PROGRESS_V1(documents_root)
% Read-only progress report for the 300 weak true-state method-unit pools.
if nargin<1,documents_root="";end
if strlength(string(documents_root))>0,cands=string(documents_root);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end;outroot="";
for c=cands(:)',o=fullfile(c,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');if isfile(fullfile(o,'study.mat')),outroot=string(o);break,end,end
assert(strlength(outroot)>0,'Fresh study not found.');d=fullfile(outroot,'fresh_v2_diagnostic_suite_v1','weak_true_state_pool_audit_v1');
if ~isfolder(d),fprintf('WEAK POOL AUDIT: 0/300 complete.\n');return,end
F=dir(fullfile(d,'worker_*_of_*_pool.csv'));K=table();
for i=1:numel(F),Q=readtable(fullfile(F(i).folder,F(i).name),'TextType','string');if ~isempty(Q),K=[K;Q(:,{'method','scenario_id','replicate'})];end,end %#ok<AGROW>
if isempty(K),n=0;nt=0;nr=0;else,K=unique(K,'rows');n=height(K);nt=sum(string(K.method)=="TFKM");nr=sum(string(K.method)=="TRKM");end
fprintf('\n============================================================\n');fprintf('WEAK TRUE-STATE POOL AUDIT PROGRESS\n');fprintf('============================================================\n');fprintf('ALL : %d/300 (%.1f%%)\n',n,100*n/300);fprintf('TFKM: %d/150\n',nt);fprintf('TRKM: %d/150\n',nr);
end
