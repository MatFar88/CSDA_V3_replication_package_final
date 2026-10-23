function P=inspect_hier_stable4_progress_v1(output_root,scope)
%INSPECT_HIER_STABLE4_PROGRESS_V1 Report completed 150-row challenger cubes.
if nargin<2 || strlength(string(scope))==0,scope="CLEAN_BRIDGE";end
scope=upper(string(scope));
if scope=="CLEAN_BRIDGE",scenarios=[23001 23004];elseif scope=="ALL_TFKM",scenarios=[23001 23004 23005 23007];else,error('Unknown scope.');end
root=fullfile(char(string(output_root)),'hier_stable4_full_pipeline_v1','candidates');
P=table();
for s=scenarios
    for rep=1:8
        f=fullfile(root,sprintf('hier_s%d_rep%02d.csv',s,rep));
        rows=0; valid=0;
        if isfile(f)
            T=readtable(f,'TextType','string');rows=height(T);valid=sum(logical(T.valid));
        end
        P=[P;table(s,rep,rows,valid,rows==150 && valid==150,'VariableNames',{'scenario_id','replicate','rows','valid_rows','complete'})]; %#ok<AGROW>
    end
end
disp(P);
fprintf('Complete: %d/%d jobs.\n',sum(P.complete),height(P));
end
