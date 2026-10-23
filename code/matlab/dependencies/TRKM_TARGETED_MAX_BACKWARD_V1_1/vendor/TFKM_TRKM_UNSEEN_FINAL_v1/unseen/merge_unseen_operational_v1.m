function M=merge_unseen_operational_v1(output_root)
%MERGE_UNSEEN_OPERATIONAL_V1 Merge A:H operational candidate files only.
root=char(string(output_root));
S=load(fullfile(root,'study.mat'),'manifest'); manifest=S.manifest;
M=table();
for i=1:height(manifest)
    ses=string(manifest.session(i));
    f=fullfile(root,'candidates',sprintf('unseen_candidates_%s.csv',char(ses)));
    assert(isfile(f),'Missing unseen session file: %s',f);
    T=readtable(f,'TextType','string');
    assert(height(T)==double(manifest.expected_total(i)), ...
        'Session %s incomplete: %d/%d.',ses,height(T),manifest.expected_total(i));
    M=[M;T]; %#ok<AGROW>
end
M=sortrows(M,{'scenario_id','replicate','method','r','g','alpha'});
key=[double(M.scenario_id),double(M.replicate),double(M.r),double(M.g),double(M.alpha)];
method=string(M.method);
assert(height(M)==9600,'Expected 9600 total candidates, got %d.',height(M));
assert(sum(method=="REG_TFKM_HD")==4800 && sum(method=="REG_TRKM_HD")==4800, ...
    'Expected 4800 candidates per method.');
assert(size(unique([key,double(method=="REG_TRKM_HD")],'rows'),1)==height(M), ...
    'Duplicate candidate keys detected.');
assert(all(logical(M.valid)),'At least one invalid candidate exists.');
out=fullfile(root,'operational_selection','operational_candidates_all.csv');
writetable(M,out);
fprintf('Merged unseen operational cube: %d rows -> %s\n',height(M),out);
end
