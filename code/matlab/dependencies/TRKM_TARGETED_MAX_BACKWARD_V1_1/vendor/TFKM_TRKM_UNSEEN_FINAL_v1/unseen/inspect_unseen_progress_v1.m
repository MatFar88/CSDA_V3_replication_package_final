function P=inspect_unseen_progress_v1(output_root)
%INSPECT_UNSEEN_PROGRESS_V1 Compact progress table for sessions A:H.
S=load(fullfile(char(string(output_root)),'study.mat'),'manifest'); M=S.manifest;
P=M(:,{'session','replicate','expected_total'});
P.actual_rows=zeros(height(P),1); P.complete=false(height(P),1);
for i=1:height(P)
    f=fullfile(char(string(output_root)),'candidates',sprintf('unseen_candidates_%s.csv',char(P.session(i))));
    if isfile(f)
        T=readtable(f); P.actual_rows(i)=height(T);
        P.complete(i)=height(T)==P.expected_total(i);
    end
end
disp(P)
end
