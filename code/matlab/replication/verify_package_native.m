function N=verify_package_native(root)
%VERIFY_PACKAGE_NATIVE Verify paths and hashes without Python or fitting.
f=fullfile(root,'documentation','SHA256SUMS.csv');assert(isfile(f),'Missing integrity manifest.');
T=readtable(f,'TextType','string');
for i=1:height(T)
    p=fullfile(root,strrep(char(T.path(i)),'/',filesep));
    assert(isfile(p),'Missing packaged file: %s',p);
    assert(strcmp(file_sha256(p),char(T.sha256(i))),'Changed packaged file: %s',p);
end
N=height(T);fprintf('Package integrity: %d files verified.\n',N);
end
