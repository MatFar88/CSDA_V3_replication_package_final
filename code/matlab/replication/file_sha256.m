function hex=file_sha256(filename)
%FILE_SHA256 Read-only SHA-256. Requires MATLAB with a JVM.
assert(usejava('jvm'),'The integrity verifier requires a JVM-enabled MATLAB session.');
fid=fopen(filename,'rb');assert(fid>0,'Cannot read %s.',filename);
clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1048576,'*uint8');
    if ~isempty(bytes),md.update(typecast(bytes(:),'int8'));end
end
raw=typecast(int8(md.digest()),'uint8');
hex=lower(reshape(dec2hex(raw,2).',1,[]));
end
