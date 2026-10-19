function value=paperfig_signature()
%PAPERFIG_SIGNATURE Hash numerical source files; old experiment caches are excluded.
s=paperfig_setup();
folders={'optimization','config','baselines','channels','scenario','crlb','utils'};
digest=java.security.MessageDigest.getInstance('SHA-256');
for j=1:numel(folders)
    files=dir(fullfile(s.root,folders{j},'*.m'));
    for k=1:numel(files)
        txt=[folders{j} '/' files(k).name fileread(fullfile(files(k).folder,files(k).name))];
        digest.update(unicode2native(txt,'UTF-8'));
    end
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
