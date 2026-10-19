function download_paper_results(overwrite)
%DOWNLOAD_PAPER_RESULTS Fetch and verify the versioned numerical archives.
% download_paper_results(false) refuses to overwrite existing MAT results.
if nargin<1; overwrite=false; end
root=fileparts(mfilename('fullpath'));
manifest=jsondecode(fileread(fullfile(root,'release_manifest.json')));
asset=manifest.assets(strcmp({manifest.assets.role},'current_results'));
assert(numel(asset)==1,'Expected one current-results archive.');
if ~overwrite
    existing=dir(fullfile(root,'results','paper_fig*','data.mat'));
    assert(isempty(existing), ...
        'Stored paper data already exist. Replot directly or pass true to replace them.');
end
fprintf('Downloading %.1f MB of stored data from GitHub Releases.\n',asset.bytes/1e6);
archive=[tempname '.zip'];
cleanup=onCleanup(@()remove_temporary(archive)); %#ok<NASGU>
websave(archive,asset.url,weboptions('Timeout',1800));
digest=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(archive,'rb'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
while ~feof(fid)
    digest.update(fread(fid,1024*1024,'*uint8'));
end
clear closer
actual=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
assert(strcmp(actual,asset.sha256),'Download checksum mismatch; nothing was extracted.');
unzip(archive,root);
fprintf('Verified results extracted. Run main_replot_paper_figures(1:6).\n');
end

function remove_temporary(path)
if isfile(path); delete(path); end
end
