function simtwc_save_atomic(path,data)
%SIMTWC_SAVE_ATOMIC Save one checkpoint without exposing a partial MAT file.

[folder,~,~]=fileparts(path);
if ~exist(folder,'dir'); mkdir(folder); end
tmp=[path '.tmp.mat'];
save(tmp,'-struct','data','-v7.3');
movefile(tmp,path,'f');
end
