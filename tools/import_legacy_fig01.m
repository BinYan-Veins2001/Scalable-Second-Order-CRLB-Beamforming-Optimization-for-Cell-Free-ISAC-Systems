function destination=import_legacy_fig01(source)
%IMPORT_LEGACY_FIG01 Import the existing paper plot without running a sweep.
% Select plot arrays only; retain the original experiment specification.
assert(nargin==1 && isfile(source),'Supply the original fig15_data.mat archive.');
root=fileparts(fileparts(mfilename('fullpath')));
items=whos('-file',source); names={items.name};
keep=~endsWith(names,'_runs');
data=load(source,names{keep});
assert(data.complete,'The inherited scaling archive is incomplete.');
data.publication_provenance=struct( ...
    'source_filename','fig15_data.mat', ...
    'source_version','COSTAR_TWC_Simulations_ExactHVP_PTR_20260907', ...
    'role','Inherited paper Figure 1; not rerun with the October 19 solver', ...
    'raw_runs_omitted_from_plot_archive',true);
folder=fullfile(root,'results','paper_fig01');
if ~isfolder(folder); mkdir(folder); end
destination=fullfile(folder,'data.mat');
save(destination,'-struct','data','-v7.3');
fprintf('Imported existing Figure 1 plot arrays: %s\n',destination);
end
