function s = simtwc_setup()
%SIMTWC_SETUP Add package paths and create output directories.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root),'-begin');
add_project_paths(root);
s = simtwc_config();
if ~exist(s.data_dir,'dir'); mkdir(s.data_dir); end
if ~exist(s.figure_dir,'dir'); mkdir(s.figure_dir); end
end
