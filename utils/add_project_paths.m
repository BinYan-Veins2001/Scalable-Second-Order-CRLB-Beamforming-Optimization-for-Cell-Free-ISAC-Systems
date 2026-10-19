function add_project_paths(root_dir)
%ADD_PROJECT_PATHS Add all subfolders to MATLAB path.
if nargin < 1 || isempty(root_dir)
    root_dir = pwd;
end
addpath(genpath(root_dir));
end
