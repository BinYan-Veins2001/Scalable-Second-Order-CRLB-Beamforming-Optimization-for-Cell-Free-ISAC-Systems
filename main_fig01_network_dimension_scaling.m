function archive=main_fig01_network_dimension_scaling(options)
%MAIN_FIG01_NETWORK_DIMENSION_SCALING Full scaling sweep; run explicitly only.
% The standard paper-figure runner deliberately excludes this slow sweep.
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
s=paperfig_setup();
% All seven methods and all publication grid points run by default.
% Set options.sdr_sca_max_M explicitly to impose an AP-sweep resource cap.
archive=paperfig_compute01(options);
folder=fullfile(s.results_dir,'paper_fig01'); if ~isfolder(folder); mkdir(folder); end
copyfile(archive,fullfile(folder,'data.mat'));
archive=fullfile(folder,'data.mat'); paperfig_plot01(archive);
end
