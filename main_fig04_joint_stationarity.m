function archive=main_fig04_joint_stationarity(options)
%MAIN_FIG04_JOINT_STATIONARITY Independent dense audit and full-stop guard records.
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
archive=paperfig_compute04(options); o=paperfig_options(options);
if o.export; paperfig_plot04(archive); end
end
