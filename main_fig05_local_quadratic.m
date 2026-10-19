function archive=main_fig05_local_quadratic(options)
%MAIN_FIG05_LOCAL_QUADRATIC Controlled conditional local Newton experiment.
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
archive=paperfig_compute05(options); o=paperfig_options(options);
if o.export; paperfig_plot05(archive); end
end
