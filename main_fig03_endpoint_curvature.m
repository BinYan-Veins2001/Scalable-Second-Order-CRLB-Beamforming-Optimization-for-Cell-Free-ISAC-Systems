function archive=main_fig03_endpoint_curvature(options)
%MAIN_FIG03_ENDPOINT_CURVATURE Seven unchanged baseline methods, paired endpoints.
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
archive=paperfig_compute03(options); o=paperfig_options(options);
if o.export; paperfig_plot03(archive); end
end
