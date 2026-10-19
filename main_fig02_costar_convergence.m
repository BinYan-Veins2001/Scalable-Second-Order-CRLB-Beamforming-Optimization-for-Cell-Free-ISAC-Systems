function archive=main_fig02_costar_convergence(options)
%MAIN_FIG02_COSTAR_CONVERGENCE Current solver, complete KKT, horizontal Ritz history.
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
archive=paperfig_compute02(options); o=paperfig_options(options);
if o.export; paperfig_plot02(archive); end
end
