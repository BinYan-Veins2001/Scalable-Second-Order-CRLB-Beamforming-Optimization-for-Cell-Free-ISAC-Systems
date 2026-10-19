function result=main_fig06_structured_curvature_efficiency(options)
%MAIN_FIG06_STRUCTURED_CURVATURE_EFFICIENCY New paper Fig. 6, measured data.
% Run in an otherwise idle MATLAB process; a fresh -singleCompThread process
% is recommended for publication timings. This entry also enforces one
% computational thread temporarily when called in an ordinary MATLAB session.
%   main_fig06_structured_curvature_efficiency
%   main_fig06_structured_curvature_efficiency(struct('compute',false))
% Figure 1 is not evaluated by this entry point.
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
if ~isfield(options,'compute'); options.compute=true; end
if ~isfield(options,'plot'); options.plot=true; end
if options.compute
    previous_threads=maxNumCompThreads(1);
    cleanup_threads=onCleanup(@()maxNumCompThreads(previous_threads)); %#ok<NASGU>
    result=paperfig_compute_efficiency(options);
else
    a=load(fullfile(root,'results','paper_fig06','data.mat'),'result');
    result=a.result;
end
if options.plot; paperfig_plot_efficiency(result); end
end
