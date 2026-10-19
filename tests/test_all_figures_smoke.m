function report=test_all_figures_smoke()
%TEST_ALL_FIGURES_SMOKE Small-system end-to-end validation, NOT paper data.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
previous=getenv('COSTAR_SMOKE'); visibility=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@() restore(previous,visibility)); %#ok<NASGU>
setenv('COSTAR_SMOKE','1'); set(groot,'defaultFigureVisible','off');
[~,report]=main_run_all_twc_figures(1:19);
close all;
save(fullfile(root,'tests','figure_smoke_report.mat'),'report');
end

function restore(previous,visibility)
setenv('COSTAR_SMOKE',previous); set(groot,'defaultFigureVisible',visibility);
end
