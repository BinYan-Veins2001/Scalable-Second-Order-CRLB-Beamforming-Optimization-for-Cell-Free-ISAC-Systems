function archive = simtwc_run_figure(id)
%SIMTWC_RUN_FIGURE Compute, checkpoint, plot, and export one figure.

simtwc_setup();
archive=simtwc_compute_figure(id);
simtwc_plot_figure(id);
end
