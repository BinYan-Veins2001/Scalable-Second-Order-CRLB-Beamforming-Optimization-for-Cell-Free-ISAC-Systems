function fig=simtwc_plot_fig11()
%SIMTWC_PLOT_FIG11 Accuracy and time reliability from stored statistics.
s=simtwc_setup(); d=simtwc_load_data(11);
st=simtwc_ieee_style('double'); meta=simtwc_method_registry(d.methods);
nA=numel(meta);
fig=figure('Units','inches','Position',[1 1 st.width_in 2.70]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); hold(ax,'on'); handles=gobjects(nA,1);
for a=1:nA
    [sx,sy]=stairs(d.accuracy_grid_pct,d.accuracy_success_pct(a,:));
    handles(a)=simtwc_plot_method_curve(ax,sx,sy,meta(a),st,[],[],1:44:numel(sx));
end
set(ax,'XScale','log'); xlim(ax,d.accuracy_grid_pct([1 end])); ylim(ax,[-3 103]);
xlabel(ax,'Allowed paired CRLB gap (\%)','Interpreter','latex');
ylabel(ax,'Feasible success rate (\%)','Interpreter','latex');
title(ax,'(a) Accuracy reliability','FontWeight','normal'); grid(ax,'on'); box(ax,'on');
ax=nexttile(tl); hold(ax,'on');
for a=1:nA
    [sx,sy]=stairs(d.time_grid_s,d.target_reached_pct(a,:));
    simtwc_plot_method_curve(ax,sx,sy,meta(a),st,[],[],1:56:numel(sx));
end
set(ax,'XScale','log'); xlim(ax,d.time_grid_s([1 end])); ylim(ax,[-3 103]);
xlabel(ax,'Wall-clock budget, $t$ (s)','Interpreter','latex');
ylabel(ax,'Target-reaching rate (\%)','Interpreter','latex');
title(ax,'(b) Time reliability','FontWeight','normal'); grid(ax,'on'); box(ax,'on');
lg=legend(handles,{meta.name},'NumColumns',nA,'Box','off', ...
    'FontSize',st.legend_font_size-0.2);
lg.Layout.Tile='south'; simtwc_export_figure(fig,'fig11_reliability');
writetable(d.summary,fullfile(s.figure_dir,'fig11_reliability_summary.csv'));
end
