function fig = simtwc_plot_fig13()
%SIMTWC_PLOT_FIG13 Conditional local quadratic-convergence diagnostic.

simtwc_setup(); d=simtwc_load_data(13); st=simtwc_ieee_style('double');
costar=simtwc_method_registry({'costar'}); color=costar.color;
fig=figure('Units','inches','Position',[1 1 st.width_in 2.85]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); keep=d.trajectory_keep & isfinite(d.trajectory_error);
semilogy(ax,find(keep)-1,max(d.trajectory_error(keep),1e-14),'-d', ...
    'Color',color,'MarkerFaceColor',color,'LineWidth',st.line_width);
xlabel(ax,'Stored AS-RTR state index'); ylabel(ax,'Phase-orbit error, e_k');
title(ax,'(a) Reference trajectory','FontWeight','normal'); grid(ax,'on'); box(ax,'on');

ax=nexttile(tl); hold(ax,'on'); valid=isfinite(d.e0)&isfinite(d.e1)&d.e0>0&d.e1>0;
loglog(ax,d.e0(valid),d.e1(valid),'d','Color',color, ...
    'MarkerFaceColor',color,'MarkerSize',6,'LineStyle','none');
if any(valid)
    lo=min([d.e0(valid) d.e1(valid)]); hi=max([d.e0(valid) d.e1(valid)]);
    xx=logspace(log10(max(lo,1e-14)),log10(max(hi,lo*10)),100);
    loglog(ax,xx,xx,'--','Color',[0.35 0.35 0.35], ...
        'LineWidth',st.reference_line_width);
    scale=median(d.e1(valid)./max(d.e0(valid).^2,realmin),'omitnan');
    loglog(ax,xx,scale*xx.^2,':','Color',[0.10 0.35 0.70], ...
        'LineWidth',st.line_width);
end
set(ax,'XScale','log','YScale','log');
xlabel(ax,'e_k'); ylabel(ax,'e_{k+1}');
if d.local_window_satisfied
    panel_title=sprintf('(b) One-step map (slope %.2f)',d.loglog_slope);
else
    panel_title='(b) One-step map (outside local window)';
end
title(ax,panel_title,'FontWeight','normal'); grid(ax,'on'); box(ax,'on');
legend(ax,{'Measured','Linear','Quadratic'},'Location','best','Box','off', ...
    'FontSize',st.legend_font_size);

ax=nexttile(tl); ratio=d.e1./max(d.e0.^2,realmin);
loglog(ax,d.e0(valid),ratio(valid),'d-','Color',color, ...
    'MarkerFaceColor',color,'LineWidth',st.line_width);
xlabel(ax,'e_k'); ylabel(ax,'e_{k+1}/e_k^2');
title(ax,'(c) Quadratic-ratio stability','FontWeight','normal');
grid(ax,'on'); box(ax,'on');
annotation_text(ax,d,valid);
simtwc_export_figure(fig,'fig13_local_quadratic');
write_summary(d);
end

function annotation_text(ax,d,valid)
n=max(nnz(valid),1);
txt=sprintf(['Accepted %.0f%%\nFeasible %.0f%%\nActive set stable %.0f%%' ...
    '\nInterior step %.0f%%'],100*nnz(d.accepted&valid)/n, ...
    100*nnz(d.feasible&valid)/n,100*nnz(d.active_set_stable&valid)/n, ...
    100*nnz(d.interior_step&valid)/n);
text(ax,0.04,0.96,txt,'Units','normalized','VerticalAlignment','top', ...
    'BackgroundColor','w','EdgeColor',[0.75 0.75 0.75],'Margin',3);
end

function write_summary(d)
valid=isfinite(d.e0)&isfinite(d.e1);
T=table(d.e0(valid).',d.e1(valid).',d.accepted(valid).', ...
    d.feasible(valid).',d.active_set_stable(valid).', ...
    d.interior_step(valid).','VariableNames', ...
    {'e_k','e_k_plus_1','accepted','feasible','active_set_stable','interior_step'});
s=simtwc_config(); writetable(T,fullfile(s.figure_dir, ...
    'fig13_local_quadratic_summary.csv'));
end
