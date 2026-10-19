function fig = simtwc_plot_fig05()
%SIMTWC_PLOT_FIG05 COSTAR tangent-gradient, radius, and active-set traces.

simtwc_setup(); d=simtwc_load_data(5); st=simtwc_ieee_style('single');
fig=figure('Units','inches','Position',[1 1 st.width_in 6.4]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); g=max(d.tangent_grad_norm,realmin);
semilogy(ax,d.iter,g,'-o','Color',[0.80 0.12 0.20], ...
    'LineWidth',st.line_width,'MarkerSize',st.marker_size, ...
    'MarkerIndices',marker_indices(numel(g))); grid(ax,'on'); box(ax,'on');
ylabel(ax,'||g_k^{tan}||_2'); title(ax,'(a) Tangent-gradient norm','FontWeight','normal');
annotate_slopes(ax,d.iter,g);

ax=nexttile(tl);
semilogy(ax,d.iter,max(d.delta,realmin),'-s','Color',[0.00 0.45 0.70], ...
    'LineWidth',st.line_width,'MarkerSize',st.marker_size, ...
    'MarkerIndices',marker_indices(numel(d.delta)));
grid(ax,'on'); box(ax,'on'); ylabel(ax,'\Delta_k');
title(ax,'(b) Trust-region radius','FontWeight','normal');

ax=nexttile(tl);
stairs(ax,d.iter,d.active_count,'LineWidth',st.line_width, ...
    'Color',[0.00 0.62 0.45]);
hold(ax,'on'); changed=find(d.active_change_count>0);
plot(ax,d.iter(changed),d.active_count(changed),'o','Color',[0.84 0.37 0], ...
    'LineWidth',st.line_width,'MarkerFaceColor','w', ...
    'MarkerSize',st.marker_size);
grid(ax,'on'); box(ax,'on'); xlabel(ax,'RTR iteration, k'); ylabel(ax,'|A_k|');
title(ax,'(c) Active-set size','FontWeight','normal');
simtwc_export_figure(fig,'fig05_convergence_diag');
end

function idx=marker_indices(n)
idx=unique(round(linspace(1,n,min(10,n))));
end

function annotate_slopes(ax,k,g)
valid=find(isfinite(g)&g>0); if numel(valid)<8; return; end
n=numel(valid); early=valid(1:max(3,floor(n/3)));
late=valid(max(1,n-floor(n/3)+1):n);
p1=polyfit(k(early),log10(g(early)),1);
p2=polyfit(k(late),log10(g(late)),1);
text(ax,0.98,0.92,sprintf('log-slope: %.2g -> %.2g',p1(1),p2(1)), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',7,'BackgroundColor','w');
end
