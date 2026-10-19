function fig = simtwc_plot_fig01()
%SIMTWC_PLOT_FIG01 Topology and shared-scale position-CRLB maps.

simtwc_setup(); d=simtwc_load_data(1); st=simtwc_ieee_style('double');
fig=figure('Units','inches','Position',[1 1 st.width_in 3.05]);
ax1=axes(fig,'Position',[0.075 0.17 0.29 0.73]); hold(ax1,'on');
scatter(ax1,d.scenario.ap_pos(:,1),d.scenario.ap_pos(:,2),58,'^', ...
    'filled','MarkerFaceColor',[0.00 0.45 0.70],'MarkerEdgeColor','k');
scatter(ax1,d.scenario.user_pos(:,1),d.scenario.user_pos(:,2),34,'o', ...
    'filled','MarkerFaceColor',[0.00 0.62 0.45],'MarkerEdgeColor','k');
scatter(ax1,d.scenario.target_pos(1),d.scenario.target_pos(2),100,'p', ...
    'filled','MarkerFaceColor',[0.80 0.12 0.20],'MarkerEdgeColor','k');
xlabel(ax1,'x (m)'); ylabel(ax1,'y (m)');
title(ax1,'(a) Network topology','FontWeight','normal');
axis(ax1,'equal'); grid(ax1,'on'); box(ax1,'on');
positions=[d.scenario.ap_pos;d.scenario.user_pos;d.scenario.target_pos];
lo=min(positions,[],1); hi=max(positions,[],1);
padding=0.08*max(hi-lo,ones(1,2));
xlim(ax1,[lo(1)-padding(1),hi(1)+padding(1)]);
ylim(ax1,[lo(2)-padding(2),hi(2)+padding(2)]);
legend(ax1,{'AP','User','Target'},'Location','southoutside', ...
    'NumColumns',3,'FontSize',st.legend_font_size,'Box','off');

ax2=axes(fig,'Position',[0.43 0.18 0.205 0.68]);
ax3=axes(fig,'Position',[0.67 0.18 0.205 0.68]);
maps=d.crlb_maps; finite_values=maps(isfinite(maps)&maps>0);
cl=[min(finite_values),max(finite_values)];
if cl(2)<=cl(1); cl(2)=cl(1)*(1+eps); end
meta=simtwc_method_registry(d.map_methods); names={meta.name};
for k=1:2
    ax=ternary(k==1,ax2,ax3);
    contourf(ax,d.grid_x,d.grid_y,maps(:,:,k),18,'LineStyle','none');
    hold(ax,'on');
    plot(ax,d.scenario.target_pos(1),d.scenario.target_pos(2),'p', ...
        'MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','k', ...
        'Clipping','off');
    set(ax,'ColorScale','log','CLim',cl);
    axis(ax,'equal');
    xlim(ax,[min(d.grid_x(:))-0.5,max(d.grid_x(:))+0.5]);
    ylim(ax,[min(d.grid_y(:))-0.5,max(d.grid_y(:))+0.5]);
    box(ax,'on');
    xlabel(ax,'x (m)'); if k==1; ylabel(ax,'y (m)'); else; ax.YTickLabel=[]; end
    title(ax,names{k},'FontWeight','normal');
end
colormap(fig,turbo(256));
cb=colorbar(ax3); cb.Position=[0.895 0.18 0.018 0.68];
cb.Label.String='Position CRLB (m^2)';
annotation(fig,'textbox',[0.43 0.88 0.45 0.07],'String', ...
    '(b) Spatial CRLB distribution','EdgeColor','none', ...
    'HorizontalAlignment','center','FontName',st.font_name,'FontSize',8);
simtwc_export_figure(fig,'fig01_topology_crlbmap');
end

function y=ternary(c,a,b)
if c; y=a; else; y=b; end
end
