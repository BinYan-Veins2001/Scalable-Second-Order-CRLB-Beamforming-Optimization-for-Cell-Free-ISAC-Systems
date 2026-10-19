function fig = simtwc_plot_fig12()
%SIMTWC_PLOT_FIG12 Endpoint first- and second-order stationarity audit.

simtwc_setup(); d=simtwc_load_data(12); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods);

fig=figure('Units','inches','Position',[1 1 st.width_in 2.85]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
diagnostic_panel(nexttile(tl),d.kkt_residual,meta,d.spec.sosp_tol, ...
    '(a) First-order stationarity','r_{KKT}',st);
diagnostic_panel(nexttile(tl),d.normalized_negative_curvature,meta, ...
    d.spec.sosp_tol,'(b) Negative-curvature residual', ...
    '[-\lambda_{min}^{H}]_+ / max(1, norm(g))',st);
diagnostic_panel(nexttile(tl),d.sosp_residual,meta,d.spec.sosp_tol, ...
    '(c) Combined SOSP certificate','r_{SOSP}',st);
simtwc_export_figure(fig,'fig12_second_order_certificate');
write_summary(d,meta);
end

function diagnostic_panel(ax,Y,meta,tol,ttl,ylab,st)
hold(ax,'on'); nA=numel(meta); floor_value=1e-12;
for a=1:nA
    y=Y(a,:); y=y(isfinite(y));
    if isempty(y); continue; end
    yp=max(y,floor_value);
    x=a+linspace(-0.13,0.13,numel(yp));
    plot(ax,x,yp,'LineStyle','none','Marker','o','MarkerSize',3, ...
        'MarkerEdgeColor',meta(a).color,'MarkerFaceColor','none', ...
        'Color',meta(a).color);
    q=prctile(yp,[25 50 75]);
    plot(ax,[a a],q([1 3]),'-','Color',meta(a).color, ...
        'LineWidth',1.7);
    plot(ax,a,q(2),meta(a).marker,'Color',meta(a).color, ...
        'MarkerFaceColor',meta(a).color,'MarkerSize',6, ...
        'LineWidth',1.1);
end
yline(ax,tol,'--','Tolerance','Color',[0.35 0.35 0.35], ...
    'LineWidth',st.reference_line_width,'LabelHorizontalAlignment','left');
set(ax,'YScale','log','XTick',1:nA,'XTickLabel',{meta.name}, ...
    'XTickLabelRotation',32,'XLim',[0.45 nA+0.55]);
finite_y=max(Y(isfinite(Y)),floor_value);
if isempty(finite_y); finite_y=tol; end
ylim(ax,[min(floor_value,tol/10) max(10*tol,2*max(finite_y))]);
ylabel(ax,ylab); title(ax,ttl,'FontWeight','normal');
grid(ax,'on'); box(ax,'on');
end

function write_summary(d,meta)
nA=numel(meta); method=strings(nA,1); median_kkt=nan(nA,1);
median_nc=nan(nA,1); median_sosp=nan(nA,1); sosp_pass_pct=nan(nA,1);
feasible_pct=nan(nA,1);
for a=1:nA
    method(a)=string(meta(a).name);
    median_kkt(a)=median(d.kkt_residual(a,:),'omitnan');
    median_nc(a)=median(d.normalized_negative_curvature(a,:),'omitnan');
    median_sosp(a)=median(d.sosp_residual(a,:),'omitnan');
    valid=isfinite(d.sosp_residual(a,:));
    if any(valid)
        sosp_pass_pct(a)=100*mean(d.sosp_residual(a,valid)<=d.spec.sosp_tol);
    end
    feasible_pct(a)=100*mean(d.feasible(a,:));
end
T=table(method,median_kkt,median_nc,median_sosp,sosp_pass_pct,feasible_pct);
s=simtwc_config(); writetable(T,fullfile(s.figure_dir, ...
    'fig12_second_order_certificate_summary.csv'));
end
