function fig = simtwc_plot_fig02()
%SIMTWC_PLOT_FIG02 CRLB versus end-to-end wall-clock time.

simtwc_setup(); d=simtwc_load_data(2); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods);
fig=figure('Units','inches','Position',[1 1 st.width_in 3.15]);
ax=axes(fig); hold(ax,'on');
for a=1:numel(meta)
    r=d.runs{a}; if ~ok(r); continue; end
    n=numel(r.history.time_s); marks=unique(round(linspace(1,n,min(9,n))));
    simtwc_plot_method_curve(ax,r.history.time_s, ...
        max(r.history.crlb,realmin),meta(a),st,[],[],marks);
end
ic=find(strcmp(d.methods,'costar'),1);
if ~isempty(ic) && ok(d.runs{ic}) && isfield(d.runs{ic}.history,'rtr_start_time_s')
    xline(ax,d.runs{ic}.history.rtr_start_time_s,'--','RTR starts', ...
        'Color',[0.35 0.35 0.35],'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','right','LabelOrientation','horizontal', ...
        'HandleVisibility','off','LineWidth',st.reference_line_width);
end
if ~isempty(ic) && ok(d.runs{ic}) && isfield(d.runs{ic}.history,'terminal_start_time_s')
    xline(ax,d.runs{ic}.history.terminal_start_time_s,':','Accuracy stage', ...
        'Color',[0.65 0.1 0.2],'LabelVerticalAlignment','middle', ...
        'LabelHorizontalAlignment','left','LabelOrientation','horizontal', ...
        'HandleVisibility','off','LineWidth',st.reference_line_width);
end
set(ax,'YScale','log');
xlabel(ax,'Runtime (s)'); ylabel(ax,'CRLB ($\mathrm{m}^2$)','Interpreter','latex');
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','northeast','NumColumns',2,'Box','off', ...
    'FontSize',st.legend_font_size);
simtwc_export_figure(fig,'fig02_crlb_runtime');
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
