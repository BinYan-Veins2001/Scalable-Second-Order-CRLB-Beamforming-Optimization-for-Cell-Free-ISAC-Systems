function fig = simtwc_plot_fig04()
%SIMTWC_PLOT_FIG04 Equal-accuracy timing and LaTeX timing table.

simtwc_setup(); d=simtwc_load_data(4); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods); nA=numel(meta);
[target_mu,target_se]=simtwc_mean_se(d.time_to_target_s,2);
[kkt_mu,kkt_se]=simtwc_mean_se(d.matched_kkt_time_s,2);

fig=figure('Units','inches','Position',[1 1 st.width_in 3.0]);
ax=axes(fig); hold(ax,'on');
Y=[target_mu,kkt_mu]; E=[target_se,kkt_se];
b=bar(ax,1:nA,Y,'grouped');
b(1).FaceColor=[0.00 0.45 0.70]; b(2).FaceColor=[0.80 0.12 0.20];
for j=1:2
    b(j).EdgeColor=st.bar_edge_color;
    b(j).LineWidth=st.bar_edge_width;
    x=b(j).XEndPoints;
    errorbar(ax,x,Y(:,j),E(:,j),'Color',[0.25 0.25 0.25], ...
        'LineStyle','none','LineWidth',st.error_line_width, ...
        'CapSize',st.error_cap_size);
end
set(ax,'YScale','log','XTick',1:nA,'XTickLabel',{meta.name}, ...
    'XTickLabelRotation',25);
ylabel(ax,'First-hit runtime (s)'); grid(ax,'on'); box(ax,'on');
legend(ax,{'Time-to-target','Matched KKT'},'Location','northwest', ...
    'Box','off','FontSize',st.legend_font_size);
simtwc_export_figure(fig,'fig04_matched_timing');
write_table(d,meta);
end

function write_table(d,meta)
s=simtwc_config(); [gmu,gse]=simtwc_mean_se(d.paired_crlb_gap_pct,2);
tt=median(d.time_to_target_s,2,'omitnan'); tk=median(d.matched_kkt_time_s,2,'omitnan');
path=fullfile(s.figure_dir,'tab_timing.tex'); fid=fopen(path,'w');
if fid<0; error('Cannot write %s.',path); end
c=onCleanup(@()fclose(fid));
fprintf(fid,'\\begin{tabular}{lccccc}\n\\toprule\n');
fprintf(fid,'Method & CRLB gap (\\%%) & $t_{1.01}$ (s) & $t_{\\rm KKT}$ (s) & KKT reach (\\%%) & Feasible (\\%%) \\\\\n');
fprintf(fid,'\\midrule\n');
for a=1:numel(meta)
    fprintf(fid,'%s & %.3f $\\pm$ %.3f & %s & %s & %.1f & %.1f \\\\\n', ...
        meta(a).name,gmu(a),gse(a),fmt(tt(a)),fmt(tk(a)), ...
        100*d.kkt_reach_rate(a),100*d.feasibility_rate(a));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
end

function s=fmt(x)
if isfinite(x); s=sprintf('%.4g',x); else; s='--'; end
end
