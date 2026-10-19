function fig = simtwc_plot_fig08()
%SIMTWC_PLOT_FIG08 Compact and large-area AP scaling.

simtwc_setup(); d=simtwc_load_data(8); st=simtwc_ieee_style('double');
[cc,cr,cf]=collect(d.compact_runs);
[lc,lr,lf]=collect(d.large_runs);
cm=simtwc_method_registry(d.compact_methods);
lm=simtwc_method_registry(d.large_methods);

fig=figure('Units','inches','Position',[1 1 st.width_in 5.05]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl);
hcompact=line_panel(ax1,d.spec.compact.M_values,cc,cm,true, ...
    '(a) Compact: CRLB','CRLB (m^2)');
line_panel(nexttile(tl),d.spec.compact.M_values,cr,cm,true, ...
    '(b) Compact: runtime','Runtime (s)');
line_panel(nexttile(tl),d.spec.compact.M_values,100*cf,cm,false, ...
    '(c) Compact: feasibility','Feasibility rate (%)'); ylim([0 100]);
line_panel(nexttile(tl),d.spec.large.M_values,lc,lm,true, ...
    '(d) Large area: CRLB','CRLB (m^2)');
line_panel(nexttile(tl),d.spec.large.M_values,lr,lm,true, ...
    '(e) Large area: runtime','Runtime (s)');
ax=nexttile(tl);
line_panel(ax,d.spec.large.M_values,100*lf,lm,false, ...
    '(f) Large area: feasibility','Feasibility rate (%)'); ylim(ax,[0 100]);
lg=legend(ax,hcompact,{cm.name},'Location','southoutside','NumColumns',4, ...
    'Box','off','FontSize',st.legend_font_size);
lg.Layout.Tile='south';
simtwc_export_figure(fig,'fig08_ap_scaling');
end

function [crlb,runtime,feas]=collect(runs)
[nA,nM,nT]=size(runs); C=nan(nA,nM,nT); R=nan(size(C)); F=nan(size(C));
for a=1:nA
    for m=1:nM
        for t=1:nT
            q=runs{a,m,t};
            if isstruct(q) && isfield(q,'status') && strcmp(q.status,'failed')
                F(a,m,t)=0;
                continue;
            end
            if ~ok(q); continue; end
            C(a,m,t)=q.metrics.crlb; R(a,m,t)=q.runtime_median_s; F(a,m,t)=q.feasible;
        end
    end
end
crlb=mean(C,3,'omitnan'); runtime=mean(R,3,'omitnan');
feas=mean(F,3,'omitnan');
end

function h=line_panel(ax,x,Y,meta,use_log,ttl,ylab)
hold(ax,'on');
st=simtwc_ieee_style('double');
h=gobjects(numel(meta),1);
for a=1:numel(meta)
    line_style=meta(a).line;
    if strcmp(meta(a).key,'sdr_sca'); line_style='none'; end
    h(a)=simtwc_plot_method_curve( ...
        ax,x,Y(a,:),meta(a),st,[],line_style);
end
if use_log; set(ax,'YScale','log'); end
xlabel(ax,'Number of APs, M'); ylabel(ax,ylab);
title(ax,ttl,'FontWeight','normal'); grid(ax,'on'); box(ax,'on');
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
