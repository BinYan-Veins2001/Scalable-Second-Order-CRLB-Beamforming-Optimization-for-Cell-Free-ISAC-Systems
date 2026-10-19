function fig = simtwc_plot_fig03()
%SIMTWC_PLOT_FIG03 Paired CRLB gap, runtime, and unified KKT residual.

simtwc_setup(); d=simtwc_load_data(3); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods); nA=numel(meta); nT=numel(d.seeds);
crlb=nan(nA,nT); runtime=nan(nA,nT); kkt=nan(nA,nT); feasible=false(nA,nT);
for a=1:nA
    for t=1:nT
        r=d.runs{a,t}; if ~ok(r); continue; end
        crlb(a,t)=r.metrics.crlb; runtime(a,t)=r.runtime_median_s;
        kkt(a,t)=r.kkt.r_kkt; feasible(a,t)=r.feasible;
    end
end
isdr=find(strcmp(d.methods,'sdr_sca'),1); gap=nan(size(crlb));
for t=1:nT
    if feasible(isdr,t)
        valid=feasible(:,t); gap(valid,t)=100*(crlb(valid,t)/crlb(isdr,t)-1);
    end
end
[gap_mu,gap_se]=simtwc_mean_se(gap,2);
[rt_mu,rt_se]=simtwc_mean_se(runtime,2);
[kkt_mu,kkt_se]=simtwc_mean_se(kkt,2);

fig=figure('Units','inches','Position',[1 1 st.width_in 2.75]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
draw_panel(nexttile(tl),gap_mu,gap_se,meta,false, ...
    '(a) Paired CRLB gap','Gap to SDR-SCA (%)',st);
draw_panel(nexttile(tl),rt_mu,rt_se,meta,true, ...
    '(b) End-to-end runtime','Runtime (s)',st);
draw_panel(nexttile(tl),kkt_mu,kkt_se,meta,true, ...
    '(c) Unified KKT residual','r_{KKT}',st);
simtwc_export_figure(fig,'fig03_main_comparison');
end

function draw_panel(ax,mu,se,meta,use_log,ttl,ylab,st)
hold(ax,'on'); n=numel(meta); b=bar(ax,1:n,mu,0.72,'FaceColor','flat');
b.CData=vertcat(meta.color);
b.EdgeColor=st.bar_edge_color; b.LineWidth=st.bar_edge_width;
errorbar(ax,1:n,mu,se,'Color',[0.25 0.25 0.25], ...
    'LineStyle','none','LineWidth',st.error_line_width, ...
    'CapSize',st.error_cap_size);
if use_log; set(ax,'YScale','log'); end
set(ax,'XTick',1:n,'XTickLabel',{meta.name},'XTickLabelRotation',32);
ylabel(ax,ylab); title(ax,ttl,'FontWeight','normal'); grid(ax,'on'); box(ax,'on');
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
