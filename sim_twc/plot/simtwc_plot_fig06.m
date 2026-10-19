function fig = simtwc_plot_fig06()
%SIMTWC_PLOT_FIG06 Curvature-model and negative-curvature ablation.

simtwc_setup(); d=simtwc_load_data(6); st=simtwc_ieee_style('double');
nV=numel(d.variants); nR=numel(d.regimes); nT=numel(d.seeds);
crlb=nan(nV,nR,nT); kkt=nan(size(crlb)); runtime=nan(size(crlb));
for v=1:nV
    for r=1:nR
        for t=1:nT
            q=d.runs{v,r,t}; if ~ok(q); continue; end
            crlb(v,r,t)=q.metrics.crlb; kkt(v,r,t)=q.kkt.r_kkt;
            runtime(v,r,t)=q.runtime_median_s;
        end
    end
end
gap=nan(size(crlb)); ratio=nan(size(runtime));
for r=1:nR
    for t=1:nT
        gap(:,r,t)=100*(crlb(:,r,t)/crlb(1,r,t)-1);
        ratio(:,r,t)=runtime(:,r,t)/runtime(1,r,t);
    end
end
[gmu,gse]=simtwc_mean_se(gap,3);
[kmu,kse]=simtwc_mean_se(kkt,3);
[rmu,rse]=simtwc_mean_se(ratio,3);
colors=[0.80 0.12 0.20;0.95 0.55 0.20;0.00 0.45 0.70;0.45 0.65 0.85];
labels={'Exact + NC','Exact, NC off','GN + NC','GN, NC off'};
regimes={d.regimes.name};

fig=figure('Units','inches','Position',[1 1 st.width_in 2.9]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
grouped(nexttile(tl),gmu',gse',colors,regimes, ...
    '(a) CRLB change','Change from full COSTAR (%)',false,st);
grouped(nexttile(tl),kmu',kse',colors,regimes, ...
    '(b) Unified stationarity','r_{KKT}',true,st);
ax=nexttile(tl);
grouped(ax,rmu',rse',colors,regimes,'(c) Runtime ratio', ...
    'Runtime / full COSTAR',false,st);
yline(ax,1,':','Color',[0.3 0.3 0.3],'HandleVisibility','off', ...
    'LineWidth',st.reference_line_width);
lg=legend(ax,labels,'Location','southoutside','NumColumns',2,'Box','off');
lg.Layout.Tile='south';
simtwc_export_figure(fig,'fig06_curvature_ablation');
end

function grouped(ax,Y,E,colors,xlabels,ttl,ylab,use_log,st)
hold(ax,'on'); b=bar(ax,Y,'grouped');
for j=1:numel(b)
    b(j).FaceColor=colors(j,:);
    b(j).EdgeColor=st.bar_edge_color;
    b(j).LineWidth=st.bar_edge_width;
end
for j=1:numel(b)
    errorbar(ax,b(j).XEndPoints,Y(:,j),E(:,j), ...
        'Color',[0.25 0.25 0.25],'LineStyle','none', ...
        'CapSize',st.error_cap_size,'LineWidth',st.error_line_width);
end
if use_log; set(ax,'YScale','log'); end
set(ax,'XTick',1:numel(xlabels),'XTickLabel',xlabels, ...
    'XTickLabelRotation',22);
title(ax,ttl,'FontWeight','normal'); ylabel(ax,ylab); grid(ax,'on'); box(ax,'on');
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
