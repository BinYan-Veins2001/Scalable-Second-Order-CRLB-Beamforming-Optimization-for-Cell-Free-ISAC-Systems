function fig = simtwc_plot_fig09()
%SIMTWC_PLOT_FIG09 Power sweep and COSTAR active-constraint composition.

simtwc_setup(); d=simtwc_load_data(9); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods); nA=numel(meta); nP=numel(d.power_dBm);
nT=numel(d.seeds); C=nan(nA,nP,nT); R=nan(size(C));
for a=1:nA
    for p=1:nP
        for t=1:nT
            q=d.runs{a,p,t}; if ~ok(q); continue; end
            C(a,p,t)=q.metrics.crlb; R(a,p,t)=q.runtime_median_s;
        end
    end
end
[cmu,cse]=simtwc_mean_se(C,3); [rmu,rse]=simtwc_mean_se(R,3);
asinr=mean(d.active_sinr,2,'omitnan'); apow=mean(d.active_power,2,'omitnan');

fig=figure('Units','inches','Position',[1 1 st.width_in 2.85]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); method_curves(ax,d.power_dBm,cmu,cse,meta,true);
ylabel(ax,'CRLB (m^2)'); title(ax,'(a) Sensing accuracy','FontWeight','normal');
legend(ax,{meta.name},'Location','best','Box','off', ...
    'FontSize',st.legend_font_size);
ax=nexttile(tl); method_curves(ax,d.power_dBm,rmu,rse,meta,true);
ylabel(ax,'Runtime (s)'); title(ax,'(b) End-to-end runtime','FontWeight','normal');
ax=nexttile(tl); b=bar(ax,d.power_dBm,[asinr apow],0.7,'stacked');
for j=1:numel(b)
    b(j).EdgeColor=st.bar_edge_color;
    b(j).LineWidth=st.bar_edge_width;
end
xlabel(ax,'P_{max} (dBm)'); ylabel(ax,'Active constraints');
title(ax,'(c) COSTAR active set','FontWeight','normal'); grid(ax,'on'); box(ax,'on');
legend(ax,{'SINR','AP power'},'Location','northwest','Box','off', ...
    'FontSize',st.legend_font_size);
simtwc_export_figure(fig,'fig09_power_sweep');
end

function method_curves(ax,x,mu,se,meta,use_log)
hold(ax,'on');
st=simtwc_ieee_style('double');
for a=1:numel(meta)
    simtwc_plot_method_curve(ax,x,mu(a,:),meta(a),st,se(a,:));
end
if use_log; set(ax,'YScale','log'); end
xlabel(ax,'P_{max} (dBm)'); grid(ax,'on'); box(ax,'on');
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
