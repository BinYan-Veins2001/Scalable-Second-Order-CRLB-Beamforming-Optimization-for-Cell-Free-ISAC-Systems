function fig = simtwc_plot_fig07()
%SIMTWC_PLOT_FIG07 Feasible physical CRLB, close competitors, time and QoS.

simtwc_setup(); d=simtwc_load_data(7); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods);
nA=numel(meta); nG=numel(d.gamma_dB); nT=numel(d.seeds);
crlb=nan(nA,nG,nT); runtime=nan(size(crlb)); feasible=false(size(crlb));
for a=1:nA
    for g=1:nG
        for t=1:nT
            q=d.runs{a,g,t}; if ~ok(q); continue; end
            crlb(a,g,t)=q.metrics.crlb; runtime(a,g,t)=q.runtime_median_s;
            feasible(a,g,t)=q.feasible;
        end
    end
end
isdr=find(strcmp(d.methods,'sdr_sca'),1); gap=nan(size(crlb));
for g=1:nG
    for t=1:nT
        if feasible(isdr,g,t)
            valid=feasible(:,g,t);
            gap(valid,g,t)=100*(crlb(valid,g,t)/crlb(isdr,g,t)-1);
        end
    end
end
[gmu,gse]=simtwc_mean_se(gap,3);
[rmu,rse]=simtwc_mean_se(runtime,3);
fr=100*mean(feasible,3);
crlb(~feasible)=NaN;
[cmu,cse]=simtwc_mean_se(crlb,3);

fig=figure('Units','inches','Position',[1 1 st.width_in 3.75]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl); h=curves(ax1,d.gamma_dB,cmu,cse,meta,true);
ylabel(ax1,'Feasible CRLB ($\mathrm{m}^2$)','Interpreter','latex');
title(ax1,'(a) Feasible CRLB','FontWeight','normal');
ix=find(ismember(d.methods,{'costar','mrfp_sca','sdr_sca'}));
ax2=nexttile(tl); curves(ax2,d.gamma_dB,gmu(ix,:),gse(ix,:),meta(ix),false);
ylabel(ax2,'Paired gap to SDR-SCA (\%)','Interpreter','latex');
title(ax2,'(b) Close-competitor gap','FontWeight','normal');
ax3=nexttile(tl); curves(ax3,d.gamma_dB,rmu,rse,meta,true);
ylabel(ax3,'Runtime (s)'); title(ax3,'(c) End-to-end runtime','FontWeight','normal');
ax4=nexttile(tl); curves(ax4,d.gamma_dB,fr,zeros(size(fr)),meta,false);
ylabel(ax4,'Feasibility rate (\%)','Interpreter','latex');
ylim(ax4,[-3 103]); yticks(ax4,0:20:100);
title(ax4,'(d) Feasibility','FontWeight','normal');
xlabel(ax3,'Target SINR threshold, $\Gamma$ (dB)','Interpreter','latex');
xlabel(ax4,'Target SINR threshold, $\Gamma$ (dB)','Interpreter','latex');
lg=legend(h,{meta.name},'NumColumns',numel(meta),'Box','off', ...
    'FontSize',st.legend_font_size-0.2);
lg.Layout.Tile='south';
simtwc_export_figure(fig,'fig07_sinr_sweep');
end

function h=curves(ax,x,mu,se,meta,use_log)
hold(ax,'on');
st=simtwc_ieee_style('double');
h=gobjects(numel(meta),1);
for a=1:numel(meta)
    h(a)=simtwc_plot_method_curve(ax,x,mu(a,:),meta(a),st,se(a,:));
end
if use_log; set(ax,'YScale','log'); end
grid(ax,'on'); box(ax,'on');
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
