function fig = simtwc_plot_fig10()
%SIMTWC_PLOT_FIG10 Active-set threshold sensitivity and parameter table.

simtwc_setup(); d=simtwc_load_data(10); st=simtwc_ieee_style('double');
nQ=numel(d.tau); nT=numel(d.seeds);
C=nan(nQ,nT); K=nan(size(C)); J=nan(size(C)); I=nan(size(C));
for q=1:nQ
    for t=1:nT
        r=d.runs{q,t}; if ~ok(r); continue; end
        C(q,t)=r.metrics.crlb; K(q,t)=r.kkt.r_kkt;
        J(q,t)=r.active_jitter; I(q,t)=r.rtr_iterations;
    end
end
[cmu,cse]=simtwc_mean_se(C,2); [kmu,kse]=simtwc_mean_se(K,2);
[jmu,jse]=simtwc_mean_se(J,2); [imu,ise]=simtwc_mean_se(I,2);

fig=figure('Units','inches','Position',[1 1 st.width_in 2.9]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); yyaxis(ax,'left');
metric_curve(ax,d.tau,cmu,cse,'-','o',[0.00 0.45 0.70],st);
set(ax,'XScale','log'); ylabel(ax,'Final CRLB (m^2)');
% Avoid exaggerating changes of a few parts per million by autoscaling.
baseline=mean(cmu,'omitnan');
ylim(ax,[min(cmu)*0.98 max(cmu)*1.02]);
if st.show_curve_se; ylim(ax,[min(cmu-cse)*0.98 max(cmu+cse)*1.02]); end
ytickformat(ax,'%.4f');
yyaxis(ax,'right');
metric_curve(ax,d.tau,kmu,kse,'--','s',[0.80 0.12 0.20],st);
set(ax,'YScale','log'); ylabel(ax,'r_{KKT}');
xlabel(ax,'\tau_\Gamma = \tau_P'); title(ax,'(a) Accuracy and stationarity','FontWeight','normal');
grid(ax,'on'); box(ax,'on');

ax=nexttile(tl); yyaxis(ax,'left');
metric_curve(ax,d.tau,jmu,jse,'-','d',[0.84 0.37 0.00],st);
set(ax,'XScale','log'); ylabel(ax,'Working-set transitions');
yyaxis(ax,'right');
metric_curve(ax,d.tau,imu,ise,'--','^',[0.00 0.62 0.45],st);
ylabel(ax,'RTR iterations'); xlabel(ax,'\tau_\Gamma = \tau_P');
title(ax,'(b) Active-set stability','FontWeight','normal'); grid(ax,'on'); box(ax,'on');
simtwc_export_figure(fig,'fig10_hyperparam');
write_table(d,cmu,cse,kmu,kse,jmu,jse,imu,ise);
end

function metric_curve(ax,x,y,se,line_style,marker,color,st)
hold(ax,'on');
meta=struct('color',color,'line',line_style,'marker',marker,'name','');
simtwc_plot_method_curve(ax,x,y,meta,st,se);
ax.YColor=color;
end

function write_table(d,cmu,cse,kmu,kse,jmu,jse,imu,ise)
s=simtwc_config(); path=fullfile(s.figure_dir,'tab_params.tex');
fid=fopen(path,'w'); if fid<0; error('Cannot write %s.',path); end
c=onCleanup(@()fclose(fid));
fprintf(fid,'\\begin{tabular}{ccccc}\n\\toprule\n');
fprintf(fid,'$\\tau_{\\Gamma}=\\tau_P$ & CRLB & $r_{\\rm KKT}$ & Transitions & Iterations \\\\\n');
fprintf(fid,'\\midrule\n');
for q=1:numel(d.tau)
    fprintf(fid,'%.0e & %.4g $\\pm$ %.2g & %.3g $\\pm$ %.2g & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n', ...
        d.tau(q),cmu(q),cse(q),kmu(q),kse(q),jmu(q),jse(q),imu(q),ise(q));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
