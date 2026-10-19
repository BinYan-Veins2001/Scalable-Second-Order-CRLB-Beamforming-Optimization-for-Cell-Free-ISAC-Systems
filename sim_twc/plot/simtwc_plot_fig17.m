function fig = simtwc_plot_fig17()
%SIMTWC_PLOT_FIG17 Endpoint multipliers and corresponding constraints.

simtwc_setup(); d=simtwc_load_data(17); st=simtwc_ieee_style('single');
fig=figure('Units','inches','Position',[1 1 st.width_in 4.05]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

plot_family(nexttile(tl),d.lambda(d.sinr_indices,:),d.sinr_margin_dB, ...
    'User index, $u$','SINR margin (dB)','(a) SINR constraints',st,true,true);
plot_family(nexttile(tl),d.lambda(d.power_indices,:), ...
    d.constraints(d.power_indices,:),'AP index, $m$', ...
    'Normalized constraint, $c_j$', ...
    '(b) Per-AP power constraints',st,false,true);
simtwc_export_figure(fig,'fig17_endpoint_multipliers');
write_summary(d);
end

function plot_family(ax,L,margin,xname,right_label,panel_title,st,physical_margin,show_left_label)
muL=mean(L,2,'omitnan');
muMargin=mean(margin,2,'omitnan');
if physical_margin
    muMargin(abs(muMargin)<1e-10)=0;
end
x=1:size(L,1); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

yyaxis(ax,'left');
bar(ax,x,muL,0.68,'FaceColor',[0.00 0.45 0.70], ...
    'EdgeColor',[0.00 0.30 0.48],'LineWidth',0.55);
if show_left_label
    ylabel(ax,'Lagrange multiplier, $\lambda_j$','Interpreter','latex');
end
ax.YAxis(1).Color=[0.00 0.45 0.70];

yyaxis(ax,'right');
plot(ax,x,muMargin,'-d','Color',[0.80 0.12 0.20], ...
    'MarkerFaceColor','w','LineWidth',st.line_width, ...
    'MarkerSize',st.marker_size);
yline(ax,0,':','Color',[0.30 0.30 0.30], ...
    'LineWidth',st.reference_line_width,'HandleVisibility','off');
ylabel(ax,right_label,'Interpreter','latex');
ax.YAxis(2).Color=[0.80 0.12 0.20];
if physical_margin; set_margin_limits(ax,muMargin,zeros(size(muMargin))); end
xlabel(ax,xname,'Interpreter','latex'); xlim(ax,[0.35 numel(x)+0.65]); xticks(ax,x);
title(ax,panel_title,'FontWeight','normal');
end

function set_margin_limits(ax,mu,se)
lower=min(mu-se,[],'omitnan'); upper=max(mu+se,[],'omitnan');
span=max([abs(lower),abs(upper)]);
if ~isfinite(span) || span<1e-8
    ylim(ax,[-0.05 0.05]);
elseif lower>=-1e-10
    ylim(ax,[-0.08*span 1.12*span]);
else
    pad=0.08*max(upper-lower,span);
    ylim(ax,[lower-pad upper+pad]);
end
end

function write_summary(d)
nT=size(d.lambda,2); J=size(d.lambda,1);
kind=[repmat("SINR",numel(d.sinr_indices),1); ...
    repmat("AP power",numel(d.power_indices),1)];
local_index=[(1:numel(d.sinr_indices)).';(1:numel(d.power_indices)).'];
constraint_index=(1:J).';
mean_lambda=mean(d.lambda,2,'omitnan');
mean_constraint=mean(d.constraints,2,'omitnan');
mean_complementarity=mean(d.complementarity,2,'omitnan');
active_rate=100*mean(d.active,2,'omitnan');
mean_physical_margin=[mean(d.sinr_margin_dB,2,'omitnan'); ...
    mean(d.power_margin_percent,2,'omitnan')];
T=table(constraint_index,kind,local_index,mean_lambda,mean_constraint, ...
    mean_physical_margin,mean_complementarity,active_rate,repmat(nT,J,1), ...
    'VariableNames',{'constraint_index','type','local_index', ...
    'mean_lambda','mean_normalized_constraint','mean_physical_margin', ...
    'mean_abs_lambda_c','active_rate_percent','requested_trials'});
s=simtwc_config();
writetable(T,fullfile(s.figure_dir,'fig17_endpoint_multipliers_summary.csv'));
end
