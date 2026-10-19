function fig = simtwc_plot_fig14()
%SIMTWC_PLOT_FIG14 Separate U and Nt dimension-scaling results.

simtwc_setup(); d=simtwc_load_data(14); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods);
[uc,use]=simtwc_mean_se(d.user_crlb,3);
[ur,urse]=simtwc_mean_se(d.user_runtime_s,3);
[ac,ase]=simtwc_mean_se(d.antenna_crlb,3);
[ar,arse]=simtwc_mean_se(d.antenna_runtime_s,3);

fig=figure('Units','inches','Position',[1 1 st.width_in 5.05]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); h=method_panel(ax,d.U_values,uc,use,meta,true, ...
    'Number of users, U','CRLB (m^2)','(a) User scaling: CRLB',st);
method_panel(nexttile(tl),d.U_values,ur,urse,meta,true, ...
    'Number of users, U','Runtime (s)','(b) User scaling: runtime',st);
work_panel(nexttile(tl),d.U_values,d.user_hvp_calls, ...
    d.user_tcg_iterations,'Number of users, U','(c) COSTAR work',st);
method_panel(nexttile(tl),d.Nt_values,ac,ase,meta,true, ...
    'Antennas per AP, N_t=N_r','CRLB (m^2)', ...
    '(d) Antenna scaling: CRLB',st);
method_panel(nexttile(tl),d.Nt_values,ar,arse,meta,true, ...
    'Antennas per AP, N_t=N_r','Runtime (s)', ...
    '(e) Antenna scaling: runtime',st);
axlast=nexttile(tl);
work_panel(axlast,d.Nt_values,d.antenna_hvp_calls, ...
    d.antenna_tcg_iterations,'Antennas per AP, N_t=N_r', ...
    '(f) COSTAR work',st);
lg=legend(axlast,h,{meta.name},'Orientation','horizontal','Box','off', ...
    'FontSize',st.legend_font_size,'NumColumns',numel(meta));
lg.Layout.Tile='south';
simtwc_export_figure(fig,'fig14_dimension_scaling');
write_summary(d,meta);
end

function h=method_panel(ax,x,mu,se,meta,use_log,xlab,ylab,ttl,st)
hold(ax,'on'); h=gobjects(numel(meta),1);
for a=1:numel(meta)
    err=[]; if st.show_curve_se; err=se(a,:); end
    h(a)=simtwc_plot_method_curve(ax,x,mu(a,:),meta(a),st,err);
end
if use_log; set(ax,'YScale','log'); end
xlabel(ax,xlab); ylabel(ax,ylab); title(ax,ttl,'FontWeight','normal');
grid(ax,'on'); box(ax,'on');
end

function work_panel(ax,x,H,I,xlab,ttl,st)
hold(ax,'on'); h=mean(H,[1 3],'omitnan'); q=mean(I,[1 3],'omitnan');
plot(ax,x,squeeze(h),'-d','Color',[0.80 0.12 0.20], ...
    'MarkerFaceColor',[0.80 0.12 0.20],'LineWidth',st.line_width);
plot(ax,x,squeeze(q),'--o','Color',[0.00 0.45 0.70], ...
    'MarkerFaceColor','w','LineWidth',st.line_width);
set(ax,'YScale','log'); xlabel(ax,xlab); ylabel(ax,'Count');
title(ax,ttl,'FontWeight','normal'); grid(ax,'on'); box(ax,'on');
legend(ax,{'Exact HVP calls','tCG iterations'},'Location','best', ...
    'Box','off','FontSize',st.legend_font_size);
end

function write_summary(d,meta)
T=[summary_rows('U',d.U_values,d.user_dimension,d.user_crlb, ...
    d.user_runtime_s,d.user_feasible,d.user_kkt,d.user_hvp_calls, ...
    d.user_tcg_iterations,meta); ...
   summary_rows('Nt',d.Nt_values,d.antenna_dimension,d.antenna_crlb, ...
    d.antenna_runtime_s,d.antenna_feasible,d.antenna_kkt, ...
    d.antenna_hvp_calls,d.antenna_tcg_iterations,meta)];
s=simtwc_config(); writetable(T,fullfile(s.figure_dir, ...
    'fig14_dimension_scaling_summary.csv'));
end

function T=summary_rows(sweep,values,dim,C,R,F,K,H,I,meta)
nA=numel(meta); nX=numel(values); n=nA*nX;
sweep_name=repmat(string(sweep),n,1); value=nan(n,1); variable_dimension=nan(n,1);
method=strings(n,1); mean_crlb=nan(n,1); mean_runtime_s=nan(n,1);
feasibility_pct=nan(n,1); kkt_success_pct=nan(n,1);
mean_hvp_calls=nan(n,1); mean_tcg_iterations=nan(n,1); row=0;
tol=simtwc_config().kkt_tol;
for a=1:nA
    for j=1:nX
        row=row+1; value(row)=values(j); variable_dimension(row)=dim(j);
        method(row)=string(meta(a).name);
        mean_crlb(row)=mean(C(a,j,:),'all','omitnan');
        mean_runtime_s(row)=mean(R(a,j,:),'all','omitnan');
        feasibility_pct(row)=100*mean(F(a,j,:),'all');
        kval=squeeze(K(a,j,:)); valid=isfinite(kval);
        if any(valid); kkt_success_pct(row)=100*mean(kval(valid)<=tol); end
        mean_hvp_calls(row)=mean(H(a,j,:),'all','omitnan');
        mean_tcg_iterations(row)=mean(I(a,j,:),'all','omitnan');
    end
end
T=table(sweep_name,value,variable_dimension,method,mean_crlb, ...
    mean_runtime_s,feasibility_pct,kkt_success_pct,mean_hvp_calls, ...
    mean_tcg_iterations);
end
