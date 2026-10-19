function fig = simtwc_plot_fig15()
%SIMTWC_PLOT_FIG15 Two-row, three-column network scaling comparison.
%
% Top row: mean CRLB. Bottom row: mean end-to-end runtime. Columns vary M,
% Nt, and U, respectively. The same Fig. 8 method identity is used globally.

simtwc_setup(); d=simtwc_load_data(15); st=simtwc_ieee_style('double');
meta=simtwc_method_registry(d.methods);

% Enforce the current resource policy at plot time as well. Thus, lowering
% the cap hides archived SDR-SCA points without deleting their checkpoints.
apC=d.ap_crlb; apR=d.ap_runtime_s;
isdr=find(strcmp(d.methods,'sdr_sca'),1);
if ~isempty(isdr)
    apC(isdr,~d.ap_sdr_allowed,:)=NaN;
    apR(isdr,~d.ap_sdr_allowed,:)=NaN;
end
[apCmu,apCse]=simtwc_mean_se(apC,3);
[apRmu,apRse]=simtwc_mean_se(apR,3);
[ntCmu,ntCse]=simtwc_mean_se(d.antenna_crlb,3);
[ntRmu,ntRse]=simtwc_mean_se(d.antenna_runtime_s,3);
[uCmu,uCse]=simtwc_mean_se(d.user_crlb,3);
[uRmu,uRse]=simtwc_mean_se(d.user_runtime_s,3);

fig=figure('Units','inches','Position',[1 1 st.width_in 3.45]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); h=method_panel(ax,d.M_values,apCmu,apCse,meta,st, ...
    '(a)','AP number, $M$','','CRLB ($\mathrm{m}^2$)');
method_panel(nexttile(tl),d.Nt_values,ntCmu,ntCse,meta,st, ...
    '(b)','Transmit antennas, $N_t$','','');
method_panel(nexttile(tl),d.U_values,uCmu,uCse,meta,st, ...
    '(c)','User number, $U$','','');
method_panel(nexttile(tl),d.M_values,apRmu,apRse,meta,st, ...
    '(d)','','Number of APs, $M$','Runtime (s)');
method_panel(nexttile(tl),d.Nt_values,ntRmu,ntRse,meta,st, ...
    '(e)','','Transmit antennas per AP, $N_t$','');
axlast=nexttile(tl);
method_panel(axlast,d.U_values,uRmu,uRse,meta,st, ...
    '(f)','','Number of users, $U$','');
lg=legend(axlast,h,{meta.name},'Location','southoutside', ...
    'NumColumns',numel(meta),'Orientation','horizontal','Box','off', ...
    'FontSize',st.legend_font_size-0.2);
lg.Layout.Tile='south';
simtwc_export_figure(fig,'fig15_network_dimension_scaling');
write_summary(d,meta,apC,apR);
end

function h=method_panel(ax,x,mu,se,meta,st,tag,ttl,xlab,ylab)
hold(ax,'on'); h=gobjects(numel(meta),1);
for a=1:numel(meta)
    err=[]; if st.show_curve_se; err=se(a,:); end
    h(a)=simtwc_plot_method_curve(ax,x,mu(a,:),meta(a),st,err);
end
set(ax,'YScale','log','XTick',x);
if ~isempty(xlab); xlabel(ax,xlab,'Interpreter','latex'); end
if ~isempty(ylab); ylabel(ax,ylab,'Interpreter','latex'); end
if ~isempty(ttl); title(ax,ttl,'FontWeight','normal','Interpreter','latex'); end
text(ax,0.035,0.94,tag,'Units','normalized','VerticalAlignment','top', ...
    'BackgroundColor','w','Margin',1,'FontSize',st.axis_font_size);
grid(ax,'on'); box(ax,'on');
end

function write_summary(d,meta,apC,apR)
T=[summary_rows('M',d.M_values,apC,apR,d.ap_feasible,d.ap_kkt,meta); ...
   summary_rows('Nt',d.Nt_values,d.antenna_crlb,d.antenna_runtime_s, ...
        d.antenna_feasible,d.antenna_kkt,meta); ...
   summary_rows('U',d.U_values,d.user_crlb,d.user_runtime_s, ...
        d.user_feasible,d.user_kkt,meta)];
s=simtwc_config(); writetable(T,fullfile(s.figure_dir, ...
    'fig15_network_dimension_scaling_summary.csv'));
end

function T=summary_rows(sweep,values,C,R,F,K,meta)
nA=numel(meta); nX=numel(values); n=nA*nX;
sweep_name=repmat(string(sweep),n,1); value=nan(n,1); method=strings(n,1);
sample_count=zeros(n,1); mean_crlb=nan(n,1); mean_runtime_s=nan(n,1);
feasibility_pct=nan(n,1); kkt_success_pct=nan(n,1); row=0;
tol=simtwc_config().kkt_tol;
for a=1:nA
    for j=1:nX
        row=row+1; value(row)=values(j); method(row)=string(meta(a).name);
        c=squeeze(C(a,j,:)); r=squeeze(R(a,j,:));
        valid=isfinite(c)&isfinite(r); sample_count(row)=nnz(valid);
        if any(valid)
            mean_crlb(row)=mean(c(valid)); mean_runtime_s(row)=mean(r(valid));
            f=squeeze(F(a,j,:)); feasibility_pct(row)=100*mean(f(valid));
            k=squeeze(K(a,j,:));
            kvalid=valid&isfinite(k);
            if any(kvalid); kkt_success_pct(row)=100*mean(k(kvalid)<=tol); end
        end
    end
end
T=table(sweep_name,value,method,sample_count,mean_crlb,mean_runtime_s, ...
    feasibility_pct,kkt_success_pct);
end
