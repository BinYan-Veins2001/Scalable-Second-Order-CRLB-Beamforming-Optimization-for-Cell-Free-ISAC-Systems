function fig=paperfig_plot_efficiency(result)
%PAPERFIG_PLOT_EFFICIENCY Paper Fig. 6, consistent IEEE typography/markers.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root),'-begin');
if nargin<1
    a=load(fullfile(root,'results','paper_fig06','data.mat'),'result');
    result=a.result;
elseif ischar(result) || isstring(result)
    a=load(result,'result'); result=a.result;
end
st=simtwc_ieee_style('double');
fig=figure('Units','inches','Position',[1 1 st.width_in 2.65],'Color','w');
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
rows=result.rows; panels={'AP','User'}; colors={[.80 .12 .20],[.42 .28 .68]};
for p=1:2
    ax=nexttile(layout); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    select=strcmp({rows.panel},panels{p}); xx=unique([rows(select).x]);
    structured=zeros(size(xx)); dense=structured;
    for j=1:numel(xx)
        r=rows(select & [rows.x]==xx(j));
        structured(j)=median([r.spsc_median_s]); dense(j)=median([r.dense_median_s]);
    end
    h1=semilogy(ax,xx,structured,'-d','Color',colors{1},'MarkerFaceColor','w', ...
        'LineWidth',st.line_width,'MarkerSize',st.marker_size);
    h2=semilogy(ax,xx,dense,'--^','Color',colors{2},'MarkerFaceColor','w', ...
        'LineWidth',st.line_width,'MarkerSize',st.marker_size);
    set(ax,'YScale','log','XTick',xx);
    span=max(xx)-min(xx); padding=max(0.04*span,0.1);
    xlim(ax,[min(xx)-padding,max(xx)+padding]);
    ylabel(ax,'Curvature evaluation time (s)','Interpreter','latex');
    if p==1
        xlabel(ax,'Number of APs, $M$','Interpreter','latex');
        title(ax,sprintf('(a) $U=%d$, $N_t=N_r=%d$',result.options.fixed_U,result.options.Nt), ...
            'Interpreter','latex','FontWeight','normal','FontSize',st.axis_font_size);
    else
        xlabel(ax,'Number of users, $U$','Interpreter','latex');
        title(ax,sprintf('(b) $M=%d$, $N_t=N_r=%d$',result.options.fixed_M,result.options.Nt), ...
            'Interpreter','latex','FontWeight','normal','FontSize',st.axis_font_size);
    end
    legend(ax,[h1 h2],{'SPSC','Dense tangent EVD'},'Location','northwest', ...
        'Box','off','FontSize',st.legend_font_size);
end
if exist('paperfig_export','file')==2
    paperfig_export(fig,'fig06_structured_curvature_efficiency');
else
    simtwc_export_figure(fig,'fig06_structured_curvature_efficiency');
end
end
