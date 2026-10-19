function simtwc_export_figure(fig,basename)
%SIMTWC_EXPORT_FIGURE Save publication vector PDF and EPS files.

s=simtwc_config();
if ~exist(s.figure_dir,'dir'); mkdir(s.figure_dir); end
base=fullfile(s.figure_dir,basename);
set(fig,'Color','w','Renderer','painters');
st=simtwc_ieee_style();
axs=findall(fig,'Type','axes');
set(axs,'GridColor',st.grid_color,'MinorGridColor',st.grid_color, ...
    'XMinorGrid','off','YMinorGrid','off');
exportgraphics(fig,[base '.pdf'],'ContentType','vector', ...
    'BackgroundColor','white');
print(fig,[base '.eps'],'-depsc','-painters');
exportgraphics(fig,[base '.png'],'Resolution',200,'BackgroundColor','white');
end
