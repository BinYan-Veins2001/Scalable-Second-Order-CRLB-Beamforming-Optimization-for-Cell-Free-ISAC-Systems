function paperfig_export(fig,stem)
%PAPERFIG_EXPORT Editable MATLAB plus vector and raster publication outputs.
s=paperfig_setup(); st=simtwc_ieee_style();
set(fig,'Color','w','Renderer','painters');
axs=findall(fig,'Type','axes');
set(axs,'GridColor',st.grid_color,'MinorGridColor',st.grid_color, ...
    'XMinorGrid','off','YMinorGrid','off');
base=fullfile(s.figure_dir,stem);
savefig(fig,[base '.fig']);
exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
print(fig,[base '.eps'],'-depsc','-painters');
exportgraphics(fig,[base '.png'],'Resolution',200,'BackgroundColor','white');
end
