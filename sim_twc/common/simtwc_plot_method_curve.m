function h = simtwc_plot_method_curve(ax,x,y,meta,st,se,line_style,marker_indices)
%SIMTWC_PLOT_METHOD_CURVE Draw one method using the Fig. 8 visual identity.
%
% The main curve always carries the method color, line style, marker,
% line width, and marker size.  When standard errors are supplied, they are
% rendered first as lighter, thinner whiskers so that the main line retains
% the clean appearance of the AP-scaling figure.

if nargin<6; se=[]; end
if nargin<7 || isempty(line_style); line_style=meta.line; end
if nargin<8; marker_indices=[]; end

x=x(:).'; y=y(:).';
if ~isempty(se) && st.show_curve_se
    se=se(:).';
    error_color=0.58*meta.color+0.42*[1 1 1];
    errorbar(ax,x,y,se,'LineStyle','none','Marker','none', ...
        'Color',error_color,'LineWidth',st.error_line_width, ...
        'CapSize',st.error_cap_size,'HandleVisibility','off');
end

args={'Color',meta.color,'LineStyle',line_style, ...
    'LineWidth',st.line_width,'Marker',meta.marker, ...
    'MarkerSize',st.marker_size,'DisplayName',meta.name};
if ~isempty(marker_indices)
    args=[args,{'MarkerIndices',marker_indices}];
end
h=plot(ax,x,y,args{:});
end
