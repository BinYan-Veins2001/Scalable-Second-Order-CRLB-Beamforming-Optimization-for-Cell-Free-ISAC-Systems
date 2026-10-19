function fig=paperfig_plot03(archive)
%PAPERFIG_PLOT03 Feasible-endpoint horizontal Ritz-curvature bars.

s=paperfig_setup(); if nargin<1; archive=fullfile(s.results_dir,'paper_fig03','data.mat'); end
width_mode='single';
single_column=strcmpi(width_mode,'single');
if ~single_column && ~strcmpi(width_mode,'double')
    error('Fig. 3 width mode must be single or double.');
end
d=load(archive); st=simtwc_ieee_style(width_mode);
d.feasible=false(size(d.runs)); d.lambda_min_horizontal=nan(size(d.runs)); d.kkt_residual=nan(size(d.runs));
for j=1:numel(d.runs)
    if isempty(d.runs{j}); continue; end
    d.feasible(j)=d.runs{j}.feasible;
    d.lambda_min_horizontal(j)=d.audits{j}.lambda_min_horizontal;
    d.kkt_residual(j)=d.audits{j}.kkt;
end
meta=simtwc_method_registry(d.methods);
nA=numel(meta); transition=1e-6;
median_curvature=nan(nA,1); feasible_count=zeros(nA,1);
median_kkt=nan(nA,1);
for a=1:nA
    valid=d.feasible(a,:) & isfinite(d.lambda_min_horizontal(a,:));
    feasible_count(a)=nnz(valid);
    if any(valid)
        median_curvature(a)=median(d.lambda_min_horizontal(a,valid));
        median_kkt(a)=median(d.kkt_residual(a,valid));
    end
end

if single_column; height=2.65; else; height=2.8; end
fig=figure('Units','inches','Position',[1 1 st.width_in height]);
ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
bar_values=simtwc_curvature_axis(median_curvature,transition);
b=bar(ax,1:nA,bar_values,0.65,'FaceColor','flat', ...
    'EdgeColor',st.bar_edge_color,'LineWidth',st.bar_edge_width);
for a=1:nA
    b.CData(a,:)=meta(a).color;
end
b.FaceAlpha=0.90;
yline(ax,0,'-','Color',[0.28 0.28 0.28], ...
    'LineWidth',st.reference_line_width);
finite_bars=bar_values(isfinite(bar_values));
if isempty(finite_bars)
    error('Fig. 3 has no feasible endpoint with finite curvature.');
end
span=max(0.6,max(finite_bars)-min(finite_bars));
ylim(ax,[min(0,min(finite_bars))-0.10*span, ...
    max(0,max(finite_bars))+0.10*span]);
simtwc_curvature_axis([],transition,ax);
set(ax,'XTick',1:nA,'XTickLabel',{meta.name},'XLim',[0.45 nA+0.55]);
if single_column
    set(ax,'Position',[0.18 0.28 0.78 0.68], ...
        'XTickLabelRotation',37,'FontSize',st.axis_font_size);
    ylabel(ax,'$\hat\lambda_{\min}^{\mathrm H}$','Interpreter','latex');
    basename='fig03_endpoint_curvature';
else
    set(ax,'XTickLabelRotation',0);
    ylabel(ax,'$\hat\lambda_{\min}^{\mathrm H}$ (signed log)', ...
        'Interpreter','latex');
    xlabel(ax,'Algorithm','Interpreter','latex');
    basename='fig03_endpoint_curvature';
end
paperfig_export(fig,basename);

s=paperfig_setup();
method=string({meta.name}).';
feasible_percent=100*feasible_count/numel(d.seeds);
T=table(method,feasible_count,feasible_percent, ...
    median_curvature,median_kkt);
writetable(T,fullfile(s.results_dir,'paper_fig03', ...
    'summary.csv'));
end
