function fig = paperfig_plot04(archive)
%PAPERFIG_PLOT04 Single-column joint-stationarity CDF for Theorem 3.
s=paperfig_setup();
if nargin<1; archive=fullfile(s.results_dir,'paper_fig04','data.mat'); end
d=load(archive); o=d.options; o.regimes=d.regimes;
assert(iscell(d.trials) && ~isempty(d.trials),'Fig. 4 trials are missing.');
assert(size(d.trials,1)==numel(o.regimes),'Fig. 4 regime count is inconsistent.');
st=simtwc_ieee_style('single');
colors={[0.00 0.45 0.70],[0.90 0.48 0.08],[0.42 0.28 0.68], ...
    [0.00 0.58 0.43],[0.80 0.12 0.20]};
styles={'-o','--s','-.d',':^','-v'};

% Noncompleted or infeasible trials remain in the denominator at +Inf.
q=inf(size(d.trials));
for ir=1:size(d.trials,1)
    for is=1:size(d.trials,2)
        trial=d.trials{ir,is};
        if isempty(trial) || ~strcmp(trial.status,'completed'); continue; end
        if trial.audit.primal<=trial.solution.terminal_options.feas_tol && ...
                isfinite(trial.audit.joint_ratio) && trial.audit.joint_ratio>=0
            q(ir,is)=trial.audit.joint_ratio;
        end
    end
end
positive=q(isfinite(q) & q>0);
if isempty(positive); xmin=1e-3; else; xmin=max(realmin,min(positive)/2); end
xmin=min(xmin,0.1);
xmax=max(2,max([1,positive(:).'])*1.5);

fig=figure('Units','inches','Position',[1 1 st.width_in 2.6],'Color','w');
ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
handles=gobjects(size(q,1),1); labels=cell(size(q,1),1);
for ir=1:size(q,1)
    samples=q(ir,isfinite(q(ir,:)));
    xx=unique([xmin,samples(samples>xmin),1,xmax]);
    yy=arrayfun(@(x) 100*sum(q(ir,:)<=x)/size(q,2),xx);
    handles(ir)=stairs(ax,xx,yy,styles{1+mod(ir-1,numel(styles))}, ...
        'Color',colors{1+mod(ir-1,numel(colors))},'MarkerFaceColor','w', ...
        'LineWidth',st.line_width,'MarkerSize',4);
    labels{ir}=o.regimes(ir).name;
end
threshold=xline(ax,1,':','Color',[0.25 0.25 0.25], ...
    'LineWidth',st.reference_line_width);
set(ax,'XScale','log','FontName',st.font_name,'FontSize',st.axis_font_size, ...
    'LineWidth',0.7,'TickLabelInterpreter','latex','Layer','top');
xlim(ax,[xmin*0.8,xmax*1.25]); ylim(ax,[-3 104]);
exponents=ceil(log10(xmin)):floor(log10(xmax));
if numel(exponents)>7; exponents=exponents(1:2:end); end
set(ax,'XTick',unique([10.^exponents,1]),'YTick',0:20:100);
xlabel(ax,'Joint stopping ratio $q_{\mathrm{stop}}$','Interpreter','latex');
ylabel(ax,'Empirical CDF (\%)','Interpreter','latex');
legend(ax,[handles;threshold],[labels;{'$q_{\mathrm{stop}}=1$'}], ...
    'Interpreter','latex','Location','northwest','Box','off', ...
    'FontSize',st.legend_font_size);
set(fig,'UserData',struct('joint_ratio',q,'attainment_rate',mean(q<=1,2), ...
    'epsilon_KKT',o.epsilon_KKT,'epsilon_H',o.epsilon_H));
drawnow;
paperfig_export(fig,'fig04_joint_stationarity');
fprintf('[Paper Fig4] Single-panel PDF/EPS/PNG exported to %s\n',s.figure_dir);
for ir=1:size(q,1)
    fprintf('[Paper Fig4] %s: %d/%d jointly attained.\n', ...
        o.regimes(ir).name,sum(q(ir,:)<=1),size(q,2));
end
end
