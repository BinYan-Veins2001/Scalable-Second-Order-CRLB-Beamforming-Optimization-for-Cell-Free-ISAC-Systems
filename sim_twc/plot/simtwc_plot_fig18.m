function fig = simtwc_plot_fig18()
%SIMTWC_PLOT_FIG18 Controlled second-order local-rate verification.

simtwc_setup(); d=simtwc_load_data(18); st=simtwc_ieee_style('single');
trials=d.trials(cellfun(@trial_ok,d.trials));
if isempty(trials)
    error('Fig. 18 has no theorem-admissible local-rate trials.');
end

blue=[0.00 0.36 0.68]; orange=[0.88 0.33 0.10];
teal=[0.00 0.55 0.50]; gray=[0.38 0.38 0.38];
[r0,r1]=pooled_rate_data(trials);
C=r1./r0.^2; Cmed=median(C,'omitnan');

fig=figure('Units','inches','Position',[1 1 st.width_in 5.35]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

%% (a) Quadratic one-step stationarity map
ax=nexttile(tl,1); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
[x,order]=sort(r0); y=r1(order);
hMeasured=loglog(ax,x,y,'o-','Color',blue,'MarkerFaceColor','w', ...
    'MarkerEdgeColor',blue,'MarkerSize',5.2,'LineWidth',1.35);
xx=logspace(log10(min(x)/1.18),log10(max(x)*1.18),240);
hQuadratic=loglog(ax,xx,Cmed*xx.^2,'-','Color',orange,'LineWidth',2.0);
xmid=sqrt(min(x)*max(x)); ymid=Cmed*xmid^2;
hLinear=loglog(ax,xx,(ymid/xmid)*xx,'--','Color',gray,'LineWidth',1.2);
set(ax,'XScale','log','YScale','log');
xlim(ax,[xx(1) xx(end)]);
xlabel(ax,'$r_k=\|\mathbf{g}^{\mathrm{tan}}_k\|_2$', ...
    'Interpreter','latex');
ylabel(ax,'$r_{k+1}=\|\mathbf{g}^{\mathrm{tan}}_{k+1}\|_2$', ...
    'Interpreter','latex');
title(ax,'(a) Local stationarity contraction','FontWeight','normal');
legend(ax,[hMeasured hQuadratic hLinear], ...
    {'COSTAR (exact HVP)','$C r_k^2$','$C_{\mathrm{lin}}r_k$'}, ...
    'Interpreter','latex','Location','northwest','Box','off', ...
    'FontSize',st.legend_font_size);
fitText={sprintf('$p=%.3f$',d.pooled_slope), ...
    sprintf('$95\\%%\\;\\mathrm{CI}=[%.3f,\\,%.3f]$', ...
    d.pooled_slope_ci95(1),d.pooled_slope_ci95(2))};
text(ax,0.97,0.07,fitText,'Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'Color',blue,'FontWeight','bold','BackgroundColor','w','Margin',2, ...
    'Interpreter','latex');

%% (b) Stability of the quadratic constant
ax=nexttile(tl,2); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
xnorm=r0/max(r0); [xnorm,order]=sort(xnorm); Cnorm=C(order)/Cmed;
hRatio=plot(ax,xnorm,Cnorm,'d-','Color',teal, ...
    'MarkerFaceColor',[0.78 0.94 0.91],'MarkerEdgeColor',teal, ...
    'MarkerSize',5.0,'LineWidth',1.45);
hUnity=yline(ax,1,'--','Color',gray,'LineWidth',1.2);
pad=max(0.06,0.18*(max(Cnorm)-min(Cnorm)));
ylim(ax,[min(Cnorm)-pad max(Cnorm)+pad]); xlim(ax,[min(xnorm)*0.96 1.02]);
xlabel(ax,'$r_k/r_{\max}$','Interpreter','latex');
ylabel(ax,'$C_k/\mathrm{median}(C_k)$','Interpreter','latex');
title(ax,'(b) Quadratic-ratio stability','FontWeight','normal');
legend(ax,[hRatio hUnity], ...
    {'$C_k/\mathrm{median}(C_k)$','$1$'}, ...
    'Interpreter','latex','Location','northeast','Box','off', ...
    'FontSize',st.legend_font_size);
text(ax,0.04,0.94,sprintf('$\\mathrm{CV}=%.1f\\%%$', ...
    100*std(C)/mean(C)), ...
    'Units','normalized','VerticalAlignment','top','Color',teal, ...
    'FontWeight','bold','BackgroundColor','w','Margin',2, ...
    'Interpreter','latex');
text(ax,0.04,0.08,'$C_k=r_{k+1}/r_k^2$', ...
    'Units','normalized','VerticalAlignment','bottom','Color',gray, ...
    'BackgroundColor','w','Margin',2,'Interpreter','latex');

simtwc_export_figure(fig,'fig18_quadratic_convergence');
write_summaries(d,trials);
end

function [r0,r1]=pooled_rate_data(trials)
r0=[]; r1=[];
for seed_idx=1:numel(trials)
    q=trials{seed_idx}; mask=q.fit_mask;
    r0=[r0 q.tangent_residual_before(mask)]; %#ok<AGROW>
    r1=[r1 q.tangent_residual_after(mask)]; %#ok<AGROW>
end
end

function write_summaries(d,trials)
seed=[]; r0=[]; r1=[]; quadraticRatio=[]; forcingRatio=[];
accepted=[]; feasible=[]; stable=[]; interior=[]; exactHvp=[]; used=[];
for seed_idx=1:numel(trials)
    q=trials{seed_idx}; n=numel(q.fit_mask);
    seed=[seed;repmat(q.seed,n,1)]; %#ok<AGROW>
    r0=[r0;q.tangent_residual_before(:)]; %#ok<AGROW>
    r1=[r1;q.tangent_residual_after(:)]; %#ok<AGROW>
    quadraticRatio=[quadraticRatio; ...
        (q.tangent_residual_after./max(q.tangent_residual_before.^2,realmin)).']; %#ok<AGROW>
    forcingRatio=[forcingRatio;q.forcing_ratio(:)]; %#ok<AGROW>
    accepted=[accepted;q.accepted(:)]; feasible=[feasible;q.feasible(:)]; %#ok<AGROW>
    stable=[stable;q.active_set_stable(:)]; interior=[interior;q.interior_step(:)]; %#ok<AGROW>
    exactHvp=[exactHvp;q.exact_hvp(:)]; used=[used;q.fit_mask(:)]; %#ok<AGROW>
end
sampleTable=table(seed,r0,r1,quadraticRatio,forcingRatio,accepted, ...
    feasible,stable,interior,exactHvp,used,'VariableNames', ...
    {'seed','r_k','r_k_plus_1','quadratic_ratio','forcing_ratio', ...
    'accepted','feasible','active_set_stable','interior_step', ...
    'exact_hvp','used_in_fit'});

seed=cellfun(@(q)q.seed,trials).';
slope=cellfun(@(q)q.slope,trials).';
slopeCiLow=cellfun(@(q)q.slope_ci95(1),trials).';
slopeCiHigh=cellfun(@(q)q.slope_ci95(2),trials).';
quadraticConstant=cellfun(@(q)q.quadratic_constant,trials).';
quadraticCv=cellfun(@(q)q.quadratic_cv,trials).';
modelAdvantage=cellfun(@(q)q.model_advantage,trials).';
qualifiedRate=cellfun(@(q)q.qualified_rate,trials).';
seedTable=table(seed,slope,slopeCiLow,slopeCiHigh,quadraticConstant, ...
    quadraticCv,modelAdvantage,qualifiedRate);

s=simtwc_config();
writetable(sampleTable,fullfile(s.figure_dir, ...
    'fig18_quadratic_convergence_summary.csv'));
writetable(seedTable,fullfile(s.figure_dir, ...
    'fig18_quadratic_convergence_per_seed.csv'));

audit=trials{1}.hvp_audit;
hvpTable=table(audit.fd_steps(:),audit.fd_relative_error(:), ...
    'VariableNames',{'fd_displacement','relative_hvp_error'});
writetable(hvpTable,fullfile(s.figure_dir, ...
    'fig18_exact_hvp_consistency.csv'));

fid=fopen(fullfile(s.figure_dir,'fig18_quadratic_convergence_fit.txt'),'w');
if fid>=0
    fprintf(fid,'successful_realizations=%d\n',d.successful_seed_count);
    fprintf(fid,'pooled_order=%.10f\n',d.pooled_slope);
    fprintf(fid,'pooled_order_ci95=[%.10f, %.10f]\n', ...
        d.pooled_slope_ci95(1),d.pooled_slope_ci95(2));
    fprintf(fid,'quadratic_ratio_cv=%.10f\n', ...
        std(quadraticRatio)/mean(quadraticRatio));
    fprintf(fid,'quadratic_vs_linear_rmse_ratio=%.10f\n', ...
        median(modelAdvantage));
    fprintf(fid,'qualified_steps=%d/%d\n',nnz(used),numel(used));
    fprintf(fid,'minimum_hvp_discrepancy=%.10e\n', ...
        min(audit.fd_relative_error));
    fprintf(fid,'self_adjoint_error=%.10e\n',audit.self_adjoint_error);
    fclose(fid);
end
end

function yes=trial_ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok')&& ...
    isfield(q,'fit_mask')&&nnz(q.fit_mask)>=3;
end
