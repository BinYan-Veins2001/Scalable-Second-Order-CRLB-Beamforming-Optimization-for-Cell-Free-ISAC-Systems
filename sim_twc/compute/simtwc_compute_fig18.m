function archive = simtwc_compute_fig18()
%SIMTWC_COMPUTE_FIG18 Controlled local quadratic-rate audit.
%
% A high-accuracy COSTAR endpoint defines the identified local geometry.
% Feasible perturbations are restored to the same constraint surface and
% advanced by one exact-HVP AS-RTR step with the theorem's O(r_k^2)
% inexact-Newton forcing rule. Only accepted, feasible, active-set-stable,
% interior steps enter the rate fit. A separate directional audit checks
% the exact HVP against central differences and self-adjointness.

s=simtwc_setup(); spec=s.fig18;
data=simtwc_initialize_data(18,spec,{'costar'});
path=simtwc_archive_path(18);
seeds=spec.seeds(:).';
if ~isfield(data,'trials') || numel(data.trials)~=numel(seeds)
    data.trials=cell(1,numel(seeds));
end

for seed_idx=1:numel(seeds)
    if trial_ok(data.trials{seed_idx}); continue; end
    seed=seeds(seed_idx);
    fprintf('[Fig18] Seed %d (%d/%d): reference and local audit\n', ...
        seed,seed_idx,numel(seeds));
    try
        data.trials{seed_idx}=run_trial(spec,s,seed);
    catch ME
        data.trials{seed_idx}=struct('status','failed','seed',seed, ...
            'message',ME.message,'identifier',ME.identifier);
        fprintf('    FAILED seed %d: %s\n',seed,ME.message);
    end
    simtwc_save_atomic(path,data);
end

[data.pooled_slope,data.pooled_intercept_by_seed, ...
    data.pooled_slope_se,data.pooled_slope_ci95,data.pooled_seed_ids] = ...
    fixed_effect_order(data.trials);
data.seeds=seeds;
data.successful_seed_count=nnz(cellfun(@trial_ok,data.trials));
data.attempted_seed_count=nnz(cellfun(@is_attempted,data.trials));
data.complete=data.attempted_seed_count==numel(seeds);
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data);
archive=path;
end

function trial=run_trial(spec,s,seed)
p=simtwc_make_problem(spec,seed);
ropts=reference_options(spec,seed);
reference=simtwc_try_run('costar',p,true,struct('rtr_opts',ropts));
if ~record_ok(reference)
    error('High-accuracy COSTAR reference failed: %s',reference.message);
end
reference=tighten_reference(p,reference,spec);

Wstar=reference.Wmat;
[Wbase,base_error,base_index]=feasible_side_base( ...
    p,reference,Wstar,spec);
Wbase=phase_align(Wbase,Wstar);
[direction,active_reference]=identified_tangent_direction( ...
    p,Wstar,Wbase,spec,seed);
targets=sort(real(spec.local_test_errors(:)),'descend');

opts=reference.solution.options;
opts.phase2_variant='core_kkt_triggered';
opts.max_iter=1;
opts.delta0=spec.one_step_delta;
opts.delta_max=spec.one_step_delta;
opts.curvature_policy='none';
opts.stationarity_tol=0;
opts.grad_floor=0;
opts.min_actual_decrease=0;
opts.pred_tol=0;
opts.hessian_model='exact_lagrangian';
opts.sinr_active_tol=spec.sinr_active_tol;
opts.power_active_tol=spec.power_active_tol;
opts.cg_max_iter=spec.local_cg_max_iter;
opts.cg_tol=spec.local_cg_tol;
opts.cg_forcing_power=2;
opts.cg_forcing_constant=spec.local_forcing_constant;
opts.cg_absolute_floor=spec.local_cg_absolute_floor;
opts.project_phase_directions=true;
opts.enable_tcg_negative_curvature=false;
opts.sinr_retraction_margin=0;
opts.sinr_retraction_max_iter=8;
opts.sinr_retraction_regularization=1e-14;

one_step=cell(1,numel(targets));
for point_idx=1:numel(targets)
    W0=identified_manifold_start( ...
        p,Wstar,direction,active_reference,targets(point_idx),spec);
    q0=simtwc_unified_kkt(p,W0);
    if q0.r_c>s.feasibility_tol
        one_step{point_idx}=struct('status','failed', ...
            'message',sprintf('Local start is infeasible: r_c=%.3e.',q0.r_c));
        continue;
    end
    gtan0=local_tangent_gradient_norm(p,W0,opts);
    fprintf('    point %d/%d, e0=%.3e\n', ...
        point_idx,numel(targets),targets(point_idx));
    try
        [sol,hist]=solve_active_rtr_phase2( ...
            p.cfg,p.channels,p.smat,matrix_to_wcell(p.cfg,W0),opts);
        W1=vectorize_wcell_to_matrix(p.cfg,sol.wcell);
        q1=simtwc_unified_kkt(p,W1);
        accepted=~isempty(hist.accepted) && logical(hist.accepted(1));
        interior=~isempty(hist.step_norm) && ~isempty(hist.delta) && ...
            hist.step_norm(1)<0.99*hist.delta(1);
        stable=size(hist.active_membership,2)>=2 && ...
            isequal(hist.active_membership(:,1), ...
            hist.active_membership(:,end));
        exact_hvp=isfield(hist,'work') && ...
            hist.work.exact_hvp_calls>0 && hist.work.fd_hvp_calls==0;
        linear_residual=first_or_nan(hist.tcg_final_res_norm);
        one_step{point_idx}=struct('status','ok', ...
            'target_error',targets(point_idx), ...
            'e0',quotient_distance(W0,Wstar), ...
            'e1',quotient_distance(W1,Wstar), ...
            'accepted',accepted,'feasible',q1.r_c<=s.feasibility_tol, ...
            'active_set_stable',stable,'interior_step',interior, ...
            'exact_hvp',exact_hvp,'gtan_before',gtan0, ...
            'linear_residual',linear_residual, ...
            'forcing_ratio',linear_residual/max(gtan0^2,realmin), ...
            'step_norm',first_or_nan(hist.step_norm), ...
            'trust_radius',first_or_nan(hist.delta), ...
            'rho',first_or_nan(hist.ratio), ...
            'kkt_before',q0.r_kkt,'kkt_after',q1.r_kkt, ...
            'model_stationarity_before',hist.stationarity(1), ...
            'model_stationarity_after',hist.stationarity(end), ...
            'tangent_residual_before',hist.tangent_grad_norm(1), ...
            'tangent_residual_after',hist.tangent_grad_norm(end));
    catch ME
        one_step{point_idx}=struct('status','failed', ...
            'message',ME.message,'identifier',ME.identifier);
    end
end

trial=struct('status','ok','seed',seed,'reference',reference, ...
    'reference_Wmat',Wstar,'base_Wmat',Wbase,'base_error',base_error, ...
    'base_index',base_index,'target_errors',targets, ...
    'one_step',{one_step});
trial.hvp_audit=exact_hvp_audit(p,Wstar,spec,seed);
trial=finalize_trial(trial);
if nnz(trial.fit_mask)<3
    trial.status='failed';
    trial.message='Fewer than three theorem-admissible local steps.';
end
end

function reference=tighten_reference(p,reference,spec)
opts=reference.solution.options;
opts.phase2_variant='core_kkt_triggered';
opts.max_iter=25;
opts.stationarity_tol=1e-13;
opts.grad_floor=0;
opts.curvature_policy='none';
opts.min_actual_decrease=0;
opts.pred_tol=0;
opts.delta0=spec.one_step_delta;
opts.delta_max=spec.one_step_delta;
opts.hessian_model='exact_lagrangian';
opts.cg_max_iter=spec.local_cg_max_iter;
opts.cg_forcing_power=2;
opts.cg_forcing_constant=spec.local_forcing_constant;
opts.cg_absolute_floor=spec.local_cg_absolute_floor;
opts.project_phase_directions=true;
opts.enable_tcg_negative_curvature=false;
opts.sinr_retraction_margin=0;
opts.sinr_retraction_max_iter=8;
opts.sinr_retraction_regularization=1e-14;
try
    [sol,hist]=solve_active_rtr_phase2( ...
        p.cfg,p.channels,p.smat,matrix_to_wcell(p.cfg,reference.Wmat),opts);
    W=vectorize_wcell_to_matrix(p.cfg,sol.wcell);
    q=simtwc_unified_kkt(p,W);
    if q.r_kkt<reference.kkt.r_kkt
        reference.entry_Wmat=reference.Wmat;
        reference.refinement_history=hist;
        reference.Wmat=W;
        reference.wcell=sol.wcell;
        reference.metrics=sol.metrics;
        reference.kkt=q;
        reference.solution=sol;
    end
catch
    % Retain the independently computed COSTAR endpoint if the optional
    % high-accuracy polishing pass cannot improve it.
end
end

function trial=finalize_trial(trial)
n=numel(trial.one_step);
fields={'e0','e1','gtan_before','linear_residual','forcing_ratio', ...
    'step_norm','trust_radius','rho','kkt_before','kkt_after', ...
    'model_stationarity_before','model_stationarity_after', ...
    'tangent_residual_before','tangent_residual_after'};
for field_idx=1:numel(fields)
    trial.(fields{field_idx})=nan(1,n);
end
flags={'accepted','feasible','active_set_stable','interior_step','exact_hvp'};
for flag_idx=1:numel(flags)
    trial.(flags{flag_idx})=false(1,n);
end
for point_idx=1:n
    q=trial.one_step{point_idx};
    if ~record_ok(q); continue; end
    for field_idx=1:numel(fields)
        trial.(fields{field_idx})(point_idx)=q.(fields{field_idx});
    end
    for flag_idx=1:numel(flags)
        trial.(flags{flag_idx})(point_idx)=logical(q.(flags{flag_idx}));
    end
end
trial.fit_mask=trial.accepted & trial.feasible & ...
    trial.active_set_stable & trial.interior_step & trial.exact_hvp & ...
    isfinite(trial.tangent_residual_before) & ...
    isfinite(trial.tangent_residual_after) & ...
    trial.tangent_residual_before>0 & trial.tangent_residual_after>0 & ...
    trial.tangent_residual_after<trial.tangent_residual_before;

trial.slope=NaN; trial.intercept=NaN; trial.slope_se=NaN;
trial.slope_ci95=[NaN NaN]; trial.quadratic_constant=NaN;
trial.quadratic_cv=NaN; trial.rmse_linear=NaN;
trial.rmse_quadratic=NaN; trial.model_advantage=NaN;
if nnz(trial.fit_mask)>=2
    r0=trial.tangent_residual_before(trial.fit_mask);
    r1=trial.tangent_residual_after(trial.fit_mask);
    [trial.slope,trial.intercept,trial.slope_se]=loglog_fit(r0,r1);
    trial.slope_ci95=trial.slope+[-1 1]*1.96*trial.slope_se;
    C=r1./r0.^2;
    trial.quadratic_constant=median(C,'omitnan');
    trial.quadratic_cv=std(C)/max(mean(C),realmin);
    trial.rmse_linear=fixed_order_rmse(r0,r1,1);
    trial.rmse_quadratic=fixed_order_rmse(r0,r1,2);
    trial.model_advantage=trial.rmse_linear/ ...
        max(trial.rmse_quadratic,realmin);
end
trial.qualified_rate=mean(trial.fit_mask);
end

function [Wbase,err,idx]=feasible_side_base(p,reference,Wstar,spec)
Z=reference.history.z_snapshots;
n=size(Z,2); errors=nan(1,n); feasible=false(1,n);
for state_idx=1:n
    W=real_vector_to_complex_matrix( ...
        Z(:,state_idx),p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
    errors(state_idx)=quotient_distance(W,Wstar);
    q=simtwc_unified_kkt(p,W);
    feasible(state_idx)=q.r_c<=simtwc_config().feasibility_tol;
end
candidates=find(feasible & errors>1.05*max(spec.local_test_errors) & ...
    errors<=spec.local_error_max);
if isempty(candidates)
    candidates=find(feasible & errors>1.05*max(spec.local_test_errors));
end
if isempty(candidates)
    error('No feasible trajectory side is available for the local test.');
end
% The closest feasible point supplies a trajectory-consistent orientation;
% the actual local starts are reconstructed on the identified constraint
% surface below.
[err,pos]=min(errors(candidates)); idx=candidates(pos);
Wbase=real_vector_to_complex_matrix( ...
    Z(:,idx),p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
end

function [direction,active]=identified_tangent_direction( ...
    p,Wstar,Wbase,spec,seed)
zstar=complex_matrix_to_real_vector(Wstar);
[c,G]=real_constraints(p,Wstar);
active=select_active(c,p.cfg.sys.U,spec);
Q=horizontal_normal_basis(Wstar,G(:,active));
raw=complex_matrix_to_real_vector(Wbase-Wstar);
direction=raw-Q*(Q'*raw);
if norm(direction)<=1e-10*max(1,norm(raw))
    stream=RandStream('mt19937ar','Seed',seed+18000);
    raw=randn(stream,numel(zstar),1);
    direction=raw-Q*(Q'*raw);
end
% Remove components in the locally flat solution set.  Applying the
% exact reduced Lagrangian Hessian maps the trajectory orientation into
% the identifiable horizontal range used by the local-rate theorem.
g=real_objective_gradient(p,Wstar);
lambda=zeros(numel(c),1);
if any(active); lambda(active)=lsqnonneg(G(:,active),-g); end
V=real_vector_to_complex_matrix(direction,size(Wstar,1),size(Wstar,2));
Hdirection=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
    p.cfg,p.channels,p.smat,Wstar,V,lambda));
direction=Hdirection-Q*(Q'*Hdirection);
direction=direction/max(norm(direction),realmin);
end

function W0=identified_manifold_start( ...
    p,Wstar,direction,active,target_error,spec)
zstar=complex_matrix_to_real_vector(Wstar);
scale=target_error*max(1,norm(Wstar,'fro'));
z=zstar+scale*direction;
target=-spec.local_manifold_margin*ones(nnz(active),1);
for correction_iter=1:spec.local_manifold_corrections
    W=real_vector_to_complex_matrix(z,size(Wstar,1),size(Wstar,2));
    [c,G]=real_constraints(p,W);
    A=G(:,active);
    residual=c(active)-target;
    if norm(residual,inf)<=spec.local_manifold_tolerance; break; end
    correction=lsqminnorm(A',-residual,1e-10);
    correction_cap=max(spec.local_correction_cap_fraction*scale,1e-10);
    if norm(correction)>correction_cap
        correction=correction*(correction_cap/norm(correction));
    end
    base_norm=norm(residual);
    accepted=false;
    for backtrack=0:spec.local_correction_backtracks
        alpha=2^(-backtrack);
        ztrial=z+alpha*correction;
        Wtrial=real_vector_to_complex_matrix( ...
            ztrial,size(Wstar,1),size(Wstar,2));
        ctrial=compute_constraints_grad( ...
            p.cfg,p.channels,p.smat,Wtrial);
        if norm(ctrial(active)-target)<base_norm
            z=ztrial; accepted=true; break;
        end
    end
    if ~accepted; break; end
end
W0=real_vector_to_complex_matrix(z,size(Wstar,1),size(Wstar,2));
W0=phase_align(W0,Wstar);
end

function audit=exact_hvp_audit(p,Wstar,spec,seed)
z=complex_matrix_to_real_vector(Wstar);
[c,G]=real_constraints(p,Wstar);
active=select_active(c,p.cfg.sys.U,spec);
g=real_objective_gradient(p,Wstar);
lambda=zeros(numel(c),1);
if any(active); lambda(active)=lsqnonneg(G(:,active),-g); end
Q=horizontal_normal_basis(Wstar,G(:,active));
stream=RandStream('mt19937ar','Seed',seed+19000);
raw=randn(stream,numel(z),1); v=raw-Q*(Q'*raw);
v=v/max(norm(v),realmin);
V=real_vector_to_complex_matrix(v,size(Wstar,1),size(Wstar,2));
Hexact=exact_lagrangian_hess_vec( ...
    p.cfg,p.channels,p.smat,Wstar,V,lambda);
hexact=complex_matrix_to_real_vector(Hexact);
steps=spec.hvp_fd_steps(:).'; errors=nan(size(steps));
for step_idx=1:numel(steps)
    h=steps(step_idx);
    gp=lagrangian_gradient(p,z+h*v,lambda);
    gm=lagrangian_gradient(p,z-h*v,lambda);
    hfd=(gp-gm)/(2*h);
    errors(step_idx)=norm(hfd-hexact)/max(norm(hexact),realmin);
end
raw2=randn(stream,numel(z),1); u=raw2-Q*(Q'*raw2);
u=u/max(norm(u),realmin);
U=real_vector_to_complex_matrix(u,size(Wstar,1),size(Wstar,2));
Hu=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
    p.cfg,p.channels,p.smat,Wstar,U,lambda));
symmetry=abs(u'*hexact-v'*Hu)/ ...
    max([realmin,abs(u'*hexact),abs(v'*Hu)]);
audit=struct('fd_steps',steps,'fd_relative_error',errors, ...
    'self_adjoint_error',symmetry,'active_count',nnz(active), ...
    'lambda',lambda,'exact_hvp',hexact);
end

function value=lagrangian_gradient(p,z,lambda)
W=real_vector_to_complex_matrix(z,p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
[~,gcomplex]=compute_crlb_objective_grad(p.cfg,p.smat,W);
[~,Gcomplex]=compute_constraints_grad(p.cfg,p.channels,p.smat,W);
value=complex_matrix_to_real_vector(gcomplex);
for constraint_idx=1:numel(lambda)
    value=value+lambda(constraint_idx)* ...
        complex_matrix_to_real_vector(Gcomplex(:,:,constraint_idx));
end
end

function g=real_objective_gradient(p,W)
[~,gcomplex]=compute_crlb_objective_grad(p.cfg,p.smat,W);
g=complex_matrix_to_real_vector(gcomplex);
end

function [c,G]=real_constraints(p,W)
[c,Gcomplex]=compute_constraints_grad(p.cfg,p.channels,p.smat,W);
dim=2*numel(W); G=zeros(dim,numel(c));
for constraint_idx=1:numel(c)
    G(:,constraint_idx)= ...
        complex_matrix_to_real_vector(Gcomplex(:,:,constraint_idx));
end
end

function active=select_active(c,U,spec)
active=false(size(c));
active(1:U)=c(1:U)>=-spec.sinr_active_tol;
active(U+1:end)=c(U+1:end)>=-spec.power_active_tol;
end

function Q=horizontal_normal_basis(W,A)
[N,U]=size(W); phase=zeros(2*N*U,U);
for user_idx=1:U
    q=zeros(N,U); q(:,user_idx)=1i*W(:,user_idx);
    phase(:,user_idx)=complex_matrix_to_real_vector(q);
end
Q=orth([A phase],1e-10);
end

function opts=reference_options(spec,seed)
opts=struct('phase2_variant','core_kkt_triggered', ...
    'max_iter',spec.reference_max_iter, ...
    'stationarity_tol',spec.reference_stationarity_tol, ...
    'curvature_tol',spec.reference_curvature_tol, ...
    'cg_tol',spec.reference_cg_tol, ...
    'cg_max_iter',spec.reference_cg_max_iter, ...
    'lanczos_max_iter',spec.reference_lanczos_max_iter, ...
    'delta_min',1e-12,'min_actual_decrease',1e-15,'pred_tol',1e-15, ...
    'curvature_check_stationarity',max(10*spec.reference_stationarity_tol, ...
    spec.reference_curvature_tol),'hessian_model','exact_lagrangian', ...
    'sinr_active_tol',spec.sinr_active_tol, ...
    'power_active_tol',spec.power_active_tol,'seed',seed+4000);
opts.terminal_options=struct('epsilon_KKT',spec.reference_stationarity_tol, ...
    'epsilon_H',spec.reference_curvature_tol,'max_iter',spec.reference_max_iter, ...
    'cg_tol',spec.reference_cg_tol,'cg_max_iter',spec.reference_cg_max_iter);
end

function gtan=local_tangent_gradient_norm(p,W,opts)
[~,g_complex]=compute_crlb_objective_grad(p.cfg,p.smat,W);
[c,Gcomplex]=compute_constraints_grad(p.cfg,p.channels,p.smat,W);
g=complex_matrix_to_real_vector(g_complex);
G=zeros(numel(g),numel(c));
for constraint_idx=1:numel(c)
    G(:,constraint_idx)= ...
        complex_matrix_to_real_vector(Gcomplex(:,:,constraint_idx));
end
active=false(size(c));
active(1:p.cfg.sys.U)=c(1:p.cfg.sys.U)>=-opts.sinr_active_tol;
active(p.cfg.sys.U+1:end)=c(p.cfg.sys.U+1:end)>=-opts.power_active_tol;
if any(active)
    rank_tol=1e-10;
    if isfield(opts,'rank_tol'); rank_tol=opts.rank_tol; end
    Q=orth(G(:,active),rank_tol);
    g=g-Q*(Q'*g);
end
gtan=norm(g);
end

function [slope,intercepts,slope_se,ci95,seed_ids]=fixed_effect_order(trials)
x=[]; y=[]; group=[]; seed_ids=[];
for trial_idx=1:numel(trials)
    q=trials{trial_idx};
    if ~trial_ok(q); continue; end
    mask=q.fit_mask;
    seed_ids(end+1)=q.seed; %#ok<AGROW>
    group_id=numel(seed_ids);
    x=[x;log10(q.tangent_residual_before(mask)).']; %#ok<AGROW>
    y=[y;log10(q.tangent_residual_after(mask)).']; %#ok<AGROW>
    group=[group;group_id*ones(nnz(mask),1)]; %#ok<AGROW>
end
slope=NaN; intercepts=[]; slope_se=NaN; ci95=[NaN NaN];
if numel(x)<3 || numel(seed_ids)<1; return; end
D=zeros(numel(x),numel(seed_ids));
for row_idx=1:numel(x); D(row_idx,group(row_idx))=1; end
X=[D x]; beta=X\y;
intercepts=beta(1:end-1); slope=beta(end);
degrees=max(numel(y)-rank(X),1);
residual=y-X*beta;
covariance=(residual'*residual/degrees)*pinv(X'*X);
slope_se=sqrt(max(covariance(end,end),0));
ci95=slope+[-1 1]*1.96*slope_se;
end

function [slope,intercept,slope_se]=loglog_fit(e0,e1)
x=log10(e0(:)); y=log10(e1(:));
X=[ones(size(x)) x]; beta=X\y;
intercept=beta(1); slope=beta(2); slope_se=NaN;
if numel(x)>2
    residual=y-X*beta;
    covariance=(residual'*residual/(numel(x)-2))*pinv(X'*X);
    slope_se=sqrt(max(covariance(2,2),0));
end
end

function value=fixed_order_rmse(e0,e1,order)
x=log10(e0(:)); y=log10(e1(:));
intercept=mean(y-order*x);
value=sqrt(mean((y-intercept-order*x).^2));
end

function Wa=phase_align(W,Wstar)
Wa=W;
for user_idx=1:size(W,2)
    alpha=Wstar(:,user_idx)'*W(:,user_idx);
    if abs(alpha)>0
        Wa(:,user_idx)=W(:,user_idx)*exp(-1i*angle(alpha));
    end
end
end

function d=quotient_distance(W,Wstar)
d=norm(phase_align(W,Wstar)-Wstar,'fro')/max(1,norm(Wstar,'fro'));
end

function x=first_or_nan(v)
if isempty(v); x=NaN; else; x=v(1); end
end

function yes=record_ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end

function yes=trial_ok(q)
yes=record_ok(q)&&isfield(q,'fit_mask')&&nnz(q.fit_mask)>=3;
end

function yes=is_attempted(q)
yes=isstruct(q)&&isfield(q,'status');
end
