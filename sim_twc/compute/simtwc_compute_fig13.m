function archive = simtwc_compute_fig13()
%SIMTWC_COMPUTE_FIG13 Conditional local quadratic-convergence diagnostic.
%
% A tightly solved COSTAR endpoint defines a phase-orbit reference. Feasible
% iterates from the same trajectory are restarted for exactly one AS-RTR
% step. The measured map e_k -> e_{k+1} tests the local quadratic regime;
% it is a numerical diagnostic, not a proof of the theorem.

s=simtwc_setup(); spec=s.fig13;
data=simtwc_initialize_data(13,spec,{'costar'}); path=simtwc_archive_path(13);
p=simtwc_make_problem(spec,spec.seed);
ropts=reference_options(spec);

if ~isfield(data,'reference') || ~ok(data.reference)
    fprintf('[Fig13] Computing high-accuracy COSTAR reference, seed=%d\n',spec.seed);
    data.reference=simtwc_try_run('costar',p,true,struct('rtr_opts',ropts));
    simtwc_save_atomic(path,data);
end
if ~ok(data.reference)
    error('Fig. 13 reference solve failed: %s',data.reference.message);
end

Wstar=data.reference.Wmat;
[trajectory_error,trajectory_keep]=trajectory_distances( ...
    p,data.reference.history.z_snapshots,Wstar);
candidate_idx=select_candidates(trajectory_error,trajectory_keep,spec);
if ~isfield(data,'one_step') || numel(data.one_step)~=numel(candidate_idx)
    data.one_step=cell(1,numel(candidate_idx));
end

local_opts=ropts;
local_opts.max_iter=1;
local_opts.delta0=spec.one_step_delta;
local_opts.delta_max=spec.one_step_delta;
local_opts.curvature_policy='none';
local_opts.stationarity_tol=0;
local_opts.grad_floor=0;
local_opts.min_actual_decrease=0;
local_opts.pred_tol=0;
for j=1:numel(candidate_idx)
    if ok(data.one_step{j}); continue; end
    idx=candidate_idx(j);
    z0=data.reference.history.z_snapshots(:,idx);
    W0=real_vector_to_complex_matrix(z0,p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
    fprintf('[Fig13] One-step restart %d/%d, e0=%.3e\n', ...
        j,numel(candidate_idx),trajectory_error(idx));
    try
        [sol,hist]=solve_active_rtr_phase2( ...
            p.cfg,p.channels,p.smat,matrix_to_wcell(p.cfg,W0),local_opts);
        W1=vectorize_wcell_to_matrix(p.cfg,sol.wcell);
        kkt1=simtwc_unified_kkt(p,W1);
        accepted=~isempty(hist.accepted) && logical(hist.accepted(1));
        interior=~isempty(hist.step_norm) && ~isempty(hist.delta) && ...
            hist.step_norm(1)<0.99*hist.delta(1);
        stable=size(hist.active_membership,2)>=2 && ...
            isequal(hist.active_membership(:,1),hist.active_membership(:,end));
        data.one_step{j}=struct('status','ok','source_index',idx, ...
            'e0',trajectory_error(idx),'e1',quotient_distance(W1,Wstar), ...
            'accepted',accepted,'feasible',kkt1.r_c<=s.feasibility_tol, ...
            'active_set_stable',stable,'interior_step',interior, ...
            'step_norm',first_or_nan(hist.step_norm), ...
            'trust_radius',first_or_nan(hist.delta), ...
            'kkt_after',kkt1.r_kkt,'history',hist);
    catch ME
        data.one_step{j}=struct('status','failed','source_index',idx, ...
            'message',ME.message);
    end
    simtwc_save_atomic(path,data);
end

data.seed=spec.seed; data.reference_Wmat=Wstar;
data.trajectory_error=trajectory_error;
data.trajectory_keep=trajectory_keep;
data.candidate_indices=candidate_idx;
data.local_window_satisfied=all( ...
    trajectory_error(candidate_idx)<=spec.local_error_max);
[data.e0,data.e1,data.accepted,data.feasible,data.active_set_stable, ...
    data.interior_step]=extract(data.one_step);
fit_mask=data.accepted & data.feasible & data.active_set_stable & ...
    data.interior_step & isfinite(data.e0) & isfinite(data.e1) & ...
    data.e0>0 & data.e1>0 & data.e1<data.e0;
data.fit_mask=fit_mask;
data.loglog_slope=NaN; data.loglog_intercept=NaN;
if nnz(fit_mask)>=2
    coeff=polyfit(log10(data.e0(fit_mask)),log10(data.e1(fit_mask)),1);
    data.loglog_slope=coeff(1); data.loglog_intercept=coeff(2);
end
data.complete=all(cellfun(@ok,data.one_step));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function opts=reference_options(spec)
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
    'power_active_tol',spec.power_active_tol,'seed',spec.seed+4000);
opts.terminal_options=struct('epsilon_KKT',spec.reference_stationarity_tol, ...
    'epsilon_H',spec.reference_curvature_tol,'max_iter',spec.reference_max_iter, ...
    'cg_tol',spec.reference_cg_tol,'cg_max_iter',spec.reference_cg_max_iter);
end

function [err,keep]=trajectory_distances(p,Z,Wstar)
n=size(Z,2); err=nan(1,n); keep=false(1,n);
for k=1:n
    W=real_vector_to_complex_matrix(Z(:,k),p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
    err(k)=quotient_distance(W,Wstar);
    q=simtwc_unified_kkt(p,W);
    keep(k)=q.r_c<=simtwc_config().feasibility_tol;
    if k>1 && norm(Z(:,k)-Z(:,k-1))<=1e-13*max(1,norm(Z(:,k)))
        keep(k)=false;
    end
end
end

function idx=select_candidates(err,keep,spec)
valid=find(keep & isfinite(err) & err>spec.local_error_floor & ...
    err<=spec.local_error_max);
if isempty(valid)
    error('No feasible nonterminal trajectory point is available for Fig. 13.');
end
[~,order]=sort(err(valid)); valid=valid(order);
n=min(spec.local_points,numel(valid));
pick=unique(round(linspace(1,numel(valid),n)),'stable');
idx=valid(pick);
end

function d=quotient_distance(W,Wstar)
Wa=W;
for u=1:size(W,2)
    alpha=Wstar(:,u)'*W(:,u);
    if abs(alpha)>0; Wa(:,u)=W(:,u)*exp(-1i*angle(alpha)); end
end
d=norm(Wa-Wstar,'fro')/max(1,norm(Wstar,'fro'));
end

function [e0,e1,a,f,s,i]=extract(runs)
n=numel(runs); e0=nan(1,n); e1=nan(1,n);
a=false(1,n); f=false(1,n); s=false(1,n); i=false(1,n);
for k=1:n
    q=runs{k}; if ~ok(q); continue; end
    e0(k)=q.e0; e1(k)=q.e1; a(k)=q.accepted; f(k)=q.feasible;
    s(k)=q.active_set_stable; i(k)=q.interior_step;
end
end

function x=first_or_nan(v)
if isempty(v); x=NaN; else; x=v(1); end
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
