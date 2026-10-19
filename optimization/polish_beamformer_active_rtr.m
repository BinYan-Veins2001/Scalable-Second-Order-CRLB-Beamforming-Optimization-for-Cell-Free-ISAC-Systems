function [sol, hist] = polish_beamformer_active_rtr( ...
    cfg, channels, smat, init_wcell, opts)
%POLISH_BEAMFORMER_ACTIVE_RTR Active-set Riemannian trust-region polish.
%
% Active AP power balls are treated as a product of complex spheres. Nearly
% active SINR inequalities are added to the critical tangent space. Each
% subproblem is solved by Steihaug-Toint tCG using finite-difference Hessian-
% vector products of the original constrained Lagrangian.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(opts);

N = cfg.sys.M * cfg.sys.Nt;
U = cfg.sys.U;
dim = 2*N*U;
Wmat = vectorize_wcell_to_matrix(cfg,init_wcell);
Wmat = project_per_ap_power_matrix(cfg,Wmat);
z = complex_matrix_to_real_vector(Wmat);
Delta = opts.delta0;

hist = init_history();
sol.success = false;
sol.message = 'maximum RTR iterations reached';
lambda_min = NaN;
negative_curvature = inf;

timer = tic;
for iter = 0:opts.max_iter
    [f,grad_f,c,Gz,crlb] = derivatives(z);
    active = select_active_constraints(c);
    lambda = estimate_multipliers(Gz,grad_f,active);
    lag_grad = grad_f + Gz*lambda;
    Q = active_normal_basis(Gz,active);
    project_tangent = @(v) v - Q*(Q'*v);
    rgrad = project_tangent(grad_f);
    stationarity = norm(lag_grad) / max(1,norm(grad_f));
    primal = max([0;c]);

    Hv_raw = @(v) lagrangian_hess_vec(z,v,lambda);
    Hv_tangent = @(v) project_tangent(Hv_raw(project_tangent(v)));

    do_curvature_check = iter == 0 || mod(iter,opts.curvature_check_every) == 0 || ...
        stationarity <= opts.curvature_check_stationarity || iter == opts.max_iter;
    curvature_info = struct('ritz_vector',[],'iter',0,'ritz_values',[],'T',[]);
    if do_curvature_check
        lanczos_opts.max_iter = min(opts.lanczos_max_iter,dim);
        lanczos_opts.tol = opts.lanczos_tol;
        lanczos_opts.seed = opts.seed + iter;
        [lambda_min,curvature_info] = estimate_min_curvature_lanczos( ...
            Hv_tangent,dim,lanczos_opts);
        negative_curvature = max(0,-lambda_min);
    end

    hist = store_history(hist,iter,f,crlb,toc(timer),c,stationarity,primal, ...
        lambda_min,negative_curvature,Delta,z,active,false,NaN,0,'diagnostic');

    if primal <= opts.feas_tol && stationarity <= opts.stationarity_tol && ...
            negative_curvature <= opts.curvature_tol
        sol.success = true;
        sol.message = 'KKT and weak second-order stationarity reached';
        break;
    end
    if iter == opts.max_iter
        break;
    end

    tcg_opts.max_iter = min(opts.cg_max_iter,dim);
    tcg_opts.tol = opts.cg_tol;
    tcg_opts.nc_tol = opts.cg_nc_tol;
    tcg_opts.flat_tol = opts.cg_flat_tol;
    tcg_opts.verbose = false;

    if norm(rgrad) > opts.grad_floor
        [step,cgstats] = trust_region_tcg(rgrad,Hv_tangent,Delta,tcg_opts);
        step_type = 'tcg';
    elseif negative_curvature > opts.curvature_tol && ...
            ~isempty(curvature_info.ritz_vector)
        direction = project_tangent(curvature_info.ritz_vector);
        direction = direction/max(norm(direction),realmin);
        step = Delta*direction;
        if grad_f'*step > 0
            step = -step;
        end
        cgstats = empty_cgstats('lanczos_negative_curvature');
        cgstats.negative_curvature = true;
        step_type = 'lanczos_nc';
    else
        sol.success = primal <= opts.feas_tol && stationarity <= opts.stationarity_tol;
        sol.message = 'no first- or second-order RTR direction found';
        break;
    end

    step = project_tangent(step);
    Hstep = Hv_tangent(step);
    predicted = -(grad_f'*step + 0.5*step'*Hstep);
    ztrial = retract_active_power(z,step,c);
    ztrial = correct_sinr_retraction(ztrial,c,Delta);
    [ftrial,~,ctrial,~] = derivatives(ztrial);
    actual = f-ftrial;
    feasible_trial = max([0;ctrial]) <= opts.feas_tol;
    if predicted > opts.pred_tol && isfinite(actual)
        ratio = actual/predicted;
    else
        ratio = -inf;
    end
    accepted = feasible_trial && actual > opts.min_actual_decrease && ...
        ratio >= opts.eta_accept;

    boundary_step = norm(step) >= opts.boundary_tol*Delta;
    if ~feasible_trial || ratio < opts.eta_shrink
        Delta = max(opts.gamma_shrink*Delta,opts.delta_min);
    elseif ratio > opts.eta_expand && boundary_step
        Delta = min(opts.gamma_expand*Delta,opts.delta_max);
    end
    if accepted
        z = ztrial;
    end

    hist.accepted(end) = accepted;
    hist.ratio(end) = ratio;
    hist.step_norm(end) = norm(step);
    hist.step_type{end} = step_type;
    hist.cg_iter(end) = cgstats.iter;
    hist.cg_negative_curvature(end) = cgstats.negative_curvature;

    if Delta <= opts.delta_min && ~accepted
        sol.message = 'trust-region radius reached minimum';
        break;
    end
end
runtime = toc(timer);

Wmat = real_vector_to_complex_matrix(z,N,U);
wcell = matrix_to_wcell(cfg,Wmat);
metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
sol.wcell = wcell;
sol.Wmat = Wmat;
sol.z = z;
sol.metrics = metrics;
sol.runtime = runtime;
sol.lambda = lambda;
sol.lambda_min = lambda_min;
sol.negative_curvature = negative_curvature;
sol.stationarity = stationarity;
sol.primal_residual = primal;
sol.delta = Delta;
sol.hist = hist;

    function [f,grad_f,c,Gz,crlbq] = derivatives(zq)
        Wq = real_vector_to_complex_matrix(zq,N,U);
        [f,g_complex,gout] = compute_crlb_objective_grad(cfg,smat,Wq);
        crlbq = gout.crlb;
        [c,Gcomplex] = compute_constraints_grad(cfg,channels,smat,Wq);
        grad_f = complex_matrix_to_real_vector(g_complex);
        Gz = zeros(dim,numel(c));
        for jj = 1:numel(c)
            Gz(:,jj) = complex_matrix_to_real_vector(Gcomplex(:,:,jj));
        end
    end

    function active = select_active_constraints(cq)
        active = false(size(cq));
        active(1:U) = cq(1:U) >= -opts.sinr_active_tol;
        active(U+1:end) = cq(U+1:end) >= -opts.power_active_tol;
    end

    function Q = active_normal_basis(Gq,active)
        if any(active)
            Q = orth(Gq(:,active),opts.rank_tol);
        else
            Q = zeros(dim,0);
        end
    end

    function lambdaq = estimate_multipliers(Gq,gq,active)
        lambdaq = zeros(size(Gq,2),1);
        if any(active)
            lambdaq(active) = lsqnonneg(Gq(:,active),-gq);
        end
    end

    function Hv = lagrangian_hess_vec(zq,v,lambdaq)
        vnorm = norm(v);
        if vnorm <= realmin
            Hv = zeros(size(v));
            return;
        end
        h = opts.fd_rel_step*max(1,norm(zq))/vnorm;
        [~,gp,~,Gp] = derivatives(zq+h*v);
        [~,gm,~,Gm] = derivatives(zq-h*v);
        Hv = ((gp+Gp*lambdaq)-(gm+Gm*lambdaq))/(2*h);
        Hv = real(Hv(:));
    end

    function znew = retract_active_power(zq,step,cq)
        Wq = real_vector_to_complex_matrix(zq+step,N,U);
        Nt = cfg.sys.Nt;
        for mm = 1:cfg.sys.M
            rows = (mm-1)*Nt+(1:Nt);
            block = Wq(rows,:);
            power = real(sum(abs(block(:)).^2));
            is_active_power = cq(U+mm) >= -opts.power_active_tol;
            if is_active_power && power > 0
                block = block*sqrt(cfg.power.Pmax/power);
            elseif power > cfg.power.Pmax
                block = block*sqrt(cfg.power.Pmax/power);
            end
            Wq(rows,:) = block;
        end
        znew = complex_matrix_to_real_vector(Wq);
    end

    function zq = correct_sinr_retraction(zq,cbase,trust_radius)
        % Newton-like normal correction for the intersection of the active
        % power spheres and near-active SINR surfaces. This is part of the
        % retraction, not a separate optimization or line-search step.
        for correction_iter = 1:opts.sinr_retraction_max_iter
            [~,~,cq,Gq] = derivatives(zq);
            idx = find(cq(1:U) > -opts.sinr_retraction_margin);
            if isempty(idx)
                return;
            end

            active_power = cbase(U+1:end) >= -opts.power_active_tol;
            if any(active_power)
                Qpower = orth(Gq(:,U+find(active_power)),opts.rank_tol);
                project_power_tangent = @(v) v-Qpower*(Qpower'*v);
            else
                project_power_tangent = @(v) v;
            end
            A = project_power_tangent(Gq(:,idx));
            residual = cq(idx)+opts.sinr_retraction_margin;
            normal_matrix = A'*A + opts.sinr_retraction_regularization*eye(numel(idx));
            correction = -A*(normal_matrix\residual);
            correction = project_power_tangent(correction);
            correction_cap = opts.sinr_retraction_max_fraction*trust_radius;
            if norm(correction) > correction_cap
                correction = correction*(correction_cap/norm(correction));
            end
            zq = retract_active_power(zq,correction,cbase);
        end
    end
end

function hist = init_history()
hist.iter=[]; hist.objective=[]; hist.primal_residual=[];
hist.crlb=[]; hist.time_s=[];
hist.stationarity=[]; hist.lambda_min=[]; hist.negative_curvature=[];
hist.delta=[]; hist.snapshot_iter=[]; hist.z_snapshots=[];
hist.active_count=[]; hist.accepted=[]; hist.ratio=[]; hist.step_norm=[];
hist.step_type={}; hist.cg_iter=[]; hist.cg_negative_curvature=[];
hist.max_constraint=[];
end

function hist = store_history(hist,iter,f,crlb,time_s,c,stationarity,primal,lambda_min, ...
    negative_curvature,Delta,z,active,accepted,ratio,step_norm,step_type)
hist.iter(end+1)=iter; hist.objective(end+1)=f;
hist.crlb(end+1)=crlb; hist.time_s(end+1)=time_s;
hist.primal_residual(end+1)=primal; hist.stationarity(end+1)=stationarity;
hist.lambda_min(end+1)=lambda_min;
hist.negative_curvature(end+1)=negative_curvature;
hist.delta(end+1)=Delta; hist.snapshot_iter(end+1)=iter;
hist.z_snapshots(:,end+1)=z; hist.active_count(end+1)=nnz(active);
hist.accepted(end+1)=accepted; hist.ratio(end+1)=ratio;
hist.step_norm(end+1)=step_norm; hist.step_type{end+1}=step_type;
hist.cg_iter(end+1)=0; hist.cg_negative_curvature(end+1)=false;
hist.max_constraint(end+1)=max(c);
end

function stats = empty_cgstats(reason)
stats.iter=0; stats.reason=reason; stats.negative_curvature=false;
stats.near_zero_curvature=false; stats.hit_boundary=true;
stats.min_curvature_est=NaN; stats.final_res_norm=NaN;
stats.model_value=NaN; stats.model_decrease=NaN; stats.step_norm=0;
end

function opts = fill_defaults(opts)
if ~isfield(opts,'max_iter'); opts.max_iter=120; end
if ~isfield(opts,'delta0'); opts.delta0=5e-2; end
if ~isfield(opts,'delta_max'); opts.delta_max=5e-1; end
if ~isfield(opts,'delta_min'); opts.delta_min=1e-9; end
if ~isfield(opts,'eta_accept'); opts.eta_accept=0.05; end
if ~isfield(opts,'eta_shrink'); opts.eta_shrink=0.25; end
if ~isfield(opts,'eta_expand'); opts.eta_expand=0.75; end
if ~isfield(opts,'gamma_shrink'); opts.gamma_shrink=0.5; end
if ~isfield(opts,'gamma_expand'); opts.gamma_expand=1.8; end
if ~isfield(opts,'boundary_tol'); opts.boundary_tol=0.9; end
if ~isfield(opts,'pred_tol'); opts.pred_tol=1e-12; end
if ~isfield(opts,'min_actual_decrease'); opts.min_actual_decrease=1e-12; end
if ~isfield(opts,'sinr_active_tol'); opts.sinr_active_tol=1e-3; end
if ~isfield(opts,'power_active_tol'); opts.power_active_tol=1e-5; end
if ~isfield(opts,'sinr_retraction_margin'); opts.sinr_retraction_margin=1e-5; end
if ~isfield(opts,'sinr_retraction_max_iter'); opts.sinr_retraction_max_iter=4; end
if ~isfield(opts,'sinr_retraction_regularization'); opts.sinr_retraction_regularization=1e-12; end
if ~isfield(opts,'sinr_retraction_max_fraction'); opts.sinr_retraction_max_fraction=1.0; end
if ~isfield(opts,'rank_tol'); opts.rank_tol=1e-10; end
if ~isfield(opts,'fd_rel_step'); opts.fd_rel_step=1e-5; end
if ~isfield(opts,'cg_max_iter'); opts.cg_max_iter=30; end
if ~isfield(opts,'cg_tol'); opts.cg_tol=1e-3; end
if ~isfield(opts,'cg_nc_tol'); opts.cg_nc_tol=1e-8; end
if ~isfield(opts,'cg_flat_tol'); opts.cg_flat_tol=1e-12; end
if ~isfield(opts,'grad_floor'); opts.grad_floor=1e-10; end
if ~isfield(opts,'lanczos_max_iter'); opts.lanczos_max_iter=20; end
if ~isfield(opts,'lanczos_tol'); opts.lanczos_tol=1e-8; end
if ~isfield(opts,'curvature_check_every'); opts.curvature_check_every=10; end
if ~isfield(opts,'curvature_check_stationarity'); opts.curvature_check_stationarity=2e-3; end
if ~isfield(opts,'feas_tol'); opts.feas_tol=1e-6; end
if ~isfield(opts,'stationarity_tol'); opts.stationarity_tol=1e-3; end
if ~isfield(opts,'curvature_tol'); opts.curvature_tol=1e-5; end
if ~isfield(opts,'seed'); opts.seed=2026; end
end
