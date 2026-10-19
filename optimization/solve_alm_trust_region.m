function [sol, hist] = solve_alm_trust_region(cfg, channels, smat, init_wcell, init_state)
%SOLVE_ALM_TRUST_REGION Safeguarded projected trust-region ALM.
%
% Core design:
%   - Original beamforming-vector variables, no SDR lifting.
%   - Real-domain variable z = [real(vec(Wmat)); imag(vec(Wmat))].
%   - SINR constraints are handled by shifted inequality ALM.
%   - Per-AP power constraints are handled by exact projection.
%   - Inner ALM subproblem is solved by safeguarded trust-region steps.
%   - Trust-region subproblem is solved by Steihaug-Toint truncated CG.
%   - Hessian-vector products are computed by finite differences of gradients.
%
% Important safeguards in this version:
%   1) Negative curvature in tCG is detected only for sufficiently negative
%      Rayleigh curvature, reducing false detections from finite-difference
%      noise.
%   2) The tCG model is built in the joint tangent space of active per-AP
%      power balls and nearly active SINR constraints. This reduces both
%      power-projection distortion and shifted-ALM active-set crossings.
%   3) Every TR candidate is compared with a projected-gradient fallback
%      candidate. If the TR model is inaccurate, the algorithm can still
%      make a stable first-order projected step.
%   4) A candidate is accepted only if it gives actual ALM decrease. This
%      protects against a poor or locally invalid Hessian model.
%   5) The best feasible physical-CRLB iterate is retained as an incumbent.

cfg = fill_tr_defaults(cfg);
solver_timer = tic;

% compute_crlb_objective_grad reads cfg.alm.objective_type in this package.
% Mirror the TR objective settings into cfg.alm for compatibility.
if ~isfield(cfg,'alm'); cfg.alm = struct(); end
cfg.alm.objective_type = cfg.tr.objective_type;
cfg.alm.objective_scale_mode = cfg.tr.objective_scale_mode;
cfg.alm.det_floor = cfg.tr.det_floor;

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M * Nt;
dim = 2 * N * U;

if nargin < 4 || isempty(init_wcell)
    switch lower(cfg.tr.init_mode)
        case 'identity'
            init_wcell = initialize_beamformer_identity(cfg);
        case 'zf'
            init_wcell = build_shared_theorem_aware_zf_initializer( ...
                cfg,channels,smat);
        case 'random'
            init_wcell = random_initial_beamformers_tr(cfg);
        otherwise
            error('Unknown cfg.tr.init_mode: %s', cfg.tr.init_mode);
    end
end

Wmat0 = vectorize_wcell_to_matrix(cfg, init_wcell);
Wmat0 = project_per_ap_power_matrix(cfg, Wmat0);
z = complex_matrix_to_real_vector(Wmat0);

% Objective scaling: use initial objective for raw ratio, no scaling for log.
[F0, ~] = compute_crlb_objective_grad(cfg, smat, Wmat0);
if strcmpi(cfg.tr.objective_scale_mode, 'initial') && ...
        ~strcmpi(cfg.tr.objective_type, 'log_crlb_ratio')
    obj_scale = max(abs(F0), 1e-12);
else
    obj_scale = 1;
end

mu = zeros(U,1);
rho = cfg.tr.rho0;
Delta = cfg.tr.delta0;
if nargin >= 5 && ~isempty(init_state)
    if isfield(init_state, 'mu') && numel(init_state.mu) == U
        mu = max(0, real(init_state.mu(:)));
    end
    if isfield(init_state, 'rho') && isfinite(init_state.rho)
        rho = min(max(real(init_state.rho), cfg.tr.rho0), cfg.tr.rho_max);
    end
    if isfield(init_state, 'delta') && isfinite(init_state.delta)
        Delta = min(max(real(init_state.delta), cfg.tr.delta_min), ...
            cfg.tr.delta_max);
    end
end

prev_outer_sinr_violation = inf;
prev_outer_crlb = inf;
outer_stall_count = 0;
global_iter = 0;
best_feasible_z = [];
best_feasible_crlb = inf;
best_feasible_iter = nan;
best_feasible_outer = nan;

hist = init_tr_history();
hist.objective_type = cfg.tr.objective_type;
hist.constraint_mode = 'SINR-ALM + AP-power-projection';
hist.hessian_mode = cfg.tr.hessian_mode;
hist.subproblem_solver = 'safeguarded tangent-space Steihaug-Toint truncated CG';

sol.success = false;
sol.message = 'max_outer reached';

[Linit, ginit, info_init] = value_grad_z(z);
pg_init = projected_grad_mapping_norm(cfg, z, ginit, N, U, cfg.tr.pg_step);
if cfg.tr.return_best_feasible && ...
        info_init.max_sinr_violation <= cfg.tr.incumbent_feas_tol && ...
        info_init.max_power_violation <= cfg.tr.incumbent_feas_tol
    best_feasible_z = z;
    best_feasible_crlb = info_init.crlb;
    best_feasible_iter = 0;
    best_feasible_outer = 0;
end

if cfg.tr.verbose
    fprintf(['[PTR-ALM] init_mode=%s, objective=%s, dim(z)=%d, ' ...
             'delta0=%.3e, deltaMax=%.3e, rho0=%.3e\n'], ...
        cfg.tr.init_mode, cfg.tr.objective_type, dim, Delta, cfg.tr.delta_max, rho);

    fprintf(['[PTR-ALM init k=0] AL=%.3e, CRLB=%.3e, sinrViol=%.3e, ' ...
             'powerViol=%.3e, pg=%.3e, grad=%.3e, minSINR=%.2f dB, maxPower=%.2f dBm\n'], ...
        Linit, info_init.crlb, info_init.max_sinr_violation, ...
        info_init.max_power_violation, pg_init, norm(ginit), ...
        pow2db(max(info_init.raw.min_sinr, realmin)), ...
        w2dbm(max(info_init.raw.max_ap_power, realmin)));
end

% Store k=0 explicitly for fair convergence plotting.
dummy_stats = make_dummy_cgstats('initial');
hist = store_tr_history(hist, 0, 0, 0, Linit, info_init, ...
    norm(ginit), pg_init, rho, Delta, NaN, NaN, NaN, 0, false, ...
    dummy_stats, 'initial', NaN, NaN);
hist.time_s(end+1) = toc(solver_timer);
hist = maybe_store_iterate(hist, cfg, z, 0);
hist.tangent_active_count(end+1) = 0;
hist.tangent_sinr_count(end+1) = 0;

for outer = 1:cfg.tr.max_outer

    if cfg.tr.verbose
        fprintf('\nPTR-ALM outer %02d/%02d, rho=%.3e, Delta=%.3e\n', ...
            outer, cfg.tr.max_outer, rho, Delta);
    end

    prev_inner_crlb = inf;
    inner_stall_count = 0;

    for inner = 1:cfg.tr.max_inner

        [Lval, grad_z, info] = value_grad_z(z);
        grad_norm = norm(grad_z);
        pg_norm = projected_grad_mapping_norm(cfg, z, grad_z, N, U, cfg.tr.pg_step);

        if isfinite(prev_inner_crlb)
            inner_crlb_relchg = abs(info.crlb - prev_inner_crlb) / max(1, abs(prev_inner_crlb));
        else
            inner_crlb_relchg = inf;
        end

        if inner_crlb_relchg <= cfg.tr.inner_crlb_rel_tol
            inner_stall_count = inner_stall_count + 1;
        else
            inner_stall_count = 0;
        end

        if cfg.tr.verbose && (inner == 1 || mod(inner, cfg.tr.print_every) == 0)
            fprintf(['  inner %03d: AL=%.3e, CRLB=%.3e, sinrViol=%.3e, ' ...
                     'powerViol=%.3e, pg=%.3e, grad=%.3e, Delta=%.3e\n'], ...
                inner, Lval, info.crlb, info.max_sinr_violation, ...
                info.max_power_violation, pg_norm, grad_norm, Delta);
        end

        if inner >= cfg.tr.min_inner && pg_norm <= cfg.tr.inner_pg_tol
            if cfg.tr.verbose
                fprintf('  inner stopped by projected-gradient tolerance at inner %d.\n', inner);
            end
            break;
        end

        if inner >= cfg.tr.min_inner && inner_stall_count >= cfg.tr.inner_patience && ...
                pg_norm <= cfg.tr.inner_plateau_pg_tol
            if cfg.tr.verbose
                fprintf('  inner stopped by CRLB plateau and small projected gradient at inner %d.\n', inner);
            end
            break;
        end

        if inner >= cfg.tr.min_inner && cfg.tr.stop_inner_on_feasible && ...
                info.max_sinr_violation <= cfg.tr.feas_tol
            if cfg.tr.verbose
                fprintf('  inner stopped after reaching SINR feasibility.\n');
            end
            break;
        end

        if Delta <= cfg.tr.delta_min
            if cfg.tr.verbose
                fprintf('  inner stopped because Delta <= delta_min.\n');
            end
            break;
        end

        % Raw HVP around the current point. The structured mode exploits the
        % CRLB/SINR quadratic structure and is both cheaper and less noisy
        % than full-dimensional finite differences.
        if strcmpi(cfg.tr.hessian_mode, 'structured')
            Hvfun_raw = @(v) structured_hess_vec_z(z, v);
        else
            grad_fun = @(zz) value_grad_only_z(zz);
            Hvfun_raw = @(v) hess_vec_fd_wrapper(grad_fun, z, v, cfg.tr.hv_eps);
        end

        % Joint tangent model for active AP power balls and nearly active
        % SINR inequalities. The SINR normals prevent a large local-Hessian
        % step from crossing a shifted-ALM kink in first order.
        Wcurr = real_vector_to_complex_matrix(z, N, U);
        Wcurr = project_per_ap_power_matrix(cfg, Wcurr);

        if cfg.tr.use_tangent_projection
            [Qactive,tangent_active_count,tangent_sinr_count] = ...
                build_tr_active_normal_basis( ...
                cfg,channels,smat,Wcurr,N,U);
            Ptan = @(v) v-Qactive*(Qactive'*v);
            grad_model = Ptan(grad_z);
            Hvfun_tcg = @(v) Ptan(Hvfun_raw(Ptan(v)));

            % If the tangent gradient is numerically tiny, tCG cannot help;
            % projected-gradient fallback will still be built below.
            if norm(grad_model) <= cfg.tr.tangent_grad_floor * max(1, grad_norm)
                grad_model = zeros(size(grad_z));
            end
        else
            tangent_active_count = 0;
            tangent_sinr_count = 0;
            grad_model = grad_z;
            Hvfun_tcg = Hvfun_raw;
        end

        % Candidate A: safeguarded TR/tCG step.
        tcg_opts.max_iter = min(cfg.tr.cg_max_iter, dim);
        tcg_opts.tol = cfg.tr.cg_tol;
        tcg_opts.nc_tol = cfg.tr.cg_nc_tol;
        tcg_opts.flat_tol = cfg.tr.cg_flat_tol;
        tcg_opts.verbose = false;

        allow_tcg = outer >= cfg.tr.tcg_start_outer && ...
            info.max_sinr_violation <= cfg.tr.tcg_start_sinr_violation;

        if allow_tcg && norm(grad_model) > 0
            [s_raw, cgstats] = trust_region_tcg(grad_model, Hvfun_tcg, Delta, tcg_opts);
            cand_tr = build_tr_candidate(cfg, z, s_raw, Lval, grad_z, Hvfun_raw, ...
                N, U, @value_grad_z);
        else
            s_raw = zeros(size(z));
            if ~allow_tcg
                cgstats = make_dummy_cgstats('tcg_delayed_until_feasible');
            else
                cgstats = make_dummy_cgstats('zero_tangent_gradient');
            end
            cand_tr = make_empty_candidate('trust_region');
        end

        % Candidate B: projected-gradient fallback. This is the robustness
        % guard that prevents PTR-ALM from performing much worse than the
        % projected first-order method when the Hessian model is inaccurate.
        if cfg.tr.use_pg_fallback
            cand_pg = build_projected_gradient_candidate(cfg, z, grad_z, Lval, ...
                info.max_sinr_violation, Delta, N, U, @value_grad_z);
        else
            cand_pg = make_empty_candidate('projected_gradient');
        end

        % Select the best actual ALM-decreasing candidate.
        selected = make_empty_candidate('none');
        selected.L = Lval;
        selected.info = info;
        selected.z = z;
        selected.step_norm = 0;

        if cand_tr.valid && cand_tr.L < selected.L
            selected = cand_tr;
        end
        if cand_pg.valid && cand_pg.L < selected.L
            selected = cand_pg;
        end

        accepted = selected.valid && selected.L < Lval - cfg.tr.min_actual_decrease;

        if accepted
            z = selected.z;
            candidate_type = selected.type;
        else
            candidate_type = 'rejected';
        end

        % Evaluate current accepted point for history consistency.
        [Lrecord, grecord, current_info] = value_grad_z(z);
        grad_record_norm = norm(grecord);
        pg_record_norm = projected_grad_mapping_norm(cfg, z, grecord, N, U, cfg.tr.pg_step);
        if cfg.tr.return_best_feasible && ...
                current_info.max_sinr_violation <= ...
                cfg.tr.incumbent_feas_tol && ...
                current_info.max_power_violation <= ...
                cfg.tr.incumbent_feas_tol && ...
                current_info.crlb < best_feasible_crlb
            best_feasible_z = z;
            best_feasible_crlb = current_info.crlb;
            best_feasible_iter = global_iter+1;
            best_feasible_outer = outer;
            hist.incumbent_iter(end+1) = global_iter+1;
            hist.incumbent_outer(end+1) = outer;
            hist.incumbent_crlb(end+1) = current_info.crlb;
        end

        % Trust-region radius update based on the TR model quality. If the
        % PG fallback wins, shrink mildly; if no candidate is accepted,
        % shrink normally.
        boundary_step = norm(s_raw) >= cfg.tr.boundary_tol * max(Delta, realmin);
        if strcmp(candidate_type, 'trust_region')
            if cand_tr.ratio < cfg.tr.eta_shrink || ~isfinite(cand_tr.ratio)
                Delta = max(cfg.tr.gamma_shrink * Delta, cfg.tr.delta_min);
            elseif cand_tr.ratio > cfg.tr.eta_expand && boundary_step
                Delta = min(cfg.tr.gamma_expand * Delta, cfg.tr.delta_max);
            end
        elseif strcmp(candidate_type, 'projected_gradient')
            % A PG fallback win does not necessarily mean the TR model is bad;
            % early feasibility-restoring PG steps can simply be much larger.
            % If the TR candidate itself has good model agreement and hit the
            % boundary, expand the radius so second-order steps can compete
            % in later iterations. Otherwise keep the radius nearly unchanged
            % instead of steadily shrinking it to delta_min.
            if cand_tr.evaluated && isfinite(cand_tr.ratio) && ...
                    cand_tr.ratio > cfg.tr.eta_expand && boundary_step
                Delta = min(cfg.tr.gamma_expand * Delta, cfg.tr.delta_max);
            elseif cand_tr.evaluated && ...
                    (~isfinite(cand_tr.ratio) || ...
                    cand_tr.ratio < cfg.tr.eta_shrink)
                Delta = max(cfg.tr.gamma_shrink * Delta, cfg.tr.delta_min);
            else
                Delta = min(max(cfg.tr.gamma_pg_success * Delta, cfg.tr.delta_min), cfg.tr.delta_max);
            end
        else
            Delta = max(cfg.tr.gamma_shrink * Delta, cfg.tr.delta_min);
        end

        global_iter = global_iter + 1;
        hist = store_tr_history(hist, global_iter, outer, inner, Lrecord, current_info, ...
            grad_record_norm, pg_record_norm, rho, Delta, cand_tr.ratio, ...
            cand_tr.pred, cand_tr.ared, selected.step_norm, accepted, cgstats, ...
            candidate_type, cand_pg.ared, cand_pg.step_norm);
        hist.time_s(end+1) = toc(solver_timer);
        hist = maybe_store_iterate(hist, cfg, z, global_iter);
        hist.tangent_active_count(end+1) = tangent_active_count;
        hist.tangent_sinr_count(end+1) = tangent_sinr_count;

        if cfg.tr.verbose && (mod(inner, cfg.tr.print_every) == 0 || ...
                cgstats.negative_curvature || strcmp(candidate_type,'projected_gradient'))
            fprintf(['    step: type=%s, acc=%d, TRratio=%.3e, TRared=%.3e, ' ...
                     'PGared=%.3e, step=%.3e, cg=%d, nc=%d, reason=%s\n'], ...
                candidate_type, accepted, cand_tr.ratio, cand_tr.ared, ...
                cand_pg.ared, selected.step_norm, cgstats.iter, ...
                cgstats.negative_curvature, cgstats.reason);
        end

        % Keep the pre-step CRLB as the baseline for the next iteration.
        % The next loop evaluates the accepted point as its current iterate;
        % comparing it with current_info.crlb would always give zero change
        % and falsely trigger the CRLB-plateau stop after a few iterations.
        prev_inner_crlb = info.crlb;

        if inner >= cfg.tr.min_inner && selected.step_norm <= cfg.tr.step_tol
            if cfg.tr.verbose
                fprintf('  inner stopped by small accepted effective step at inner %d.\n', inner);
            end
            break;
        end
    end

    % ============================================================== 
    % Outer SINR multiplier/rho update and stopping check
    % ============================================================== 

    Wmat = real_vector_to_complex_matrix(z, N, U);
    Wmat = project_per_ap_power_matrix(cfg, Wmat);
    z = complex_matrix_to_real_vector(Wmat);

    [c_all, ~, raw] = compute_constraints_grad(cfg, channels, smat, Wmat);
    c_sinr = c_all(1:U);
    outer_sinr_violation = max(max(c_sinr,0));
    if numel(c_all) >= U + M
        outer_power_violation = max(max(c_all(U+1:U+M),0));
    else
        outer_power_violation = max(0, max(raw.max_ap_power) / max(cfg.power.Pmax,realmin) - 1);
    end
    outer_violation = max(outer_sinr_violation, outer_power_violation);

    % Multipliers only for SINR constraints.
    mu = max(0, mu + rho * c_sinr);

    [Lcheck, grad_check_z, outinfo] = value_grad_z(z);
    outer_pg_norm = projected_grad_mapping_norm(cfg, z, grad_check_z, N, U, cfg.tr.pg_step);
    outer_grad_norm = norm(grad_check_z);
    current_crlb = outinfo.crlb;
    if cfg.tr.return_best_feasible && ...
            outer_sinr_violation <= cfg.tr.incumbent_feas_tol && ...
            outer_power_violation <= cfg.tr.incumbent_feas_tol && ...
            current_crlb < best_feasible_crlb
        best_feasible_z = z;
        best_feasible_crlb = current_crlb;
        best_feasible_iter = global_iter;
        best_feasible_outer = outer;
        hist.incumbent_iter(end+1) = global_iter;
        hist.incumbent_outer(end+1) = outer;
        hist.incumbent_crlb(end+1) = current_crlb;
    end

    if isfinite(prev_outer_crlb)
        outer_crlb_relchg = abs(current_crlb - prev_outer_crlb) / max(1, abs(prev_outer_crlb));
    else
        outer_crlb_relchg = inf;
    end

    lambda_min_est = NaN;
    neg_curv_violation = NaN;
    if cfg.tr.enable_curvature_diagnostic
        grad_fun = @(zz) value_grad_only_z(zz);
        Hvfun_raw = @(v) hess_vec_fd_wrapper(grad_fun, z, v, cfg.tr.hv_eps);
        lanczos_opts.max_iter = cfg.tr.lanczos_max_iter;
        lanczos_opts.tol = cfg.tr.lanczos_tol;
        [lambda_min_est, ~] = estimate_min_curvature_lanczos(Hvfun_raw, dim, lanczos_opts);
        neg_curv_violation = max(0, -lambda_min_est);
    end

    hist.outer_violation(end+1) = outer_violation;
    hist.outer_sinr_violation(end+1) = outer_sinr_violation;
    hist.outer_power_violation(end+1) = outer_power_violation;
    hist.outer_crlb(end+1) = current_crlb;
    hist.outer_rho(end+1) = rho;
    hist.outer_delta(end+1) = Delta;
    hist.outer_grad_norm(end+1) = outer_grad_norm;
    hist.outer_pg_norm(end+1) = outer_pg_norm;
    hist.outer_crlb_relchg(end+1) = outer_crlb_relchg;
    hist.outer_alm_obj(end+1) = Lcheck;
    hist.outer_lambda_min_est(end+1) = lambda_min_est;
    hist.outer_neg_curv_violation(end+1) = neg_curv_violation;

    stop_feas = outer_violation <= cfg.tr.feas_tol;
    stop_stat = outer_pg_norm <= cfg.tr.outer_pg_tol;
    stop_crlb = outer_crlb_relchg <= cfg.tr.outer_crlb_rel_tol;
    if cfg.tr.enable_curvature_diagnostic
        stop_curv = neg_curv_violation <= cfg.tr.neg_curv_tol;
    else
        stop_curv = true;
    end

    if stop_feas && stop_crlb
        outer_stall_count = outer_stall_count + 1;
    else
        outer_stall_count = 0;
    end
    hist.outer_stall_count(end+1) = outer_stall_count;

    if cfg.tr.verbose
        if cfg.tr.enable_curvature_diagnostic
            fprintf(['  outer end: sinrViol=%.3e, powerViol=%.3e, pg=%.3e, ' ...
                     'CRLB=%.3e, relCRLB=%.3e, lambdaMin=%.3e, rho=%.3e, Delta=%.3e\n'], ...
                outer_sinr_violation, outer_power_violation, outer_pg_norm, ...
                current_crlb, outer_crlb_relchg, lambda_min_est, rho, Delta);
        else
            fprintf(['  outer end: sinrViol=%.3e, powerViol=%.3e, pg=%.3e, ' ...
                     'CRLB=%.3e, relCRLB=%.3e, rho=%.3e, Delta=%.3e\n'], ...
                outer_sinr_violation, outer_power_violation, outer_pg_norm, ...
                current_crlb, outer_crlb_relchg, rho, Delta);
        end
    end

    if outer >= cfg.tr.min_outer && stop_feas && stop_stat && stop_curv
        sol.success = true;
        sol.message = 'feasibility + projected stationarity + curvature reached';
        break;
    end

    if outer >= cfg.tr.min_outer && stop_feas && cfg.tr.stop_on_feasible_only
        sol.success = true;
        sol.message = 'feasibility reached';
        break;
    end

    if outer >= cfg.tr.min_outer && stop_feas && ~cfg.tr.require_stationarity_for_stop && ...
            outer_stall_count >= cfg.tr.outer_patience
        sol.success = true;
        sol.message = 'feasibility + CRLB plateau reached';
        break;
    end

    if outer_sinr_violation > cfg.tr.feas_tol
        if outer == 1 || ~isfinite(prev_outer_sinr_violation)
            rho = min(cfg.tr.rho_growth * rho, cfg.tr.rho_max);
        elseif outer_sinr_violation > cfg.tr.rho_update_threshold * prev_outer_sinr_violation
            rho = min(cfg.tr.rho_growth * rho, cfg.tr.rho_max);
        end
    end

    prev_outer_sinr_violation = outer_sinr_violation;
    prev_outer_crlb = current_crlb;
end

Wmat = real_vector_to_complex_matrix(z, N, U);
Wmat = project_per_ap_power_matrix(cfg, Wmat);
wcell = matrix_to_wcell(cfg, Wmat);
metrics = evaluate_beamforming_solution(cfg, channels, smat, wcell);

returned_best_incumbent = false;
if cfg.tr.return_best_feasible && ~isempty(best_feasible_z)
    Wbest = real_vector_to_complex_matrix(best_feasible_z,N,U);
    Wbest = project_per_ap_power_matrix(cfg,Wbest);
    wbest = matrix_to_wcell(cfg,Wbest);
    metrics_best = evaluate_beamforming_solution( ...
        cfg,channels,smat,wbest);
    current_feasible = all(metrics.sinr >= cfg.comm.Gamma * ...
        (1-cfg.tr.incumbent_feas_tol)) && ...
        all(metrics.ap_power <= cfg.power.Pmax * ...
        (1+cfg.tr.incumbent_feas_tol));
    best_is_feasible = all(metrics_best.sinr >= cfg.comm.Gamma * ...
        (1-cfg.tr.incumbent_feas_tol)) && ...
        all(metrics_best.ap_power <= cfg.power.Pmax * ...
        (1+cfg.tr.incumbent_feas_tol));
    if best_is_feasible && ...
            (~current_feasible || metrics_best.crlb < metrics.crlb)
        Wmat = Wbest;
        wcell = wbest;
        metrics = metrics_best;
        returned_best_incumbent = true;
        sol.message = [sol.message '; returned best feasible incumbent'];
    end
end

sol.wcell = wcell;
sol.Wmat = Wmat;
sol.z = complex_matrix_to_real_vector(Wmat);
sol.initial_Wmat = Wmat0;
sol.initial_wcell = matrix_to_wcell(cfg,Wmat0);
sol.initial_z = complex_matrix_to_real_vector(Wmat0);
sol.mu = mu;
sol.rho = rho;
sol.delta = Delta;
sol.state = struct('mu', mu, 'rho', rho, 'delta', Delta);
sol.obj_scale = obj_scale;
sol.metrics = metrics;
sol.hist = hist;
sol.returned_best_incumbent = returned_best_incumbent;
sol.best_feasible_crlb = best_feasible_crlb;
sol.best_feasible_iter = best_feasible_iter;
sol.best_feasible_outer = best_feasible_outer;

    function [Lval, grad_z, info] = value_grad_z(zz)
        Wquery = real_vector_to_complex_matrix(zz, N, U);
        [Lval, grad_complex, info] = sinr_augmented_lagrangian_value_grad( ...
            cfg, channels, smat, Wquery, mu, rho, obj_scale);
        grad_z = complex_matrix_to_real_vector(grad_complex);
    end

    function grad_z = value_grad_only_z(zz)
        [~, grad_z, ~] = value_grad_z(zz);
    end

    function Hv_z = structured_hess_vec_z(zz, vv)
        Wquery = real_vector_to_complex_matrix(zz, N, U);
        Vquery = real_vector_to_complex_matrix(vv, N, U);
        Hv_complex = sinr_augmented_lagrangian_hess_vec( ...
            cfg, channels, smat, Wquery, Vquery, mu, rho, obj_scale);
        Hv_z = complex_matrix_to_real_vector(Hv_complex);
    end

end

%% ========================================================================
%  Candidate builders
%% ========================================================================

function cand = build_tr_candidate(cfg, z, s_raw, Lval, grad_z, Hvfun_raw, N, U, value_grad_handle)

cand = make_empty_candidate('trust_region');

ztrial_raw = z + s_raw;
Wtrial = real_vector_to_complex_matrix(ztrial_raw, N, U);
Wtrial = project_per_ap_power_matrix(cfg, Wtrial);
ztrial = complex_matrix_to_real_vector(Wtrial);
s_eff = ztrial - z;
step_norm = norm(s_eff) / max(1, norm(z));

cand.z = ztrial;
cand.step_norm = step_norm;

if norm(s_eff) <= cfg.tr.step_tol
    return;
end

[Ltrial, ~, info_trial] = value_grad_handle(ztrial);
% The trial was evaluated even when it increases the objective and is thus
% invalid for acceptance. Radius calibration must still use its bad model
% ratio; otherwise a rejected TR model can remain permanently oversized.
cand.evaluated = true;
Hs_eff = Hvfun_raw(s_eff);
pred = -(grad_z' * s_eff + 0.5 * (s_eff' * Hs_eff));
ared = Lval - Ltrial;

if pred > cfg.tr.pred_tol && isfinite(ared) && isfinite(pred)
    ratio = ared / pred;
else
    ratio = -inf;
end

cand.L = Ltrial;
cand.info = info_trial;
cand.pred = pred;
cand.ared = ared;
cand.ratio = ratio;

if ~isfinite(Ltrial) || ~isfinite(ared)
    return;
end

actual_decrease_ok = ared > cfg.tr.min_actual_decrease;
ratio_ok = ratio >= cfg.tr.eta_accept;

cand.valid = actual_decrease_ok && ...
    (ratio_ok || cfg.tr.accept_actual_decrease_without_ratio);

end

function cand = build_projected_gradient_candidate(cfg, z, grad_z, Lval, ...
    current_sinr_violation, Delta, N, U, value_grad_handle)

cand = make_empty_candidate('projected_gradient');

if norm(grad_z) <= 0
    return;
end

direction = -grad_z;
dir_norm = norm(direction);
if dir_norm > cfg.tr.grad_clip
    direction = direction * (cfg.tr.grad_clip / dir_norm);
    dir_norm = cfg.tr.grad_clip;
end

alpha = cfg.tr.pg_alpha0;
if cfg.tr.pg_respect_delta
    max_raw_step = cfg.tr.pg_max_step_fraction * max(Delta, cfg.tr.delta_min);
    alpha = min(alpha, max_raw_step / max(dir_norm, realmin));
end

best_L = inf;
best_z = z;
best_info = [];
best_step_norm = 0;
best_ared = -inf;
success = false;

for bt = 1:cfg.tr.pg_max_backtracking
    ztrial_raw = z + alpha * direction;
    Wtrial = real_vector_to_complex_matrix(ztrial_raw, N, U);
    Wtrial = project_per_ap_power_matrix(cfg, Wtrial);
    ztrial = complex_matrix_to_real_vector(Wtrial);
    s_eff = ztrial - z;
    step_norm = norm(s_eff) / max(1, norm(z));

    if norm(s_eff) <= cfg.tr.step_tol
        alpha = alpha * cfg.tr.pg_backtracking_beta;
        continue;
    end

    [Ltrial, ~, info_trial] = value_grad_handle(ztrial);
    if cfg.tr.clip_feasibility_crossing && ...
            current_sinr_violation > cfg.tr.feas_tol && ...
            info_trial.max_sinr_violation <= cfg.tr.feas_tol
        [ztrial, Ltrial, info_trial] = clip_to_feasibility_boundary( ...
            cfg, z, ztrial, value_grad_handle);
        s_eff = ztrial - z;
        step_norm = norm(s_eff) / max(1, norm(z));
    end
    slope_eff = grad_z' * s_eff;
    ared = Lval - Ltrial;

    if isfinite(Ltrial) && Ltrial < best_L
        best_L = Ltrial;
        best_z = ztrial;
        best_info = info_trial;
        best_step_norm = step_norm;
        best_ared = ared;
    end

    if slope_eff < 0 && isfinite(Ltrial) && ...
            Ltrial <= Lval + cfg.tr.pg_armijo_c * slope_eff && ...
            ared > cfg.tr.min_actual_decrease
        success = true;
        break;
    end

    alpha = alpha * cfg.tr.pg_backtracking_beta;
end

cand.z = best_z;
cand.L = best_L;
cand.info = best_info;
cand.step_norm = best_step_norm;
cand.ared = best_ared;
cand.pred = NaN;
cand.ratio = NaN;

if isempty(best_info) || ~isfinite(best_L)
    return;
end

cand.valid = success || ...
    (cfg.tr.pg_accept_best_decrease && best_ared > cfg.tr.min_actual_decrease);

end

function [z_feas, L_feas, info_feas] = clip_to_feasibility_boundary( ...
    cfg, z_infeas, z_feas, value_grad_handle)
% The projected AP-power set is convex, so this segment remains power feasible.
low = 0;
high = 1;
[L_feas, ~, info_feas] = value_grad_handle(z_feas);

for k = 1:cfg.tr.feasibility_crossing_bisection_iters
    mid = 0.5 * (low + high);
    z_mid = z_infeas + mid * (z_feas - z_infeas);
    [L_mid, ~, info_mid] = value_grad_handle(z_mid);
    if info_mid.max_sinr_violation <= cfg.tr.feas_tol
        high = mid;
        z_feas = z_mid;
        L_feas = L_mid;
        info_feas = info_mid;
    else
        low = mid;
    end
end
end

function cand = make_empty_candidate(type)
cand.type = type;
cand.valid = false;
cand.z = [];
cand.L = inf;
cand.info = [];
cand.pred = NaN;
cand.ared = -inf;
cand.ratio = -inf;
cand.step_norm = 0;
cand.evaluated = false;
end

%% ========================================================================
%  Hessian-vector and tangent projection helpers
%% ========================================================================

function Hv = hess_vec_fd_wrapper(grad_fun, z, v, eps_base)
[Hv, ~] = hess_vec_fd(grad_fun, z, v, eps_base);
Hv = real(Hv(:));
end

function [Qactive,active_count,sinr_count] = build_tr_active_normal_basis( ...
    cfg,channels,smat,Wref,N,U)
%BUILD_TR_ACTIVE_NORMAL_BASIS Joint tangent model for active power balls
%and nearly active SINR inequalities. Including the SINR normals prevents a
%large tCG step from crossing the shifted-ALM kink in first order.

[c_all,Gcomplex] = compute_constraints_grad(cfg,channels,smat,Wref);
active = false(numel(c_all),1);
if cfg.tr.use_sinr_tangent_projection
    active(1:U) = c_all(1:U) >= -cfg.tr.sinr_tangent_active_tol;
end
power_idx = U+(1:cfg.sys.M);
power_idx = power_idx(power_idx <= numel(c_all));
active(power_idx) = c_all(power_idx) >= -cfg.tr.power_active_tol;

idx = find(active);
Gz = zeros(2*N*U,numel(idx));
for k = 1:numel(idx)
    Gz(:,k) = complex_matrix_to_real_vector(Gcomplex(:,:,idx(k)));
end
if isempty(Gz)
    Qactive = zeros(2*N*U,0);
else
    Qactive = orth(Gz,cfg.tr.tangent_rank_tol);
end
active_count = size(Qactive,2);
sinr_count = nnz(active(1:U));
end

%% ========================================================================
%  History
%% ========================================================================

function hist = init_tr_history()
hist.iter = [];
hist.time_s = [];
hist.outer = [];
hist.inner = [];
hist.alm_obj = [];
hist.crlb = [];
hist.crlb_ratio = [];
hist.optim_obj = [];
hist.max_violation = [];
hist.max_sinr_violation = [];
hist.max_power_violation = [];
hist.min_sinr_dB = [];
hist.max_ap_power_dBm = [];
hist.grad_norm = [];
hist.pg_norm = [];
hist.rho = [];
hist.delta = [];
hist.tr_ratio = [];
hist.pred_reduction = [];
hist.actual_reduction = [];
hist.pg_actual_reduction = [];
hist.pg_step_norm = [];
hist.step_norm = [];
hist.accepted = [];
hist.candidate_type = {};
hist.cg_iter = [];
hist.cg_reason = {};
hist.negative_curvature_flag = [];
hist.near_zero_curvature_flag = [];
hist.min_curvature_est = [];
hist.outer_violation = [];
hist.outer_sinr_violation = [];
hist.outer_power_violation = [];
hist.outer_crlb = [];
hist.outer_rho = [];
hist.outer_delta = [];
hist.outer_grad_norm = [];
hist.outer_pg_norm = [];
hist.outer_crlb_relchg = [];
hist.outer_alm_obj = [];
hist.outer_lambda_min_est = [];
hist.outer_neg_curv_violation = [];
hist.outer_stall_count = [];
hist.snapshot_iter = [];
hist.z_snapshots = [];
hist.tangent_active_count = [];
hist.tangent_sinr_count = [];
hist.incumbent_iter = [];
hist.incumbent_outer = [];
hist.incumbent_crlb = [];
end

function hist = store_tr_history(hist, iter, outer, inner, Lval, info, grad_norm, pg_norm, ...
    rho, Delta, tr_ratio, pred, ared, step_norm, accepted, cgstats, candidate_type, ...
    pg_ared, pg_step_norm)

hist.iter(end+1) = iter;
hist.outer(end+1) = outer;
hist.inner(end+1) = inner;
hist.alm_obj(end+1) = Lval;
hist.crlb(end+1) = info.crlb;
hist.crlb_ratio(end+1) = info.ratio;
hist.optim_obj(end+1) = info.optim_obj;
hist.max_violation(end+1) = info.max_violation;
hist.max_sinr_violation(end+1) = info.max_sinr_violation;
hist.max_power_violation(end+1) = info.max_power_violation;
hist.min_sinr_dB(end+1) = pow2db(max(info.raw.min_sinr, realmin));
hist.max_ap_power_dBm(end+1) = w2dbm(max(info.raw.max_ap_power, realmin));
hist.grad_norm(end+1) = grad_norm;
hist.pg_norm(end+1) = pg_norm;
hist.rho(end+1) = rho;
hist.delta(end+1) = Delta;
hist.tr_ratio(end+1) = tr_ratio;
hist.pred_reduction(end+1) = pred;
hist.actual_reduction(end+1) = ared;
hist.pg_actual_reduction(end+1) = pg_ared;
hist.pg_step_norm(end+1) = pg_step_norm;
hist.step_norm(end+1) = step_norm;
hist.accepted(end+1) = accepted;
hist.candidate_type{end+1} = candidate_type;
hist.cg_iter(end+1) = cgstats.iter;
hist.cg_reason{end+1} = cgstats.reason;
hist.negative_curvature_flag(end+1) = cgstats.negative_curvature;
if isfield(cgstats,'near_zero_curvature')
    hist.near_zero_curvature_flag(end+1) = cgstats.near_zero_curvature;
else
    hist.near_zero_curvature_flag(end+1) = false;
end
hist.min_curvature_est(end+1) = cgstats.min_curvature_est;
end

function hist = maybe_store_iterate(hist, cfg, z, iter)
if ~isfield(cfg, 'diagnostics') || ...
        ~isfield(cfg.diagnostics, 'store_iterates') || ...
        ~cfg.diagnostics.store_iterates
    return;
end
if isfield(cfg.diagnostics, 'store_stride')
    stride = max(1, round(cfg.diagnostics.store_stride));
else
    stride = 1;
end
if isempty(hist.snapshot_iter) || mod(iter, stride) == 0
    hist.snapshot_iter(end+1) = iter;
    hist.z_snapshots(:,end+1) = z(:);
end
end

function stats = make_dummy_cgstats(reason)
stats.iter = 0;
stats.reason = reason;
stats.negative_curvature = false;
stats.near_zero_curvature = false;
stats.hit_boundary = false;
stats.min_curvature_est = NaN;
stats.final_res_norm = NaN;
stats.model_value = NaN;
stats.model_decrease = NaN;
stats.step_norm = 0;
end

%% ========================================================================
%  Defaults and initialization
%% ========================================================================

function cfg = fill_tr_defaults(cfg)
if ~isfield(cfg,'tr'); cfg.tr = struct(); end
if ~isfield(cfg.tr,'max_outer'); cfg.tr.max_outer = 60; end
if ~isfield(cfg.tr,'max_inner'); cfg.tr.max_inner = 60; end
if ~isfield(cfg.tr,'min_outer'); cfg.tr.min_outer = 3; end
if ~isfield(cfg.tr,'min_inner'); cfg.tr.min_inner = 3; end
if ~isfield(cfg.tr,'rho0'); cfg.tr.rho0 = 1; end
if ~isfield(cfg.tr,'rho_growth'); cfg.tr.rho_growth = 3; end
if ~isfield(cfg.tr,'rho_max'); cfg.tr.rho_max = 1e10; end
if ~isfield(cfg.tr,'rho_update_threshold'); cfg.tr.rho_update_threshold = 0.9; end
if ~isfield(cfg.tr,'feas_tol'); cfg.tr.feas_tol = 1e-4; end

% Conservative trust-region radii. The first package used delta_max=1,
% which is too large when most APs are active at Pmax and the projection
% strongly distorts the tCG step.
if ~isfield(cfg.tr,'delta0'); cfg.tr.delta0 = 1e-3; end
if ~isfield(cfg.tr,'delta_max'); cfg.tr.delta_max = 5e-2; end
if ~isfield(cfg.tr,'delta_min'); cfg.tr.delta_min = 1e-8; end
if ~isfield(cfg.tr,'eta_accept'); cfg.tr.eta_accept = 0.10; end
if ~isfield(cfg.tr,'eta_shrink'); cfg.tr.eta_shrink = 0.25; end
if ~isfield(cfg.tr,'eta_expand'); cfg.tr.eta_expand = 0.75; end
if ~isfield(cfg.tr,'gamma_shrink'); cfg.tr.gamma_shrink = 0.5; end
if ~isfield(cfg.tr,'gamma_expand'); cfg.tr.gamma_expand = 1.5; end
if ~isfield(cfg.tr,'gamma_pg_success'); cfg.tr.gamma_pg_success = 1.0; end
if ~isfield(cfg.tr,'boundary_tol'); cfg.tr.boundary_tol = 0.9; end
if ~isfield(cfg.tr,'pred_tol'); cfg.tr.pred_tol = 1e-14; end
if ~isfield(cfg.tr,'min_actual_decrease'); cfg.tr.min_actual_decrease = 1e-12; end
if ~isfield(cfg.tr,'accept_actual_decrease_without_ratio'); cfg.tr.accept_actual_decrease_without_ratio = true; end

if ~isfield(cfg.tr,'cg_max_iter'); cfg.tr.cg_max_iter = 60; end
if ~isfield(cfg.tr,'cg_tol'); cfg.tr.cg_tol = 1e-3; end
if ~isfield(cfg.tr,'cg_nc_tol'); cfg.tr.cg_nc_tol = 1e-8; end
if ~isfield(cfg.tr,'cg_flat_tol'); cfg.tr.cg_flat_tol = 1e-12; end
if ~isfield(cfg.tr,'hessian_mode'); cfg.tr.hessian_mode = 'structured'; end
if ~isfield(cfg.tr,'hv_eps'); cfg.tr.hv_eps = 1e-5; end

if ~isfield(cfg.tr,'use_tangent_projection'); cfg.tr.use_tangent_projection = true; end
if ~isfield(cfg.tr,'power_active_tol'); cfg.tr.power_active_tol = 1e-5; end
if ~isfield(cfg.tr,'tangent_grad_floor'); cfg.tr.tangent_grad_floor = 1e-12; end
if ~isfield(cfg.tr,'use_sinr_tangent_projection')
    cfg.tr.use_sinr_tangent_projection = true;
end
if ~isfield(cfg.tr,'sinr_tangent_active_tol')
    cfg.tr.sinr_tangent_active_tol = 1e-3;
end
if ~isfield(cfg.tr,'tangent_rank_tol'); cfg.tr.tangent_rank_tol = 1e-10; end
if ~isfield(cfg.tr,'return_best_feasible'); cfg.tr.return_best_feasible = true; end
if ~isfield(cfg.tr,'incumbent_feas_tol')
    cfg.tr.incumbent_feas_tol = 1e-4;
end

if ~isfield(cfg.tr,'use_pg_fallback'); cfg.tr.use_pg_fallback = true; end
if ~isfield(cfg.tr,'pg_step'); cfg.tr.pg_step = 1e-3; end
if ~isfield(cfg.tr,'pg_alpha0'); cfg.tr.pg_alpha0 = 1e-3; end
if ~isfield(cfg.tr,'pg_respect_delta'); cfg.tr.pg_respect_delta = false; end
if ~isfield(cfg.tr,'pg_max_step_fraction'); cfg.tr.pg_max_step_fraction = 1.0; end
if ~isfield(cfg.tr,'pg_max_backtracking'); cfg.tr.pg_max_backtracking = 30; end
if ~isfield(cfg.tr,'pg_backtracking_beta'); cfg.tr.pg_backtracking_beta = 0.5; end
if ~isfield(cfg.tr,'pg_armijo_c'); cfg.tr.pg_armijo_c = 1e-4; end
if ~isfield(cfg.tr,'pg_accept_best_decrease'); cfg.tr.pg_accept_best_decrease = true; end
if ~isfield(cfg.tr,'grad_clip'); cfg.tr.grad_clip = 1e3; end
if ~isfield(cfg.tr,'tcg_start_outer'); cfg.tr.tcg_start_outer = 3; end
if ~isfield(cfg.tr,'tcg_start_sinr_violation'); cfg.tr.tcg_start_sinr_violation = cfg.tr.feas_tol; end

if ~isfield(cfg.tr,'inner_pg_tol'); cfg.tr.inner_pg_tol = 1e-4; end
if ~isfield(cfg.tr,'inner_plateau_pg_tol'); cfg.tr.inner_plateau_pg_tol = 1e-2; end
if ~isfield(cfg.tr,'outer_pg_tol'); cfg.tr.outer_pg_tol = 1e-4; end
if ~isfield(cfg.tr,'inner_crlb_rel_tol'); cfg.tr.inner_crlb_rel_tol = 1e-6; end
if ~isfield(cfg.tr,'outer_crlb_rel_tol'); cfg.tr.outer_crlb_rel_tol = 1e-5; end
if ~isfield(cfg.tr,'inner_patience'); cfg.tr.inner_patience = 6; end
if ~isfield(cfg.tr,'outer_patience'); cfg.tr.outer_patience = 2; end
if ~isfield(cfg.tr,'step_tol'); cfg.tr.step_tol = 1e-10; end
if ~isfield(cfg.tr,'neg_curv_tol'); cfg.tr.neg_curv_tol = 1e-5; end
if ~isfield(cfg.tr,'enable_curvature_diagnostic'); cfg.tr.enable_curvature_diagnostic = false; end
if ~isfield(cfg.tr,'lanczos_max_iter'); cfg.tr.lanczos_max_iter = 20; end
if ~isfield(cfg.tr,'lanczos_tol'); cfg.tr.lanczos_tol = 1e-8; end
if ~isfield(cfg.tr,'require_stationarity_for_stop'); cfg.tr.require_stationarity_for_stop = false; end
if ~isfield(cfg.tr,'stop_inner_on_feasible'); cfg.tr.stop_inner_on_feasible = false; end
if ~isfield(cfg.tr,'stop_on_feasible_only'); cfg.tr.stop_on_feasible_only = false; end
if ~isfield(cfg.tr,'clip_feasibility_crossing'); cfg.tr.clip_feasibility_crossing = false; end
if ~isfield(cfg.tr,'feasibility_crossing_bisection_iters'); cfg.tr.feasibility_crossing_bisection_iters = 12; end
if ~isfield(cfg.tr,'objective_type'); cfg.tr.objective_type = 'log_crlb_ratio'; end
if ~isfield(cfg.tr,'objective_scale_mode'); cfg.tr.objective_scale_mode = 'none'; end
if ~isfield(cfg.tr,'det_floor'); cfg.tr.det_floor = 1e-30; end
if ~isfield(cfg.tr,'init_mode'); cfg.tr.init_mode = 'zf'; end
if ~isfield(cfg.tr,'identity_power_scale'); cfg.tr.identity_power_scale = cfg.opt.init_power_scale; end
if ~isfield(cfg.tr,'verbose'); cfg.tr.verbose = true; end
if ~isfield(cfg.tr,'print_every'); cfg.tr.print_every = 5; end
end

function wcell = random_initial_beamformers_tr(cfg)
M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
Pmax = cfg.power.Pmax;
wcell = cell(M,U);
for m = 1:M
    X = (randn(Nt,U) + 1j*randn(Nt,U)) / sqrt(2);
    X = X / max(norm(X,'fro'), realmin) * sqrt(0.5*Pmax);
    for u = 1:U
        wcell{m,u} = X(:,u);
    end
end
end
