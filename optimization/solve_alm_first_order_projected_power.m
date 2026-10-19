function [sol, hist] = solve_alm_first_order_projected_power(cfg, channels, smat, init_wcell)
%SOLVE_ALM_FIRST_ORDER_PROJECTED_POWER First-order ALM with AP power projection.
%
% This companion baseline is included for comparison with PTR-ALM:
%   - SINR constraints: shifted inequality ALM;
%   - AP power constraints: exact projection.

cfg = fill_proj_alm_defaults(cfg);
solver_timer = tic;

ppfo_enabled = cfg.ppfo.enabled;
if ppfo_enabled
    perturb_stream = RandStream('mt19937ar','Seed',cfg.ppfo.seed);
else
    perturb_stream = [];
end

M = cfg.sys.M; %#ok<NASGU>
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = cfg.sys.M * Nt;

caller_supplied_initializer = nargin >= 4 && ~isempty(init_wcell);
if ~caller_supplied_initializer
    switch lower(cfg.alm.init_mode)
        case 'identity'
            init_wcell = initialize_beamformer_identity(cfg);
        case 'zf'
            init_wcell = build_shared_theorem_aware_zf_initializer( ...
                cfg,channels,smat);
        otherwise
            init_wcell = random_initial_beamformers_projected(cfg);
    end
end

Wmat0 = vectorize_wcell_to_matrix(cfg, init_wcell);
Wmat0 = project_per_ap_power_matrix(cfg, Wmat0);
initial_Wmat = Wmat0;
z = complex_matrix_to_real_vector(Wmat0);

[F0, ~] = compute_crlb_objective_grad(cfg, smat, Wmat0);
if strcmpi(cfg.alm.objective_scale_mode,'initial') && ~strcmpi(cfg.alm.objective_type,'log_crlb_ratio')
    obj_scale = max(abs(F0), 1e-12);
else
    obj_scale = 1;
end

mu = zeros(U,1);
rho = cfg.alm.rho0;
prev_outer_sinr_violation = inf;
prev_outer_crlb = inf;
outer_stall_count = 0;
global_iter = 0;
perturbation_count = 0;
last_perturb_iter = -inf;
best_feasible_z = [];
best_feasible_crlb = inf;
best_feasible_iter = nan;
best_feasible_outer = nan;

hist = init_proj_alm_history();
hist.objective_type = cfg.alm.objective_type;
hist.constraint_mode = 'SINR-ALM + AP-power-projection';
hist.ppfo_enabled = ppfo_enabled;
hist.ppfo_config = cfg.ppfo;

sol.success = false;
sol.message = 'max_outer reached';

if cfg.alm.verbose
    if ppfo_enabled
        solver_name = 'PPFO-ALM';
    else
        solver_name = 'Proj-FO-ALM';
    end
    fprintf('[%s] objective=%s, dim(z)=%d\n', ...
        solver_name,cfg.alm.objective_type,numel(z));
end

for outer = 1:cfg.alm.max_outer
    if cfg.alm.verbose
        fprintf('\nProj-FO-ALM outer %02d/%02d, rho=%.3e\n', outer, cfg.alm.max_outer, rho);
    end

    prev_inner_L = inf;
    prev_inner_crlb = inf;
    inner_stall_count = 0;
    escape_active = false;
    escape_reference_z = [];
    escape_reference_L = nan;
    escape_reference_crlb = nan;
    escape_reference_sinr_violation = nan;

    for inner = 1:cfg.alm.max_inner
        [Lval, grad_z, info] = value_grad_z(z);
        grad_norm = norm(grad_z);
        pg_norm = projected_grad_mapping_norm(cfg, z, grad_z, N, U, cfg.alm.pg_step);

        if isfinite(prev_inner_L)
            inner_rel_alm = abs(Lval - prev_inner_L) / max(1, abs(prev_inner_L));
        else
            inner_rel_alm = inf;
        end
        if isfinite(prev_inner_crlb)
            inner_rel_crlb = abs(info.crlb - prev_inner_crlb) / max(1, abs(prev_inner_crlb));
        else
            inner_rel_crlb = inf;
        end
        if inner_rel_alm <= cfg.alm.inner_rel_tol && inner_rel_crlb <= cfg.alm.inner_crlb_rel_tol
            inner_stall_count = inner_stall_count + 1;
        else
            inner_stall_count = 0;
        end

        hist = store_proj_alm_history(hist, global_iter, outer, inner, Lval, info, ...
            grad_norm, pg_norm, rho, inner_rel_alm, inner_rel_crlb, inner_stall_count);
        hist.time_s(end+1) = toc(solver_timer);
        hist = maybe_store_iterate(hist, cfg, z, global_iter);

        if cfg.alm.verbose && (inner == 1 || mod(inner,cfg.alm.print_every)==0)
            fprintf(['  inner %03d: AL=%.3e, CRLB=%.3e, sinrViol=%.3e, ' ...
                     'powerViol=%.3e, pg=%.3e\n'], inner, Lval, info.crlb, ...
                info.max_sinr_violation, info.max_power_violation, pg_norm);
        end

        escape_resolved_this_iter = false;
        if escape_active && ...
                global_iter-last_perturb_iter >= cfg.ppfo.escape_window_iters
            escape_decrease = escape_reference_L-Lval;
            required_decrease = cfg.ppfo.escape_decrease_relative * ...
                max(1,abs(escape_reference_L));
            crlb_decrease = escape_reference_crlb-info.crlb;
            required_crlb_decrease = ...
                cfg.ppfo.escape_crlb_decrease_relative * ...
                max(abs(escape_reference_crlb),realmin);
            crlb_success = isfinite(crlb_decrease) && ...
                crlb_decrease >= required_crlb_decrease;
            allowed_escape_violation = max( ...
                cfg.ppfo.max_sinr_violation, ...
                escape_reference_sinr_violation);
            escape_feasible = info.max_sinr_violation <= ...
                allowed_escape_violation;
            escape_success = isfinite(escape_decrease) && ...
                escape_decrease >= required_decrease && escape_feasible && ...
                (~cfg.ppfo.escape_requires_crlb_improvement || crlb_success);
            hist.escape_check_iter(end+1) = global_iter;
            hist.escape_decrease(end+1) = escape_decrease;
            hist.escape_required_decrease(end+1) = required_decrease;
            hist.escape_crlb_decrease(end+1) = crlb_decrease;
            hist.escape_required_crlb_decrease(end+1) = ...
                required_crlb_decrease;
            hist.escape_crlb_success(end+1) = crlb_success;
            hist.escape_sinr_violation(end+1) = ...
                info.max_sinr_violation;
            hist.escape_feasible(end+1) = escape_feasible;
            hist.escape_success(end+1) = escape_success;
            hist.escape_reverted(end+1) = ~escape_success;
            escape_active = false;
            escape_resolved_this_iter = true;
            if ~escape_success
                z = escape_reference_z;
                if cfg.alm.verbose
                    fprintf(['  PPFO escape window ended without enough ' ...
                        'decrease; reverting to the pre-perturbation point.\n']);
                end
                break;
            end
        end

        % Only accepted/unperturbed iterates can become the returned
        % incumbent. Points inside an unresolved escape window are excluded.
        if ppfo_enabled && ~escape_active && ...
                info.max_sinr_violation <= cfg.ppfo.incumbent_feas_tol && ...
                info.max_power_violation <= cfg.ppfo.incumbent_feas_tol && ...
                info.crlb < best_feasible_crlb
            best_feasible_z = z;
            best_feasible_crlb = info.crlb;
            best_feasible_iter = global_iter;
            best_feasible_outer = outer;
            hist.incumbent_iter(end+1) = global_iter;
            hist.incumbent_outer(end+1) = outer;
            hist.incumbent_crlb(end+1) = info.crlb;
        end

        remaining_inner = cfg.alm.max_inner-inner;
        can_perturb = ppfo_enabled && ~escape_active && ...
            ~escape_resolved_this_iter && ...
            perturbation_count < cfg.ppfo.max_perturbations && ...
            outer >= cfg.ppfo.min_outer && ...
            inner >= cfg.alm.min_inner && ...
            remaining_inner >= cfg.ppfo.escape_window_iters && ...
            pg_norm <= cfg.ppfo.pg_trigger_tol && ...
            info.max_sinr_violation <= cfg.ppfo.max_sinr_violation && ...
            global_iter-last_perturb_iter > cfg.ppfo.escape_window_iters;

        if can_perturb
            escape_reference_z = z;
            escape_reference_L = Lval;
            escape_reference_crlb = info.crlb;
            escape_reference_sinr_violation = info.max_sinr_violation;
            [z,perturb_info] = perturb_power_tangent_space( ...
                cfg,z,N,U,perturb_stream);
            perturbation_count = perturbation_count+1;
            last_perturb_iter = global_iter;
            escape_active = true;

            hist.perturb_iter(end+1) = global_iter;
            hist.perturb_outer(end+1) = outer;
            hist.perturb_inner(end+1) = inner;
            hist.perturb_pg_norm(end+1) = pg_norm;
            hist.perturb_reference_alm(end+1) = Lval;
            hist.perturb_reference_crlb(end+1) = info.crlb;
            hist.perturb_reference_sinr_violation(end+1) = ...
                info.max_sinr_violation;
            hist.perturb_requested_radius(end+1) = ...
                perturb_info.requested_radius;
            hist.perturb_sampled_radius(end+1) = ...
                perturb_info.sampled_radius;
            hist.perturb_step_norm(end+1) = perturb_info.actual_step_norm;
            hist.perturb_active_blocks(end+1) = perturb_info.active_blocks;
            hist.perturb_tangent_residual(end+1) = ...
                perturb_info.max_tangent_residual;

            [Lval,grad_z,info] = value_grad_z(z);
            grad_norm = norm(grad_z);
            inner_stall_count = 0;

            if cfg.alm.verbose
                fprintf(['  PPFO perturb %d: pg=%.3e, radius=%.3e, ' ...
                    'activeAP=%d\n'],perturbation_count, ...
                    hist.perturb_pg_norm(end),perturb_info.actual_step_norm, ...
                    perturb_info.active_blocks);
            end
        elseif ~escape_active
            if inner >= cfg.alm.min_inner && pg_norm <= cfg.alm.inner_grad_tol
                break;
            end
            if inner >= cfg.alm.min_inner && ...
                    inner_stall_count >= cfg.alm.inner_patience
                break;
            end
        end

        prev_inner_L = Lval;
        prev_inner_crlb = info.crlb;

        grad_used = grad_z;
        if grad_norm > cfg.alm.grad_clip
            grad_used = grad_z * (cfg.alm.grad_clip / grad_norm);
        end
        direction = -grad_used;

        alpha = cfg.alm.inner_stepsize;
        success = false;
        zbest = z;
        Lbest = Lval;
        best_step_norm = 0;

        for bt = 1:cfg.alm.max_backtracking
            zraw = z + alpha * direction;
            Wtrial = real_vector_to_complex_matrix(zraw, N, U);
            Wtrial = project_per_ap_power_matrix(cfg, Wtrial);
            ztrial = complex_matrix_to_real_vector(Wtrial);
            s_eff = ztrial - z;
            slope_eff = grad_z' * s_eff;
            [Ltrial, ~, ~] = value_grad_z(ztrial);

            if slope_eff < 0 && isfinite(Ltrial) && Ltrial <= Lval + cfg.alm.armijo_c * slope_eff
                zbest = ztrial;
                best_step_norm = norm(s_eff) / max(1, norm(z));
                success = true;
                break;
            end
            if isfinite(Ltrial) && Ltrial < Lbest
                zbest = ztrial;
                Lbest = Ltrial;
                best_step_norm = norm(s_eff) / max(1, norm(z));
            end
            alpha = alpha * cfg.alm.backtracking_beta;
        end

        hist.step(end+1) = alpha;
        hist.step_norm(end+1) = best_step_norm;
        hist.line_search_success(end+1) = success;

        z = zbest;
        global_iter = global_iter + 1;

        if inner >= cfg.alm.min_inner && best_step_norm <= cfg.alm.inner_step_tol
            break;
        end
        if ~success && alpha < cfg.alm.min_stepsize
            break;
        end
    end

    % A perturbation is never allowed to straddle two ALM outer problems,
    % because the multiplier/rho update changes the inner objective.
    if escape_active
        z = escape_reference_z;
        hist.escape_check_iter(end+1) = global_iter;
        hist.escape_decrease(end+1) = nan;
        hist.escape_required_decrease(end+1) = ...
            cfg.ppfo.escape_decrease_relative*max(1,abs(escape_reference_L));
        hist.escape_crlb_decrease(end+1) = nan;
        hist.escape_required_crlb_decrease(end+1) = ...
            cfg.ppfo.escape_crlb_decrease_relative * ...
            max(abs(escape_reference_crlb),realmin);
        hist.escape_crlb_success(end+1) = false;
        hist.escape_sinr_violation(end+1) = nan;
        hist.escape_feasible(end+1) = false;
        hist.escape_success(end+1) = false;
        hist.escape_reverted(end+1) = true;
    end

    Wmat = real_vector_to_complex_matrix(z, N, U);
    Wmat = project_per_ap_power_matrix(cfg, Wmat);
    z = complex_matrix_to_real_vector(Wmat);

    [c_all, ~, raw] = compute_constraints_grad(cfg, channels, smat, Wmat);
    c_sinr = c_all(1:U);
    outer_sinr_violation = max(max(c_sinr,0));
    outer_power_violation = max(max(c_all(U+1:U+cfg.sys.M),0));
    outer_violation = max(outer_sinr_violation, outer_power_violation);

    mu = max(0, mu + rho * c_sinr);
    [Lcheck, grad_check_z, outinfo] = value_grad_z(z);
    outer_pg_norm = projected_grad_mapping_norm(cfg, z, grad_check_z, N, U, cfg.alm.pg_step);
    outer_grad_norm = norm(grad_check_z);
    current_crlb = outinfo.crlb;

    if ppfo_enabled && ...
            outer_sinr_violation <= cfg.ppfo.incumbent_feas_tol && ...
            outer_power_violation <= cfg.ppfo.incumbent_feas_tol && ...
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

    hist.outer_violation(end+1) = outer_violation;
    hist.outer_sinr_violation(end+1) = outer_sinr_violation;
    hist.outer_power_violation(end+1) = outer_power_violation;
    hist.outer_crlb(end+1) = current_crlb;
    hist.outer_rho(end+1) = rho;
    hist.outer_grad_norm(end+1) = outer_grad_norm;
    hist.outer_pg_norm(end+1) = outer_pg_norm;
    hist.outer_crlb_relchg(end+1) = outer_crlb_relchg;
    hist.outer_alm_obj(end+1) = Lcheck;

    stop_feas = outer_violation <= cfg.alm.feas_tol;
    stop_stat = outer_pg_norm <= cfg.alm.stationarity_tol;
    stop_crlb = outer_crlb_relchg <= cfg.alm.outer_crlb_rel_tol;
    if stop_feas && stop_crlb
        outer_stall_count = outer_stall_count + 1;
    else
        outer_stall_count = 0;
    end
    hist.outer_stall_count(end+1) = outer_stall_count;

    if cfg.alm.verbose
        fprintf(['  outer end: sinrViol=%.3e, powerViol=%.3e, pg=%.3e, ' ...
                 'CRLB=%.3e, relCRLB=%.3e, minSINR=%.2f dB, maxPower=%.2f dBm\n'], ...
            outer_sinr_violation, outer_power_violation, outer_pg_norm, current_crlb, ...
            outer_crlb_relchg, pow2db(max(raw.min_sinr,realmin)), w2dbm(max(raw.max_ap_power,realmin)));
    end

    if outer >= cfg.alm.min_outer && stop_feas && stop_stat
        sol.success = true;
        sol.message = 'feasibility + projected stationarity reached';
        break;
    end
    if outer >= cfg.alm.min_outer && stop_feas && ~cfg.alm.require_stationarity_for_stop && ...
            outer_stall_count >= cfg.alm.outer_patience
        sol.success = true;
        sol.message = 'feasibility + CRLB plateau reached';
        break;
    end

    if outer_sinr_violation > cfg.alm.feas_tol
        if outer == 1 || ~isfinite(prev_outer_sinr_violation)
            rho = min(cfg.alm.rho_growth * rho, cfg.alm.rho_max);
        elseif outer_sinr_violation > cfg.alm.rho_update_threshold * prev_outer_sinr_violation
            rho = min(cfg.alm.rho_growth * rho, cfg.alm.rho_max);
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
if ppfo_enabled && ~isempty(best_feasible_z)
    Wbest = real_vector_to_complex_matrix(best_feasible_z,N,U);
    Wbest = project_per_ap_power_matrix(cfg,Wbest);
    wbest = matrix_to_wcell(cfg,Wbest);
    metrics_best = evaluate_beamforming_solution( ...
        cfg,channels,smat,wbest);
    current_feasible = all(metrics.sinr >= cfg.comm.Gamma * ...
        (1-cfg.ppfo.incumbent_feas_tol)) && ...
        all(metrics.ap_power <= cfg.power.Pmax * ...
        (1+cfg.ppfo.incumbent_feas_tol));
    best_is_feasible = all(metrics_best.sinr >= cfg.comm.Gamma * ...
        (1-cfg.ppfo.incumbent_feas_tol)) && ...
        all(metrics_best.ap_power <= cfg.power.Pmax * ...
        (1+cfg.ppfo.incumbent_feas_tol));
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
sol.initial_Wmat = initial_Wmat;
sol.initial_wcell = matrix_to_wcell(cfg,initial_Wmat);
sol.initial_z = complex_matrix_to_real_vector(initial_Wmat);
sol.caller_supplied_initializer = caller_supplied_initializer;
sol.mu = mu;
sol.rho = rho;
sol.obj_scale = obj_scale;
sol.metrics = metrics;
sol.hist = hist;
sol.ppfo_enabled = ppfo_enabled;
sol.perturbation_count = perturbation_count;
sol.returned_best_incumbent = returned_best_incumbent;
sol.best_feasible_crlb = best_feasible_crlb;
sol.best_feasible_iter = best_feasible_iter;
sol.best_feasible_outer = best_feasible_outer;
hist.initial_Wmat = initial_Wmat;
hist.caller_supplied_initializer = caller_supplied_initializer;
sol.hist = hist;

    function [Lval, grad_z, info] = value_grad_z(zz)
        Wq = real_vector_to_complex_matrix(zz, N, U);
        [Lval, grad_complex, info] = sinr_augmented_lagrangian_value_grad( ...
            cfg, channels, smat, Wq, mu, rho, obj_scale);
        grad_z = complex_matrix_to_real_vector(grad_complex);
    end
end

function hist = init_proj_alm_history()
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
hist.inner_rel_alm = [];
hist.inner_rel_crlb = [];
hist.inner_stall_count = [];
hist.step = [];
hist.step_norm = [];
hist.line_search_success = [];
hist.outer_violation = [];
hist.outer_sinr_violation = [];
hist.outer_power_violation = [];
hist.outer_crlb = [];
hist.outer_rho = [];
hist.outer_grad_norm = [];
hist.outer_pg_norm = [];
hist.outer_crlb_relchg = [];
hist.outer_alm_obj = [];
hist.outer_stall_count = [];
hist.snapshot_iter = [];
hist.z_snapshots = [];
hist.perturb_iter = [];
hist.perturb_outer = [];
hist.perturb_inner = [];
hist.perturb_pg_norm = [];
hist.perturb_reference_alm = [];
hist.perturb_reference_crlb = [];
hist.perturb_reference_sinr_violation = [];
hist.perturb_requested_radius = [];
hist.perturb_sampled_radius = [];
hist.perturb_step_norm = [];
hist.perturb_active_blocks = [];
hist.perturb_tangent_residual = [];
hist.escape_check_iter = [];
hist.escape_decrease = [];
hist.escape_required_decrease = [];
hist.escape_crlb_decrease = [];
hist.escape_required_crlb_decrease = [];
hist.escape_crlb_success = [];
hist.escape_sinr_violation = [];
hist.escape_feasible = [];
hist.escape_success = [];
hist.escape_reverted = [];
hist.incumbent_iter = [];
hist.incumbent_outer = [];
hist.incumbent_crlb = [];
end

function hist = store_proj_alm_history(hist, iter, outer, inner, Lval, info, grad_norm, pg_norm, rho, inner_rel_alm, inner_rel_crlb, inner_stall_count)
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
hist.inner_rel_alm(end+1) = inner_rel_alm;
hist.inner_rel_crlb(end+1) = inner_rel_crlb;
hist.inner_stall_count(end+1) = inner_stall_count;
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

function cfg = fill_proj_alm_defaults(cfg)
if ~isfield(cfg,'alm'); cfg.alm = struct(); end
if ~isfield(cfg.alm,'max_outer'); cfg.alm.max_outer = 20; end
if ~isfield(cfg.alm,'max_inner'); cfg.alm.max_inner = 100; end
if ~isfield(cfg.alm,'min_outer'); cfg.alm.min_outer = 3; end
if ~isfield(cfg.alm,'min_inner'); cfg.alm.min_inner = 5; end
if ~isfield(cfg.alm,'rho0'); cfg.alm.rho0 = 1; end
if ~isfield(cfg.alm,'rho_growth'); cfg.alm.rho_growth = 3; end
if ~isfield(cfg.alm,'rho_max'); cfg.alm.rho_max = 1e10; end
if ~isfield(cfg.alm,'rho_update_threshold'); cfg.alm.rho_update_threshold = 0.9; end
if ~isfield(cfg.alm,'feas_tol'); cfg.alm.feas_tol = 1e-4; end
if ~isfield(cfg.alm,'stationarity_tol'); cfg.alm.stationarity_tol = 1e-3; end
if ~isfield(cfg.alm,'outer_crlb_rel_tol'); cfg.alm.outer_crlb_rel_tol = 1e-5; end
if ~isfield(cfg.alm,'outer_patience'); cfg.alm.outer_patience = 2; end
if ~isfield(cfg.alm,'require_stationarity_for_stop'); cfg.alm.require_stationarity_for_stop = false; end
if ~isfield(cfg.alm,'inner_grad_tol'); cfg.alm.inner_grad_tol = 1e-6; end
if ~isfield(cfg.alm,'inner_rel_tol'); cfg.alm.inner_rel_tol = 1e-6; end
if ~isfield(cfg.alm,'inner_crlb_rel_tol'); cfg.alm.inner_crlb_rel_tol = 1e-6; end
if ~isfield(cfg.alm,'inner_patience'); cfg.alm.inner_patience = 5; end
if ~isfield(cfg.alm,'inner_step_tol'); cfg.alm.inner_step_tol = 1e-10; end
if ~isfield(cfg.alm,'inner_stepsize'); cfg.alm.inner_stepsize = 1e-3; end
if ~isfield(cfg.alm,'min_stepsize'); cfg.alm.min_stepsize = 1e-12; end
if ~isfield(cfg.alm,'max_backtracking'); cfg.alm.max_backtracking = 25; end
if ~isfield(cfg.alm,'backtracking_beta'); cfg.alm.backtracking_beta = 0.5; end
if ~isfield(cfg.alm,'armijo_c'); cfg.alm.armijo_c = 1e-4; end
if ~isfield(cfg.alm,'grad_clip'); cfg.alm.grad_clip = 1e3; end
if ~isfield(cfg.alm,'objective_type'); cfg.alm.objective_type = 'crlb_ratio'; end
if ~isfield(cfg.alm,'objective_scale_mode'); cfg.alm.objective_scale_mode = 'initial'; end
if ~isfield(cfg.alm,'det_floor'); cfg.alm.det_floor = 1e-30; end
if ~isfield(cfg.alm,'init_mode'); cfg.alm.init_mode = 'zf'; end
if ~isfield(cfg.alm,'pg_step'); cfg.alm.pg_step = 1e-3; end
if ~isfield(cfg.alm,'verbose'); cfg.alm.verbose = true; end
if ~isfield(cfg.alm,'print_every'); cfg.alm.print_every = 20; end
if ~isfield(cfg,'ppfo'); cfg.ppfo = struct(); end
if ~isfield(cfg.ppfo,'enabled'); cfg.ppfo.enabled = false; end
if ~isfield(cfg.ppfo,'seed')
    if isfield(cfg,'seed'); cfg.ppfo.seed = cfg.seed+5000;
    else; cfg.ppfo.seed = 5000;
    end
end
if ~isfield(cfg.ppfo,'pg_trigger_tol'); cfg.ppfo.pg_trigger_tol = 0.88; end
if ~isfield(cfg.ppfo,'perturb_radius_relative')
    cfg.ppfo.perturb_radius_relative = 2e-3;
end
if ~isfield(cfg.ppfo,'min_perturb_radius')
    cfg.ppfo.min_perturb_radius = 1e-8;
end
if ~isfield(cfg.ppfo,'power_active_tol'); cfg.ppfo.power_active_tol = 1e-6; end
if ~isfield(cfg.ppfo,'escape_window_iters')
    cfg.ppfo.escape_window_iters = 20;
end
cfg.ppfo.escape_window_iters = max(1,round(cfg.ppfo.escape_window_iters));
if ~isfield(cfg.ppfo,'escape_decrease_relative')
    cfg.ppfo.escape_decrease_relative = 1e-5;
end
if ~isfield(cfg.ppfo,'escape_requires_crlb_improvement')
    cfg.ppfo.escape_requires_crlb_improvement = true;
end
if ~isfield(cfg.ppfo,'escape_crlb_decrease_relative')
    cfg.ppfo.escape_crlb_decrease_relative = 0;
end
if ~isfield(cfg.ppfo,'max_perturbations'); cfg.ppfo.max_perturbations = 4; end
cfg.ppfo.max_perturbations = max(0,round(cfg.ppfo.max_perturbations));
if ~isfield(cfg.ppfo,'min_outer'); cfg.ppfo.min_outer = 3; end
if ~isfield(cfg.ppfo,'max_sinr_violation')
    cfg.ppfo.max_sinr_violation = 1e-2;
end
if ~isfield(cfg.ppfo,'outer_patience'); cfg.ppfo.outer_patience = 4; end
if ~isfield(cfg.ppfo,'incumbent_feas_tol')
    cfg.ppfo.incumbent_feas_tol = 1e-4;
end
end

function wcell = random_initial_beamformers_projected(cfg)
M = cfg.sys.M; U = cfg.sys.U; Nt = cfg.sys.Nt; Pmax = cfg.power.Pmax;
wcell = cell(M,U);
for m = 1:M
    X = (randn(Nt,U)+1j*randn(Nt,U))/sqrt(2);
    X = X / max(norm(X,'fro'), realmin) * sqrt(0.5*Pmax);
    for u = 1:U
        wcell{m,u} = X(:,u);
    end
end
end
