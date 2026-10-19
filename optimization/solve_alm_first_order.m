function [sol, hist] = solve_alm_first_order(cfg, channels, smat, init_wcell)
%SOLVE_ALM_FIRST_ORDER Real-domain first-order ALM baseline without saddle escaping.
%
% This baseline directly optimizes the original complex beamforming vectors,
% but the descent variable is explicitly represented in the real domain:
%
%   z = [real(vec(Wmat)); imag(vec(Wmat))] in R^{2 M Nt U}.
%
% It does not use SDR, rank-one lifting, Hessian information, random
% perturbation, negative-curvature directions, or trust-region steps.
%
% Outer loop: inequality augmented-Lagrangian multiplier/rho update.
% Inner loop: real-domain backtracking gradient descent.

cfg = fill_alm_defaults(cfg);

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M * Nt;
J = U + M;

if nargin < 4 || isempty(init_wcell)
    switch lower(cfg.alm.init_mode)
        case 'identity'
            init_wcell = initialize_beamformer_identity(cfg);
        case 'zf'
            init_wcell = build_shared_theorem_aware_zf_initializer( ...
                cfg,channels,smat);
        case 'random'
            init_wcell = random_initial_beamformers(cfg);
        otherwise
            error('Unknown ALM init_mode: %s', cfg.alm.init_mode);
    end
end

Wmat0 = vectorize_wcell_to_matrix(cfg, init_wcell);
z = complex_matrix_to_real_vector(Wmat0);

% Objective scaling improves numerical balance between CRLB ratio and constraints.
% For the log-ratio objective, the logarithm already compresses scale, so
% the default behavior is no additional initial-value scaling.
[F0, ~] = compute_crlb_objective_grad(cfg, smat, Wmat0);

if strcmpi(cfg.alm.objective_scale_mode, 'initial') && ...
        ~strcmpi(cfg.alm.objective_type, 'log_crlb_ratio')
    obj_scale = max(abs(F0), 1e-12);
else
    obj_scale = 1;
end

mu = zeros(J, 1);
rho = cfg.alm.rho0;

prev_outer_violation = inf;
prev_outer_crlb = inf;

hist = init_alm_history();
hist.init_mode = cfg.alm.init_mode;
hist.domain = 'real';
hist.objective_type = cfg.alm.objective_type;

sol.success = false;
sol.message = 'max_outer reached';

global_iter = 0;

if cfg.alm.verbose
    fprintf('[FO-ALM] init_mode=%s, objective=%s, descent_domain=real, dim(z)=%d\n', ...
        cfg.alm.init_mode, cfg.alm.objective_type, numel(z));

    % Explicitly print the initial point. This is the same point stored as
    % hist.iter=0 in the first inner-loop record, but the separate line makes
    % cross-algorithm initialization comparisons unambiguous.
    [Linit, grad_init_z, info_init, ~] = value_grad_real( ...
        cfg, channels, smat, z, mu, rho, obj_scale, N, U);
    fprintf(['[FO-ALM init k=0] AL=%.3e, CRLB=%.3e, viol=%.3e, ' ...
             'grad_z=%.3e, minSINR=%.2f dB, maxPower=%.2f dBm\n'], ...
        Linit, info_init.crlb, info_init.max_violation, norm(grad_init_z), ...
        pow2db(max(info_init.raw.min_sinr, realmin)), ...
        w2dbm(max(info_init.raw.max_ap_power, realmin)));
end

for outer = 1:cfg.alm.max_outer

    if cfg.alm.verbose
        fprintf('\nALM outer %02d/%02d, rho=%.3e\n', ...
            outer, cfg.alm.max_outer, rho);
    end

    for inner = 1:cfg.alm.max_inner

        [Lval, grad_z, info, Wmat] = value_grad_real( ...
            cfg, channels, smat, z, mu, rho, obj_scale, N, U);

        grad_norm = norm(grad_z);

        hist = store_alm_history( ...
            hist, global_iter, outer, inner, Lval, info, grad_norm, rho);
        hist = maybe_store_iterate(hist, cfg, z, global_iter);

        if cfg.alm.verbose && ((inner == 1 && global_iter > 0) || ...
                mod(inner, cfg.alm.print_every) == 0)
            fprintf(['  inner %03d: AL=%.3e, CRLB=%.3e, viol=%.3e, ' ...
                     'grad_z=%.3e, minSINR=%.2f dB\n'], ...
                inner, Lval, info.crlb, info.max_violation, ...
                grad_norm, pow2db(max(info.raw.min_sinr, realmin)));
        end

        if grad_norm <= cfg.alm.inner_grad_tol
            hist.step(end+1) = 0;
            hist.line_search_success(end+1) = true;

            if cfg.alm.verbose
                fprintf('  inner stopped by gradient tolerance at inner %d.\n', inner);
            end

            break;
        end

        grad_used = grad_z;

        if grad_norm > cfg.alm.grad_clip
            grad_used = grad_z * (cfg.alm.grad_clip / grad_norm);
        end

        direction = -grad_used;
        slope = grad_z' * direction;

        if slope >= 0
            direction = -grad_z;
            slope = -grad_norm^2;
        end

        alpha = cfg.alm.inner_stepsize;
        success = false;

        zbest = z;
        Lbest = Lval;

        for bt = 1:cfg.alm.max_backtracking

            ztrial = z + alpha * direction;

            Wtrial = real_vector_to_complex_matrix(ztrial, N, U);
            Wtrial = apply_power_safeguard(cfg, smat, Wtrial);
            ztrial = complex_matrix_to_real_vector(Wtrial);

            [Ltrial, ~, ~] = value_grad_real( ...
                cfg, channels, smat, ztrial, mu, rho, obj_scale, N, U);

            if isfinite(Ltrial) && ...
                    Ltrial <= Lval + cfg.alm.armijo_c * alpha * slope
                zbest = ztrial;
                Lbest = Ltrial;
                success = true;
                break;
            end

            if isfinite(Ltrial) && Ltrial < Lbest
                zbest = ztrial;
                Lbest = Ltrial;
            end

            alpha = alpha * cfg.alm.backtracking_beta;
        end

        hist.step(end+1) = alpha;
        hist.line_search_success(end+1) = success;

        z = zbest;
        global_iter = global_iter + 1;

        if ~success && alpha < cfg.alm.min_stepsize
            if cfg.alm.verbose
                fprintf('  line search stalled at inner %d.\n', inner);
            end
            break;
        end
    end

    % ==============================================================
    % Outer multiplier/rho update and stopping check
    % ==============================================================

    Wmat = real_vector_to_complex_matrix(z, N, U);

    [c, ~, raw] = compute_constraints_grad(cfg, channels, smat, Wmat);

    % Current normalized inequality violation.
    outer_viol = max(max(c, 0));

    % Multiplier update.
    mu = max(0, mu + rho * c);

    % Check ALM stationarity after multiplier update.
    [Lcheck, grad_check_complex, outinfo] = augmented_lagrangian_value_grad( ...
        cfg, channels, smat, Wmat, mu, rho, obj_scale);

    grad_check_z = complex_matrix_to_real_vector(grad_check_complex);
    outer_grad_norm = norm(grad_check_z);

    current_crlb = outinfo.crlb;

    if isfinite(prev_outer_crlb)
        outer_crlb_relchg = abs(current_crlb - prev_outer_crlb) / ...
            max(1, abs(prev_outer_crlb));
    else
        outer_crlb_relchg = inf;
    end

    hist.outer_violation(end+1) = outer_viol;
    hist.outer_crlb(end+1) = current_crlb;
    hist.outer_rho(end+1) = rho;
    hist.outer_grad_norm(end+1) = outer_grad_norm;
    hist.outer_crlb_relchg(end+1) = outer_crlb_relchg;
    hist.outer_alm_obj(end+1) = Lcheck;

    if cfg.alm.verbose
        fprintf(['  outer end: viol=%.3e, grad=%.3e, CRLB=%.3e, ' ...
                 'CRLB-relchg=%.3e, minSINR=%.2f dB, maxPower=%.2f dBm\n'], ...
            outer_viol, outer_grad_norm, current_crlb, outer_crlb_relchg, ...
            pow2db(max(raw.min_sinr, realmin)), ...
            w2dbm(max(raw.max_ap_power, realmin)));
    end

    % Correct stopping rule:
    % Do not stop merely because the point is feasible.
    stop_feas = outer_viol <= cfg.alm.feas_tol;
    stop_stat = outer_grad_norm <= cfg.alm.stationarity_tol;
    stop_crlb = outer_crlb_relchg <= cfg.alm.outer_crlb_rel_tol;

    if outer >= cfg.alm.min_outer && stop_feas && stop_stat && stop_crlb
        sol.success = true;
        sol.message = 'feasibility + stationarity + CRLB stability reached';
        break;
    end

    % Optional legacy behavior. Disabled by default.
    if cfg.alm.stop_on_feasible_only && outer >= cfg.alm.min_outer && stop_feas
        sol.success = true;
        sol.message = 'feasibility tolerance reached';
        break;
    end

    % Penalty update:
    % Increase rho only if feasibility is still poor or not improving.
    if outer_viol > cfg.alm.feas_tol
        if outer == 1 || ~isfinite(prev_outer_violation)
            rho = min(cfg.alm.rho_growth * rho, cfg.alm.rho_max);
        elseif outer_viol > cfg.alm.rho_update_threshold * prev_outer_violation
            rho = min(cfg.alm.rho_growth * rho, cfg.alm.rho_max);
        end
    end

    prev_outer_violation = outer_viol;
    prev_outer_crlb = current_crlb;
end

Wmat = real_vector_to_complex_matrix(z, N, U);
wcell = matrix_to_wcell(cfg, Wmat);
metrics = evaluate_beamforming_solution(cfg, channels, smat, wcell);

sol.wcell = wcell;
sol.Wmat = Wmat;
sol.z = z;
sol.mu = mu;
sol.rho = rho;
sol.obj_scale = obj_scale;
sol.metrics = metrics;
sol.hist = hist;

end


function [Lval, grad_z, info, Wmat] = value_grad_real( ...
    cfg, channels, smat, z, mu, rho, obj_scale, N, U)

Wmat = real_vector_to_complex_matrix(z, N, U);

[Lval, grad_complex, info] = augmented_lagrangian_value_grad( ...
    cfg, channels, smat, Wmat, mu, rho, obj_scale);

grad_z = complex_matrix_to_real_vector(grad_complex);

end


function hist = init_alm_history()

hist.iter = [];
hist.outer = [];
hist.inner = [];

hist.alm_obj = [];
hist.crlb = [];
hist.crlb_ratio = [];
hist.optim_obj = [];
hist.max_violation = [];
hist.min_sinr_dB = [];
hist.max_ap_power_dBm = [];
hist.grad_norm = [];
hist.rho = [];

hist.step = [];
hist.line_search_success = [];

hist.outer_violation = [];
hist.outer_crlb = [];
hist.outer_rho = [];
hist.outer_grad_norm = [];
hist.outer_crlb_relchg = [];
hist.outer_alm_obj = [];
hist.snapshot_iter = [];
hist.z_snapshots = [];

end


function hist = store_alm_history( ...
    hist, iter, outer, inner, Lval, info, grad_norm, rho)

hist.iter(end+1) = iter;
hist.outer(end+1) = outer;
hist.inner(end+1) = inner;

hist.alm_obj(end+1) = Lval;
hist.crlb(end+1) = info.crlb;
hist.crlb_ratio(end+1) = info.ratio;
hist.optim_obj(end+1) = info.optim_obj;
hist.max_violation(end+1) = info.max_violation;
hist.min_sinr_dB(end+1) = pow2db(max(info.raw.min_sinr, realmin));
hist.max_ap_power_dBm(end+1) = w2dbm(max(info.raw.max_ap_power, realmin));
hist.grad_norm(end+1) = grad_norm;
hist.rho(end+1) = rho;

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


function cfg = fill_alm_defaults(cfg)

if ~isfield(cfg, 'alm')
    cfg.alm = struct();
end

if ~isfield(cfg.alm, 'max_outer')
    cfg.alm.max_outer = 20;
end

if ~isfield(cfg.alm, 'max_inner')
    cfg.alm.max_inner = 100;
end

if ~isfield(cfg.alm, 'rho0')
    cfg.alm.rho0 = 1;
end

if ~isfield(cfg.alm, 'rho_growth')
    cfg.alm.rho_growth = 3;
end

if ~isfield(cfg.alm, 'rho_max')
    cfg.alm.rho_max = 1e6;
end

if ~isfield(cfg.alm, 'rho_update_threshold')
    cfg.alm.rho_update_threshold = 0.9;
end

if ~isfield(cfg.alm, 'feas_tol')
    cfg.alm.feas_tol = 1e-4;
end

if ~isfield(cfg.alm, 'stationarity_tol')
    cfg.alm.stationarity_tol = 1e-3;
end

if ~isfield(cfg.alm, 'outer_crlb_rel_tol')
    cfg.alm.outer_crlb_rel_tol = 1e-5;
end

if ~isfield(cfg.alm, 'min_outer')
    cfg.alm.min_outer = 3;
end

if ~isfield(cfg.alm, 'stop_on_feasible_only')
    cfg.alm.stop_on_feasible_only = false;
end

if ~isfield(cfg.alm, 'inner_grad_tol')
    cfg.alm.inner_grad_tol = 1e-6;
end

if ~isfield(cfg.alm, 'inner_stepsize')
    cfg.alm.inner_stepsize = 1e-3;
end

if ~isfield(cfg.alm, 'min_stepsize')
    cfg.alm.min_stepsize = 1e-12;
end

if ~isfield(cfg.alm, 'max_backtracking')
    cfg.alm.max_backtracking = 25;
end

if ~isfield(cfg.alm, 'backtracking_beta')
    cfg.alm.backtracking_beta = 0.5;
end

if ~isfield(cfg.alm, 'armijo_c')
    cfg.alm.armijo_c = 1e-4;
end

if ~isfield(cfg.alm, 'grad_clip')
    cfg.alm.grad_clip = 1e3;
end

if ~isfield(cfg.alm, 'objective_type')
    cfg.alm.objective_type = 'crlb_ratio';
end

if ~isfield(cfg.alm, 'objective_scale_mode')
    cfg.alm.objective_scale_mode = 'initial';
end

if ~isfield(cfg.alm, 'init_mode')
    cfg.alm.init_mode = 'zf';
end

if ~isfield(cfg.alm, 'identity_power_scale')
    cfg.alm.identity_power_scale = cfg.opt.init_power_scale;
end

if ~isfield(cfg.alm, 'real_domain')
    cfg.alm.real_domain = true;
end

if ~isfield(cfg.alm, 'det_floor')
    cfg.alm.det_floor = 1e-30;
end

if ~isfield(cfg.alm, 'use_power_safeguard')
    cfg.alm.use_power_safeguard = false;
end

if ~isfield(cfg.alm, 'power_safeguard_factor')
    cfg.alm.power_safeguard_factor = 5;
end

if ~isfield(cfg.alm, 'verbose')
    cfg.alm.verbose = true;
end

if ~isfield(cfg.alm, 'print_every')
    cfg.alm.print_every = 20;
end

end


function wcell = random_initial_beamformers(cfg)

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
Pmax = cfg.power.Pmax;

wcell = cell(M, U);

for m = 1:M
    X = (randn(Nt, U) + 1j * randn(Nt, U)) / sqrt(2);
    X = X / max(norm(X, 'fro'), realmin) * sqrt(0.5 * Pmax);

    for u = 1:U
        wcell{m,u} = X(:,u);
    end
end

end
