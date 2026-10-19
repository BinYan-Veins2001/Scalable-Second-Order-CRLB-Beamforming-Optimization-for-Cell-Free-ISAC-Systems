function [wcell, hist] = restore_sinr_near_reference(cfg, channels, smat, ref_wcell, opts)
%RESTORE_SINR_NEAR_REFERENCE Projected feasibility restoration near a reference.
%
% The routine minimizes a simple smooth penalty
%
%   0.5 * sum_u max(0, c_u(w))^2 + lambda/2 * ||w-ref||^2
%
% where c_u are the normalized SINR inequalities from compute_constraints_grad.
% Per-AP power constraints are enforced by projection.

if nargin < 5
    opts = struct();
end
if ~isfield(opts,'lambda'); opts.lambda = 1e-3; end
if ~isfield(opts,'max_iter'); opts.max_iter = 800; end
if ~isfield(opts,'alpha0'); opts.alpha0 = 1e-3; end
if ~isfield(opts,'beta'); opts.beta = 0.5; end
if ~isfield(opts,'armijo_c'); opts.armijo_c = 1e-4; end
if ~isfield(opts,'max_backtracking'); opts.max_backtracking = 25; end
if ~isfield(opts,'feas_tol'); opts.feas_tol = 1e-4; end
if ~isfield(opts,'stall_patience'); opts.stall_patience = 30; end
if ~isfield(opts,'grad_clip'); opts.grad_clip = 1e3; end
if ~isfield(opts,'verbose'); opts.verbose = false; end

U = cfg.sys.U;
N = cfg.sys.M * cfg.sys.Nt;

Wref = project_per_ap_power_matrix(cfg, vectorize_wcell_to_matrix(cfg, ref_wcell));
zref = complex_matrix_to_real_vector(Wref);
z = zref;

hist.iter = [];
hist.restore_obj = [];
hist.max_sinr_violation = [];
hist.min_sinr_dB = [];
hist.crlb = [];
hist.step_norm = [];

best_z = z;
best_score = inf;
stall = 0;
prev_obj = inf;

for it = 1:opts.max_iter
    [fval, grad_z, info] = value_grad(z);
    grad_norm = norm(grad_z);

    hist.iter(end+1) = it;
    hist.restore_obj(end+1) = fval;
    hist.max_sinr_violation(end+1) = info.max_sinr_violation;
    hist.min_sinr_dB(end+1) = info.min_sinr_dB;
    hist.crlb(end+1) = info.crlb;

    score = info.max_sinr_violation * 1e6 + info.crlb;
    if score < best_score
        best_score = score;
        best_z = z;
    end

    if opts.verbose && (it == 1 || mod(it,50) == 0)
        fprintf('  restore %04d: obj=%.3e, sinrViol=%.3e, CRLB=%.3e, minSINR=%.2f dB\n', ...
            it, fval, info.max_sinr_violation, info.crlb, info.min_sinr_dB);
    end

    if info.max_sinr_violation <= opts.feas_tol && abs(fval - prev_obj) / max(1, abs(prev_obj)) <= 1e-8
        stall = stall + 1;
    else
        stall = 0;
    end
    if stall >= opts.stall_patience
        break;
    end

    grad_used = grad_z;
    if grad_norm > opts.grad_clip
        grad_used = grad_z * (opts.grad_clip / grad_norm);
    end
    direction = -grad_used;

    alpha = opts.alpha0;
    zbest_line = z;
    fbest_line = fval;
    step_norm = 0;

    for bt = 1:opts.max_backtracking
        ztrial_raw = z + alpha * direction;
        Wtrial = real_vector_to_complex_matrix(ztrial_raw, N, U);
        Wtrial = project_per_ap_power_matrix(cfg, Wtrial);
        ztrial = complex_matrix_to_real_vector(Wtrial);
        strial = ztrial - z;
        [ftrial, ~, ~] = value_grad(ztrial);
        slope = grad_z' * strial;
        if slope < 0 && isfinite(ftrial) && ftrial <= fval + opts.armijo_c * slope
            zbest_line = ztrial;
            fbest_line = ftrial;
            step_norm = norm(strial) / max(1, norm(z));
            break;
        end
        if isfinite(ftrial) && ftrial < fbest_line
            zbest_line = ztrial;
            fbest_line = ftrial;
            step_norm = norm(strial) / max(1, norm(z));
        end
        alpha = alpha * opts.beta;
    end

    hist.step_norm(end+1) = step_norm;
    z = zbest_line;
    prev_obj = fval;
end

% Prefer the final point if feasible; otherwise return the best restoration score.
[~,~,final_info] = value_grad(z);
if final_info.max_sinr_violation > opts.feas_tol
    z = best_z;
end

Wmat = real_vector_to_complex_matrix(z, N, U);
Wmat = project_per_ap_power_matrix(cfg, Wmat);
wcell = matrix_to_wcell(cfg, Wmat);

    function [fval, grad_z, info] = value_grad(zz)
        Wmat = real_vector_to_complex_matrix(zz, N, U);
        [c, G, raw] = compute_constraints_grad(cfg, channels, smat, Wmat);
        c_sinr = c(1:U);
        active = c_sinr > 0;
        fval = 0.5 * sum(c_sinr(active).^2);
        grad = zeros(size(Wmat));
        for uu = 1:U
            if active(uu)
                grad = grad + c_sinr(uu) * G(:,:,uu);
            end
        end
        dz_reg = zz - zref;
        fval = fval + 0.5 * opts.lambda * (dz_reg' * dz_reg);
        grad_z = complex_matrix_to_real_vector(grad) + opts.lambda * dz_reg;

        metrics = evaluate_beamforming_solution(cfg, channels, smat, matrix_to_wcell(cfg, Wmat));
        info.max_sinr_violation = max(max(c_sinr,0));
        info.min_sinr_dB = pow2db(max(raw.min_sinr, realmin));
        info.crlb = metrics.crlb;
    end
end
