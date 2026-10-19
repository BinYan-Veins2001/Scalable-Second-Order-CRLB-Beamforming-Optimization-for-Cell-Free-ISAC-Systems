function cert = evaluate_as_rtr_initialization_certificate( ...
    cfg,channels,smat,init_wcell,opts)
%EVALUATE_AS_RTR_INITIALIZATION_CERTIFICATE Quantify AS-RTR recoverability.
%
% The supplied beamformer is first projected onto the per-AP power balls,
% exactly as in POLISH_BEAMFORMER_ACTIVE_RTR.  At that entry point this
% routine computes
%   (i)  the minimum-norm step into the linearized SINR feasible set;
%   (ii) the Newton-like normal step actually used by AS RTR;
%   (iii) trust-region reachability indices; and
%   (iv) a conservative nonlinear one-step sufficient certificate.
%
% The optional pure-restoration experiment applies only the AS-RTR SINR
% normal correction, without an objective/tCG step.  It is diagnostic and
% does not alter the supplied point or the AS-RTR solver.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(opts);

N = cfg.sys.M*cfg.sys.Nt;
U = cfg.sys.U;
dim = 2*N*U;
Wentry = project_per_ap_power_matrix( ...
    cfg,vectorize_wcell_to_matrix(cfg,init_wcell));
zentry = complex_matrix_to_real_vector(Wentry);
[centry,Gentry,raw_entry] = constraint_derivatives(zentry);
cbase = centry;

[A,idx,Qpower] = correction_jacobian(centry,Gentry,cbase);
residual = centry(idx)+opts.sinr_retraction_margin;
[step_linear,linear_qp_success,linear_qp_kkt] = ...
    minimum_norm_inequality_step(A,residual,opts.rank_tol);

if isempty(idx)
    step_code = zeros(dim,1);
    step_newton = zeros(dim,1);
else
    normal_matrix = A'*A + ...
        opts.sinr_retraction_regularization*eye(numel(idx));
    step_code = -A*(normal_matrix\residual);
    step_code = step_code-Qpower*(Qpower'*step_code);
    step_newton = -pinv(A',opts.rank_tol)*residual;
    step_newton = step_newton-Qpower*(Qpower'*step_newton);
end

linearized_code_residual = residual+A'*step_code;
linearized_qp_residual = residual+A'*step_linear;
singular_values = svd(A,'econ');
positive_singular_values = singular_values( ...
    singular_values > opts.rank_tol*max([1;singular_values]));
if isempty(idx)
    sigma_min = inf;
    sigma_max = 0;
    jacobian_condition = 1;
    jacobian_rank = 0;
elseif isempty(positive_singular_values)
    sigma_min = 0;
    sigma_max = 0;
    jacobian_condition = inf;
    jacobian_rank = 0;
else
    sigma_min = min(positive_singular_values);
    sigma_max = max(positive_singular_values);
    jacobian_condition = sigma_max/sigma_min;
    jacobian_rank = numel(positive_singular_values);
end

step_linear_norm = norm(step_linear);
step_code_norm = norm(step_code);
step_newton_norm = norm(step_newton);
regularized_step_relative_error = norm(step_code-step_newton) / ...
    max(step_newton_norm,realmin);
single_radius = opts.sinr_retraction_max_fraction*opts.delta0;
total_radius = opts.sinr_retraction_max_iter*single_radius;

% The SINR constraints are quadratic in the real beamforming vector.  Their
% real Hessian spectral norms therefore give exact global Jacobian-Lipschitz
% constants in the ambient space.
constraint_hessian_norm = sinr_hessian_norms(idx);
if isempty(constraint_hessian_norm)
    vector_jacobian_lipschitz = 0;
else
    vector_jacobian_lipschitz = norm(constraint_hessian_norm);
end
if isfinite(sigma_min) && sigma_min > 0
    kantorovich_indicator = vector_jacobian_lipschitz * ...
        step_newton_norm/sigma_min;
else
    kantorovich_indicator = inf;
end
if isfinite(kantorovich_indicator) && kantorovich_indicator > 0 && ...
        kantorovich_indicator <= 0.5
    kantorovich_radius_bound = ...
        (1-sqrt(1-2*kantorovich_indicator))*step_newton_norm / ...
        kantorovich_indicator;
elseif kantorovich_indicator == 0
    kantorovich_radius_bound = step_newton_norm;
else
    kantorovich_radius_bound = inf;
end
entry_power_radius_slack = inf;
for power_slack_idx = 1:cfg.sys.M
    power_slack_rows = ...
        (power_slack_idx-1)*cfg.sys.Nt+(1:cfg.sys.Nt);
    entry_power_radius_slack = min(entry_power_radius_slack, ...
        sqrt(cfg.power.Pmax)-norm(Wentry(power_slack_rows,:),'fro'));
end
entry_power_radius_slack = max(entry_power_radius_slack,0);
full_row_rank = jacobian_rank == numel(idx);
inactive_sinr = setdiff((1:U)',idx(:));
inactive_sinr_upper_bound = zeros(numel(inactive_sinr),1);
if isfinite(kantorovich_radius_bound)
    inactive_hessian_norm = sinr_hessian_norms(inactive_sinr);
    for inactive_idx = 1:numel(inactive_sinr)
        user_idx = inactive_sinr(inactive_idx);
        inactive_sinr_upper_bound(inactive_idx) = centry(user_idx) + ...
            norm(Gentry(:,user_idx))*kantorovich_radius_bound + ...
            0.5*inactive_hessian_norm(inactive_idx)* ...
            kantorovich_radius_bound^2;
    end
end
inactive_sinr_safe = all(inactive_sinr_upper_bound <= opts.feas_tol);
kantorovich_tr_certificate = ~isempty(idx) && full_row_rank && ...
    kantorovich_indicator <= 0.5 && ...
    kantorovich_radius_bound <= single_radius && ...
    nnz(cbase(U+1:end) >= -opts.power_active_tol) == 0 && ...
    kantorovich_radius_bound <= entry_power_radius_slack && ...
    inactive_sinr_safe;
majorant_iterate = 0;
majorant_history = zeros(opts.sinr_retraction_max_iter+1,1);
if isfinite(kantorovich_indicator) && kantorovich_indicator <= 0.5 && ...
        sigma_min > 0
    scaled_lipschitz = vector_jacobian_lipschitz/sigma_min;
    for majorant_idx = 1:opts.sinr_retraction_max_iter
        majorant_value = 0.5*scaled_lipschitz*majorant_iterate^2 - ...
            majorant_iterate + step_newton_norm;
        majorant_derivative = ...
            scaled_lipschitz*majorant_iterate-1;
        majorant_iterate = majorant_iterate - ...
            majorant_value/majorant_derivative;
        majorant_history(majorant_idx+1) = majorant_iterate;
    end
    distance_after_finite_steps = max( ...
        kantorovich_radius_bound-majorant_iterate,0);
    active_residual_bound = (sigma_max + ...
        vector_jacobian_lipschitz*kantorovich_radius_bound) * ...
        distance_after_finite_steps;
else
    distance_after_finite_steps = inf;
    active_residual_bound = inf;
end
kantorovich_finite_step_certificate = kantorovich_tr_certificate && ...
    active_residual_bound <= ...
    opts.sinr_retraction_margin+opts.feas_tol;

[analytic_upper_bound,analytic_power_safe,retraction_error_bound] = ...
    one_step_upper_bound(zentry,centry,Gentry,step_code,cbase,idx, ...
    constraint_hessian_norm);
code_step_uncapped = step_code_norm <= single_radius*(1+1e-12);
analytic_one_step_certificate = isempty(idx) || ...
    (code_step_uncapped && analytic_power_safe && ...
    all(analytic_upper_bound <= opts.feas_tol));

% Apply the same capped Newton-like normal correction used by AS RTR, but
% without the objective/tCG step.  This isolates feasibility restoration.
zrestore = zentry;
pure_path_length = 0;
pure_step_count = 0;
pure_constraint_history = nan(opts.sinr_retraction_max_iter+1,1);
pure_constraint_history(1) = max([0;centry(1:U)]);
for correction_iter = 1:opts.sinr_retraction_max_iter
    [cq,Gq] = constraint_derivatives(zrestore);
    [Aq,idxq,Qpowerq] = correction_jacobian(cq,Gq,cbase);
    if isempty(idxq)
        pure_constraint_history(correction_iter+1:end) = ...
            pure_constraint_history(correction_iter);
        break;
    end
    residualq = cq(idxq)+opts.sinr_retraction_margin;
    normal_matrix = Aq'*Aq + ...
        opts.sinr_retraction_regularization*eye(numel(idxq));
    correction = -Aq*(normal_matrix\residualq);
    correction = correction-Qpowerq*(Qpowerq'*correction);
    if norm(correction) > single_radius
        correction = correction*(single_radius/norm(correction));
    end
    pure_path_length = pure_path_length+norm(correction);
    pure_step_count = pure_step_count+1;
    zrestore = retract_active_power(zrestore,correction,cbase);
    [cafter,~] = constraint_derivatives(zrestore);
    pure_constraint_history(correction_iter+1) = ...
        max([0;cafter(1:U)]);
end
[crestore,~,raw_restore] = constraint_derivatives(zrestore);

cert.entry_wcell = matrix_to_wcell(cfg,Wentry);
cert.entry_z = zentry;
cert.entry_constraints = centry;
cert.entry_raw = raw_entry;
cert.entry_max_sinr_violation = max([0;centry(1:U)]);
cert.entry_max_power_violation = max([0;centry(U+1:end)]);
cert.entry_min_sinr_dB = pow2db(max(raw_entry.min_sinr,realmin));
cert.entry_max_ap_power_ratio = raw_entry.max_ap_power/cfg.power.Pmax;
cert.active_sinr = idx(:);
cert.active_sinr_count = numel(idx);
cert.active_power_count = nnz(cbase(U+1:end) >= -opts.power_active_tol);
cert.jacobian_rank = jacobian_rank;
cert.sigma_min = sigma_min;
cert.sigma_max = sigma_max;
cert.jacobian_condition = jacobian_condition;
cert.step_linear = step_linear;
cert.step_linear_norm = step_linear_norm;
cert.linear_qp_success = linear_qp_success;
cert.linear_qp_kkt = linear_qp_kkt;
cert.linearized_qp_max_residual = max_or_zero(linearized_qp_residual);
cert.step_code = step_code;
cert.step_code_norm = step_code_norm;
cert.step_newton = step_newton;
cert.step_newton_norm = step_newton_norm;
cert.regularized_step_relative_error = regularized_step_relative_error;
cert.linearized_code_max_residual = ...
    max_or_zero(linearized_code_residual);
cert.single_correction_radius = single_radius;
cert.total_correction_path_budget = total_radius;
cert.fri_linear_single = safe_ratio(step_linear_norm,single_radius);
cert.fri_linear_total = safe_ratio(step_linear_norm,total_radius);
cert.fri_code_single = safe_ratio(step_code_norm,single_radius);
cert.fri_code_total = safe_ratio(step_code_norm,total_radius);
cert.constraint_hessian_norm = constraint_hessian_norm;
cert.vector_jacobian_lipschitz = vector_jacobian_lipschitz;
cert.kantorovich_indicator = kantorovich_indicator;
cert.kantorovich_radius_bound = kantorovich_radius_bound;
cert.entry_power_radius_slack = entry_power_radius_slack;
cert.full_row_rank = full_row_rank;
cert.inactive_sinr = inactive_sinr;
cert.inactive_sinr_upper_bound = inactive_sinr_upper_bound;
cert.inactive_sinr_safe = inactive_sinr_safe;
cert.kantorovich_tr_certificate = kantorovich_tr_certificate;
cert.majorant_history = majorant_history;
cert.distance_after_finite_steps = distance_after_finite_steps;
cert.active_residual_bound = active_residual_bound;
cert.kantorovich_finite_step_certificate = ...
    kantorovich_finite_step_certificate;
cert.code_step_uncapped = code_step_uncapped;
cert.analytic_power_safe = analytic_power_safe;
cert.retraction_error_bound = retraction_error_bound;
cert.analytic_upper_bound = analytic_upper_bound;
cert.analytic_max_upper_bound = max_or_negative_inf(analytic_upper_bound);
cert.analytic_one_step_certificate = analytic_one_step_certificate;
cert.pure_restoration_constraints = crestore;
cert.pure_restoration_raw = raw_restore;
cert.pure_restoration_max_sinr_violation = ...
    max([0;crestore(1:U)]);
cert.pure_restoration_feasible = ...
    max([0;crestore]) <= opts.feas_tol;
cert.pure_restoration_min_sinr_dB = ...
    pow2db(max(raw_restore.min_sinr,realmin));
cert.pure_restoration_path_length = pure_path_length;
cert.pure_restoration_step_count = pure_step_count;
cert.pure_restoration_constraint_history = pure_constraint_history;
cert.opts = opts;

    function [c,Gz,raw] = constraint_derivatives(zq)
        Wq = real_vector_to_complex_matrix(zq,N,U);
        [c,Gcomplex,raw] = compute_constraints_grad( ...
            cfg,channels,smat,Wq);
        Gz = zeros(dim,numel(c));
        for jj = 1:numel(c)
            Gz(:,jj) = complex_matrix_to_real_vector(Gcomplex(:,:,jj));
        end
    end

    function [Aq,idxq,Qpowerq] = correction_jacobian( ...
            cq,Gq,cbaseq)
        idxq = find(cq(1:U) > -opts.sinr_retraction_margin);
        active_power = cbaseq(U+1:end) >= -opts.power_active_tol;
        if any(active_power)
            Qpowerq = orth(Gq(:,U+find(active_power)),opts.rank_tol);
        else
            Qpowerq = zeros(dim,0);
        end
        Aq = Gq(:,idxq)-Qpowerq*(Qpowerq'*Gq(:,idxq));
    end

    function znew = retract_active_power(zq,step,cbaseq)
        Wq = real_vector_to_complex_matrix(zq+step,N,U);
        for retract_power_idx = 1:cfg.sys.M
            retract_rows = ...
                (retract_power_idx-1)*cfg.sys.Nt+(1:cfg.sys.Nt);
            block = Wq(retract_rows,:);
            power = real(sum(abs(block(:)).^2));
            is_active_power = cbaseq(U+retract_power_idx) >= ...
                -opts.power_active_tol;
            if is_active_power && power > 0
                block = block*sqrt(cfg.power.Pmax/power);
            elseif power > cfg.power.Pmax
                block = block*sqrt(cfg.power.Pmax/power);
            end
            Wq(retract_rows,:) = block;
        end
        znew = complex_matrix_to_real_vector(Wq);
    end

    function values = sinr_hessian_norms(user_indices)
        values = zeros(numel(user_indices),1);
        coefficient = max(1,cfg.comm.Gamma);
        for hessian_idx = 1:numel(user_indices)
            hu = channels.h_stack(:,user_indices(hessian_idx));
            values(hessian_idx) = 2*coefficient*real(hu'*hu) / ...
                max(cfg.comm.noise_power,realmin);
        end
    end

    function [upper_bound,power_safe,error_bound] = one_step_upper_bound( ...
            zq,cq,Gq,step,cbaseq,user_indices,hessian_norm)
        if isempty(user_indices)
            upper_bound = zeros(0,1);
            power_safe = true;
            error_bound = 0;
            return;
        end
        Wq = real_vector_to_complex_matrix(zq,N,U);
        Smat = real_vector_to_complex_matrix(step,N,U);
        active_power = cbaseq(U+1:end) >= -opts.power_active_tol;
        error_sq = 0;
        power_safe = true;
        for bound_power_idx = 1:cfg.sys.M
            bound_rows = ...
                (bound_power_idx-1)*cfg.sys.Nt+(1:cfg.sys.Nt);
            xnorm = norm(Wq(bound_rows,:),'fro');
            snorm = norm(Smat(bound_rows,:),'fro');
            if active_power(bound_power_idx)
                error_m = sqrt(cfg.power.Pmax+snorm^2) - ...
                    sqrt(cfg.power.Pmax);
                error_sq = error_sq+error_m^2;
            elseif xnorm+snorm > sqrt(cfg.power.Pmax)*(1+1e-12)
                power_safe = false;
            end
        end
        error_bound = sqrt(error_sq);
        linear_part = cq(user_indices)+Gq(:,user_indices)'*step;
        step_norm = norm(step);
        upper_bound = zeros(numel(user_indices),1);
        for bound_constraint_idx = 1:numel(user_indices)
            Lk = hessian_norm(bound_constraint_idx);
            gradient_norm = norm( ...
                Gq(:,user_indices(bound_constraint_idx)));
            upper_bound(bound_constraint_idx) = ...
                linear_part(bound_constraint_idx) + ...
                0.5*Lk*step_norm^2 + ...
                (gradient_norm+Lk*step_norm)*error_bound + ...
                0.5*Lk*error_bound^2;
        end
    end
end

function [step,success,kkt_residual] = minimum_norm_inequality_step( ...
    A,residual,rank_tol)
% Solve min 0.5||s||^2 subject to A'*s+residual <= 0 through its dual.
n = size(A,1);
q = size(A,2);
if q == 0
    step = zeros(n,1);
    success = true;
    kkt_residual = 0;
    return;
end
Q = A'*A;
Q = 0.5*(Q+Q');

if q <= 12
    best_norm = inf;
    best_step = zeros(n,1);
    best_lambda = zeros(q,1);
    for mask = 1:(2^q-1)
        active = logical(bitget(mask,1:q));
        Qa = Q(active,active);
        ra = residual(active);
        lambda_active = pinv(Qa,rank_tol)*ra;
        if any(lambda_active < -1e-9)
            continue;
        end
        lambda = zeros(q,1);
        lambda(active) = max(lambda_active,0);
        candidate = -A*lambda;
        violation = A'*candidate+residual;
        if max(violation) <= 1e-7*max(1,norm(residual,inf)) && ...
                norm(candidate) < best_norm
            best_norm = norm(candidate);
            best_step = candidate;
            best_lambda = lambda;
        end
    end
    step = best_step;
    success = isfinite(best_norm);
    if ~success
        step = -A*(pinv(Q,rank_tol)*residual);
        best_lambda = pinv(Q,rank_tol)*residual;
    end
else
    lambda = zeros(q,1);
    lipschitz = max(eig(Q));
    if lipschitz <= realmin
        step = zeros(n,1);
        success = false;
        kkt_residual = inf;
        return;
    end
    for iter = 1:20000
        next_lambda = max(0,lambda-(Q*lambda-residual)/lipschitz);
        if norm(next_lambda-lambda) <= 1e-11*max(1,norm(lambda))
            lambda = next_lambda;
            break;
        end
        lambda = next_lambda;
    end
    best_lambda = lambda;
    step = -A*lambda;
    success = max(A'*step+residual) <= ...
        1e-6*max(1,norm(residual,inf));
end

primal_residual = max([0;A'*step+residual]);
dual_residual = norm(min(best_lambda,Q*best_lambda-residual),inf);
complementarity = norm(best_lambda.*(A'*step+residual),inf);
kkt_residual = max([primal_residual;dual_residual;complementarity]);
end

function opts = fill_defaults(opts)
if ~isfield(opts,'delta0'); opts.delta0=5e-2; end
if ~isfield(opts,'feas_tol'); opts.feas_tol=1e-6; end
if ~isfield(opts,'sinr_retraction_margin')
    opts.sinr_retraction_margin=1e-5;
end
if ~isfield(opts,'sinr_retraction_max_iter')
    opts.sinr_retraction_max_iter=4;
end
if ~isfield(opts,'sinr_retraction_regularization')
    opts.sinr_retraction_regularization=1e-12;
end
if ~isfield(opts,'sinr_retraction_max_fraction')
    opts.sinr_retraction_max_fraction=1;
end
if ~isfield(opts,'power_active_tol'); opts.power_active_tol=1e-5; end
if ~isfield(opts,'rank_tol'); opts.rank_tol=1e-10; end
end

function value = max_or_zero(x)
if isempty(x); value = 0; else; value = max(x); end
end

function value = max_or_negative_inf(x)
if isempty(x); value = -inf; else; value = max(x); end
end

function value = safe_ratio(numerator,denominator)
if denominator > 0
    value = numerator/denominator;
elseif numerator == 0
    value = 0;
else
    value = inf;
end
end
