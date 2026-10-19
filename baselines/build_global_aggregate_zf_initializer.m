function [wcell,info] = build_global_aggregate_zf_initializer( ...
    cfg,channels,smat,opts)
%BUILD_GLOBAL_AGGREGATE_ZF_INITIALIZER Centralized no-CVX ZF warm start.
%
% Stack the user channels into H in C^(U x MNt). When H has full row rank,
% D = H^dagger satisfies H*D = I. The strict-SINR floor is
%
%             W_floor = D*diag(sqrt((1+delta)*Gamma*sigma^2)).
%
% A complex-matrix Dykstra projection redistributes this point in the null
% space of H while preserving H*W = diag(b) and enforcing every per-AP
% Frobenius-norm ball. Optional projected CRLB steps use the same convex
% projection, so every accepted initializer iterate remains exact ZF.

if nargin < 3
    smat = [];
end
if nargin < 4
    opts = struct();
end
opts = fill_defaults(opts);
if ~(isfinite(opts.sinr_margin_fraction) && ...
        opts.sinr_margin_fraction > 0)
    error('sinr_margin_fraction must be finite and strictly positive.');
end
if ~(isfinite(opts.power_margin_fraction) && ...
        opts.power_margin_fraction > 0 && ...
        opts.power_margin_fraction < 1)
    error('power_margin_fraction must lie strictly between zero and one.');
end

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;
Pmax = cfg.power.Pmax;
timer_total = tic;

H = complex(zeros(U,N));
for u = 1:U
    for m = 1:M
        H(u,block_indices(m,Nt)) = channels.h{m,u}';
    end
end

timer_direction = tic;
[left_singular,sigma_matrix,right_singular] = svd(H,'econ');
singular_values = diag(sigma_matrix);
if isempty(singular_values)
    error('The aggregate channel matrix is empty.');
end
rank_threshold = opts.rank_tolerance*max(singular_values);
aggregate_rank = nnz(singular_values > rank_threshold);
if N < U || aggregate_rank < U
    error(['Global aggregate ZF requires MNt >= U and rank(H)=U; ' ...
        'received MNt=%d, U=%d, rank(H)=%d.'],N,U,aggregate_rank);
end
D = right_singular(:,1:U)*diag(1./singular_values(1:U))* ...
    left_singular(:,1:U)';
direction_runtime = toc(timer_direction);

gamma = cfg.comm.Gamma;
if isscalar(gamma)
    gamma = repmat(gamma,U,1);
else
    gamma = gamma(:);
end
if numel(gamma) ~= U
    error('cfg.comm.Gamma must be scalar or have one entry per user.');
end
strict_gamma = gamma*(1+opts.sinr_margin_fraction);
required_amplitude = sqrt(strict_gamma*cfg.comm.noise_power);
B = diag(required_amplitude);
Wfloor = D*B;

power_cap = (1-opts.power_margin_fraction)*Pmax;
floor_ap_power = ap_power(Wfloor);
closed_form_power_feasible = all(floor_ap_power <= ...
    power_cap*(1+opts.projection_tolerance));

timer_projection = tic;
if closed_form_power_feasible
    W = Wfloor;
    projection_info = projection_diagnostics(W,0,true);
else
    [W,projection_info] = project_global_zf_power_dykstra( ...
        Wfloor,H,D,B,M,Nt,sqrt(power_cap), ...
        opts.max_projection_cycles,opts.projection_tolerance);
end
projection_runtime = toc(timer_projection);

[current_obj,current_crlb,current_grad] = objective_and_gradient(W);
objective_history = current_obj;
crlb_history = current_crlb;
accepted_history = false(0,1);
projection_cycle_history = projection_info.cycles;
radius = opts.initial_radius*sqrt(M*Pmax);

timer_refine = tic;
if projection_info.converged && opts.max_refine_iter > 0 && ~isempty(smat)
    for refine_iter = 1:opts.max_refine_iter
        grad_norm = norm(current_grad,'fro');
        if ~isfinite(grad_norm) || grad_norm <= opts.gradient_tolerance
            break;
        end
        normalized_grad = current_grad/grad_norm;
        local_radius = radius;
        accepted = false;
        last_projection_cycles = 0;
        for backtrack_iter = 1:opts.max_backtracking
            raw_trial = W-local_radius*normalized_grad;
            [trial,trial_projection] = project_global_zf_power_dykstra( ...
                raw_trial,H,D,B,M,Nt,sqrt(power_cap), ...
                opts.max_projection_cycles,opts.projection_tolerance);
            last_projection_cycles = trial_projection.cycles;
            if ~trial_projection.converged
                local_radius = local_radius*opts.radius_shrink;
                continue;
            end
            [trial_obj,trial_crlb,trial_grad] = ...
                objective_and_gradient(trial);
            sufficient_decrease = trial_obj <= current_obj - ...
                opts.objective_tolerance*max(1,abs(current_obj));
            if isfinite(trial_obj) && sufficient_decrease
                W = trial;
                current_obj = trial_obj;
                current_crlb = trial_crlb;
                current_grad = trial_grad;
                radius = min(opts.radius_expand*local_radius, ...
                    opts.max_radius*sqrt(M*Pmax));
                accepted = true;
                break;
            end
            local_radius = local_radius*opts.radius_shrink;
        end
        accepted_history(end+1,1) = accepted; %#ok<AGROW>
        objective_history(end+1,1) = current_obj; %#ok<AGROW>
        crlb_history(end+1,1) = current_crlb; %#ok<AGROW>
        projection_cycle_history(end+1,1) = ...
            last_projection_cycles; %#ok<AGROW>
        if ~accepted || local_radius <= opts.min_radius*sqrt(M*Pmax)
            break;
        end
    end
end
refine_runtime = toc(timer_refine);

wcell = matrix_to_wcell(cfg,W);
metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
response = H*W;
off_diagonal_response = response-diag(diag(response));
zf_residual_absolute = norm(off_diagonal_response,'fro');
zf_residual = zf_residual_absolute/ ...
    max(norm(diag(diag(response)),'fro'),realmin);
response_residual_absolute = norm(response-B,'fro');
response_residual = response_residual_absolute/max(norm(B,'fro'),realmin);
entry_feasible = projection_info.converged && all(metrics.sinr >= ...
    gamma*(1-opts.feasibility_tolerance)) && ...
    all(metrics.ap_power <= Pmax*(1+opts.feasibility_tolerance));

normal_certificate = struct();
if opts.evaluate_normal_certificate
    certificate_options.delta0 = opts.certificate_delta0;
    certificate_options.feas_tol = opts.certificate_feas_tol;
    certificate_options.sinr_retraction_margin = ...
        opts.certificate_sinr_retraction_margin;
    certificate_options.sinr_retraction_max_iter = ...
        opts.certificate_max_normal_steps;
    normal_certificate = evaluate_as_rtr_initialization_certificate( ...
        cfg,channels,smat,wcell,certificate_options);
end

info.method = 'Global aggregate theorem-aware ZF initializer';
info.wcell = wcell;
info.aggregate_channel = H;
info.aggregate_pseudoinverse = D;
info.aggregate_rank = aggregate_rank;
info.aggregate_dimension = [U,N];
info.singular_values = singular_values;
info.condition_number = singular_values(1)/singular_values(U);
info.required_amplitude = required_amplitude;
info.strict_gamma = strict_gamma;
info.response_target = B;
info.response = response;
info.zf_residual = zf_residual;
info.zf_residual_absolute = zf_residual_absolute;
info.response_residual = response_residual;
info.response_residual_absolute = response_residual_absolute;
info.power_floor_matrix = Wfloor;
info.floor_ap_power = floor_ap_power;
info.closed_form_power_feasible = closed_form_power_feasible;
info.power_cap = power_cap;
info.projection_info = projection_info;
info.entry_feasible = entry_feasible;
info.feasibility_entry_certificate = entry_feasible;
info.metrics = metrics;
info.objective_history = objective_history;
info.crlb_history = crlb_history;
info.accepted_history = accepted_history;
info.num_accepted_refine_steps = nnz(accepted_history);
info.projection_cycle_history = projection_cycle_history;
info.normal_certificate = normal_certificate;
info.direction_runtime = direction_runtime;
info.projection_runtime = projection_runtime;
info.refine_runtime = refine_runtime;
info.total_runtime = toc(timer_total);
info.options = opts;

    function values = ap_power(Wvalue)
        values = zeros(M,1);
        for m1 = 1:M
            values(m1) = norm(Wvalue(block_indices(m1,Nt),:),'fro')^2;
        end
    end

    function [obj,crlb,gradW] = objective_and_gradient(Wvalue)
        if isempty(smat)
            obj = NaN;
            crlb = NaN;
            gradW = complex(zeros(size(Wvalue)));
            return;
        end
        cfg_objective = cfg;
        cfg_objective.alm.objective_type = 'log_crlb_ratio';
        [obj,gradW,gout] = compute_crlb_objective_grad( ...
            cfg_objective,smat,Wvalue);
        crlb = gout.crlb;
    end

    function diagnostics = projection_diagnostics(Wvalue,cycles,converged)
        local_response_residual = norm(H*Wvalue-B,'fro')/ ...
            max(1,norm(B,'fro'));
        local_power = ap_power(Wvalue);
        local_power_violation = max([0;local_power-power_cap])/ ...
            max(1,power_cap);
        diagnostics.cycles = cycles;
        diagnostics.converged = converged;
        diagnostics.response_residual = local_response_residual;
        diagnostics.max_power_violation = local_power_violation;
        diagnostics.ap_power = local_power;
    end
end

function [W,info] = project_global_zf_power_dykstra( ...
    V,H,D,B,M,Nt,power_radius,max_cycles,tolerance)
% Dykstra projection onto H*W=B and the M disjoint AP power balls.

[N,U] = size(V);
if N ~= M*Nt || ~isequal(size(B),[U,U])
    error('Inconsistent aggregate-ZF projection geometry.');
end

W = V;
affine_correction = complex(zeros(N,U));
ball_correction = complex(zeros(Nt,U,M));
response_scale = max(1,norm(B,'fro'));
power_cap = power_radius^2;
converged = false;
response_residual = inf;
max_power_violation = inf;

for cycle = 1:max_cycles
    Y = W+affine_correction;
    Wnew = Y+D*(B-H*Y);
    affine_correction = Y-Wnew;
    W = Wnew;

    for m = 1:M
        rows = block_indices(m,Nt);
        Yblock = W(rows,:)+ball_correction(:,:,m);
        block_norm = norm(Yblock,'fro');
        if block_norm > power_radius
            Wblock = Yblock*(power_radius/block_norm);
        else
            Wblock = Yblock;
        end
        ball_correction(:,:,m) = Yblock-Wblock;
        W(rows,:) = Wblock;
    end

    response_residual = norm(H*W-B,'fro')/response_scale;
    ap_power = zeros(M,1);
    for m = 1:M
        ap_power(m) = norm(W(block_indices(m,Nt),:),'fro')^2;
    end
    max_power_violation = max([0;ap_power-power_cap])/ ...
        max(1,power_cap);
    if response_residual <= tolerance && ...
            max_power_violation <= tolerance
        converged = true;
        break;
    end
end

info.cycles = cycle;
info.converged = converged;
info.response_residual = response_residual;
info.max_power_violation = max_power_violation;
info.ap_power = ap_power;
end

function opts = fill_defaults(opts)
if ~isfield(opts,'sinr_margin_fraction')
    opts.sinr_margin_fraction = 1e-3;
end
if ~isfield(opts,'power_margin_fraction')
    opts.power_margin_fraction = 1e-3;
end
if ~isfield(opts,'rank_tolerance'); opts.rank_tolerance=1e-11; end
if ~isfield(opts,'max_projection_cycles')
    opts.max_projection_cycles=300;
end
if ~isfield(opts,'projection_tolerance')
    opts.projection_tolerance=1e-9;
end
% The paper-facing entry is feasibility-only; CRLB descent belongs to AS-RTR.
if ~isfield(opts,'max_refine_iter'); opts.max_refine_iter=0; end
if ~isfield(opts,'max_backtracking'); opts.max_backtracking=8; end
if ~isfield(opts,'initial_radius'); opts.initial_radius=0.05; end
if ~isfield(opts,'max_radius'); opts.max_radius=0.5; end
if ~isfield(opts,'min_radius'); opts.min_radius=1e-7; end
if ~isfield(opts,'radius_shrink'); opts.radius_shrink=0.5; end
if ~isfield(opts,'radius_expand'); opts.radius_expand=1.5; end
if ~isfield(opts,'gradient_tolerance'); opts.gradient_tolerance=1e-8; end
if ~isfield(opts,'objective_tolerance')
    opts.objective_tolerance=1e-9;
end
if ~isfield(opts,'feasibility_tolerance')
    opts.feasibility_tolerance=1e-4;
end
if ~isfield(opts,'evaluate_normal_certificate')
    opts.evaluate_normal_certificate=false;
end
if ~isfield(opts,'certificate_delta0'); opts.certificate_delta0=5e-2; end
if ~isfield(opts,'certificate_feas_tol')
    opts.certificate_feas_tol=1e-6;
end
if ~isfield(opts,'certificate_sinr_retraction_margin')
    opts.certificate_sinr_retraction_margin=1e-5;
end
if ~isfield(opts,'certificate_max_normal_steps')
    opts.certificate_max_normal_steps=4;
end
end
