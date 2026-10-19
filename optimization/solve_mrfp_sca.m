function [sol,hist] = solve_mrfp_sca( ...
    cfg,channels,smat,init_wcell,opts)
%SOLVE_MRFP_SCA Multi-ratio / matrix-fractional SCA baseline.
%
% The physical CRLB, apart from its constant prefactor, is
%
%   tr(G(W)^(-1))
%     = 1/(g_a-g_ab^2/g_b) + 1/(g_b-g_ab^2/g_a),
%
% where G=[g_a,g_ab;g_ab,g_b] is the 2-by-2 Fisher information matrix.
% This is a two-ratio minimization problem.  Eliminating the inverse
% quadratic-transform auxiliaries gives the two INV_POS terms used below.
%
% Direct application of a scalar quadratic transform is not convex in the
% rank-one beamforming vectors because G(W) is quadratic.  At Wk, this
% solver therefore uses the global Loewner-affine minorant
%
%   G(W) >= Gtilde_k(W)
%
% obtained from |s|^2 >= 2*real(conj(sk)*s)-|sk|^2 for each sensing scalar
% s=a_t^T*w.  Since inversion reverses the positive-semidefinite order,
%
%   tr(G(W)^(-1)) <= tr(Gtilde_k(W)^(-1)).
%
% The right-hand side is a first-order tight majorizer.  Its minimization,
% together with phase-normalized SINR constraints and per-AP power balls,
% is a vector-domain SOCP.  No lifted N-by-N covariance variable is used.
%
% Required entry condition:
%   INIT_WCELL must satisfy the final SINR and AP-power constraints.  Using
%   the same feasible initializer as COSTAR makes the comparison traceable.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(cfg,opts);

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;
Pmax = cfg.power.Pmax;
Gamma = cfg.comm.Gamma;
noise_power = cfg.comm.noise_power;
noise_standard_deviation = max(sqrt(noise_power),realmin);
beta2 = cfg.radio.effective_bandwidth_sq(:);

if isempty(init_wcell)
    error(['MRFP-SCA requires a feasible rank-one initializer. Pass the ' ...
        'same low-complexity/ZF point used by COSTAR.']);
end
if ~isfield(channels,'h_stack') || size(channels.h_stack,1) ~= N
    error('channels.h_stack is missing or inconsistent with cfg.sys.');
end
required_sensing_fields = {'at','ca','cb','cab'};
for field_index = 1:numel(required_sensing_fields)
    if ~isfield(smat,required_sensing_fields{field_index})
        error('smat.%s is required by MRFP-SCA.', ...
            required_sensing_fields{field_index});
    end
end
for m = 1:M
    geometry_matrix = hermitianize([ ...
        smat.ca(m),smat.cab(m);smat.cab(m),smat.cb(m)]);
    geometry_scale = max(norm(geometry_matrix,2),realmin);
    if min(real(eig(geometry_matrix))) < -1e-10*geometry_scale
        error(['The sensing geometry matrix C_%d is not positive ' ...
            'semidefinite, so the MRFP Loewner majorizer is invalid.'],m);
    end
end

Wk = vectorize_wcell_to_matrix(cfg,init_wcell);
Wk = phase_align_desired_signals(Wk,channels.h_stack);
initial_wcell = matrix_to_wcell(cfg,Wk);
metrics_k = evaluate_beamforming_solution( ...
    cfg,channels,smat,initial_wcell);
if ~metrics_are_feasible(metrics_k,cfg,opts.feasibility_tol)
    error(['MRFP-SCA entry point is infeasible: min SINR %.3f dB ' ...
        '(target %.3f dB), max AP power %.6g (limit %.6g).'], ...
        pow2db(max(min(metrics_k.sinr),realmin)), ...
        cfg.comm.Gamma_dB,max(metrics_k.ap_power),Pmax);
end
if ~(isfinite(metrics_k.crlb) && metrics_k.crlb > 0 && ...
        isfinite(metrics_k.det) && metrics_k.det > 0)
    error('MRFP-SCA entry point does not have a positive-definite Fisher matrix.');
end

power_scale = max(M*Pmax,realmin);

max_points = opts.max_iter+1;
crlb_history = nan(max_points,1);
ratio_history = nan(max_points,1);
min_sinr_history = nan(max_points,1);
max_power_history = nan(max_points,1);
rel_decrease_history = nan(max_points,1);
step_norm_history = nan(max_points,1);
surrogate_history = nan(max_points,1);
solve_time_history = nan(max_points,1);
cumulative_time_history = nan(max_points,1);
backtracking_history = nan(max_points,1);
accepted_history = false(max_points,1);
cvx_status_history = strings(max_points,1);
z_snapshot_history = nan(2*N*U,max_points);

crlb_history(1) = metrics_k.crlb;
ratio_history(1) = metrics_k.crlb_ratio;
min_sinr_history(1) = min(metrics_k.sinr);
max_power_history(1) = max(metrics_k.ap_power);
rel_decrease_history(1) = 0;
step_norm_history(1) = 0;
surrogate_history(1) = metrics_k.crlb_ratio;
solve_time_history(1) = 0;
cumulative_time_history(1) = 0;
backtracking_history(1) = 0;
accepted_history(1) = true;
cvx_status_history(1) = "initial";
z_snapshot_history(:,1) = complex_matrix_to_real_vector(Wk);

if opts.verbose
    fprintf(['[MRFP-SCA] entry CRLB=%.6g, minSINR=%.3f dB, ' ...
        'maxP=%.4g W\n'],metrics_k.crlb, ...
        pow2db(max(min(metrics_k.sinr),realmin)), ...
        max(metrics_k.ap_power));
end

timer_total = tic;
num_completed = 0;
last_status = "not_started";
message = 'maximum iteration budget reached';

for iteration = 1:opts.max_iter
    Wref = Wk;
    metrics_ref = metrics_k;

    sensing_reference = zeros(M,U);
    for m = 1:M
        index_m = block_indices(m,Nt);
        steering_m = smat.at{m};
        for u = 1:U
            sensing_reference(m,u) = ...
                transpose(steering_m)*Wref(index_m,u);
        end
    end

    fim_scale = max(0.5*(metrics_ref.ga+metrics_ref.gb),realmin);
    solve_timer = tic;
    try
        cvx_begin
            if opts.cvx_quiet
                cvx_quiet true
            end
            if ~isempty(opts.cvx_solver)
                cvx_solver(opts.cvx_solver)
            end

            variable W(N,U) complex
            variable kappa_ab nonnegative
            variable tau_a nonnegative
            variable tau_b nonnegative
            expression g_a_minor
            expression g_b_minor
            expression g_ab_minor
            expression received(U,U)
            expression proximal_term

            g_a_minor = 0;
            g_b_minor = 0;
            g_ab_minor = 0;
            for m = 1:M
                index_m = block_indices(m,Nt);
                steering_m = smat.at{m};
                for u = 1:U
                    sensing_value = transpose(steering_m)*W(index_m,u);
                    affine_power = ...
                        2*real(conj(sensing_reference(m,u))*sensing_value) ...
                        - abs(sensing_reference(m,u))^2;
                    normalized_weight = beta2(u)/fim_scale;
                    g_a_minor = g_a_minor + normalized_weight* ...
                        smat.ca(m)*affine_power;
                    g_b_minor = g_b_minor + normalized_weight* ...
                        smat.cb(m)*affine_power;
                    g_ab_minor = g_ab_minor + normalized_weight* ...
                        smat.cab(m)*affine_power;
                end
            end

            % Normalize received amplitudes by the noise standard
            % deviation. This leaves every SINR cone unchanged and avoids
            % SDPT3 failures caused by 1e-5-scale channel/noise entries.
            received = (channels.h_stack'*W)/noise_standard_deviation;
            proximal_term = opts.prox_weight* ...
                sum_square_abs(W(:)-Wref(:))/power_scale;

            minimize(inv_pos(tau_a)+inv_pos(tau_b)+proximal_term)
            subject to
                g_a_minor >= opts.fim_floor;
                g_b_minor >= opts.fim_floor;
                kappa_ab >= g_ab_minor;
                kappa_ab >= -g_ab_minor;
                tau_a >= opts.fim_floor;
                tau_b >= opts.fim_floor;

                % tau_a <= g_a-g_ab^2/g_b and
                % tau_b <= g_b-g_ab^2/g_a, in rotated-SOC form.
                norm([2*kappa_ab; ...
                    g_b_minor-g_a_minor+tau_a]) <= ...
                    g_b_minor+g_a_minor-tau_a;
                norm([2*kappa_ab; ...
                    g_a_minor-g_b_minor+tau_b]) <= ...
                    g_a_minor+g_b_minor-tau_b;

                % Phase-normalized downlink SINR SOC constraints.
                for u = 1:U
                    interference_users = [1:u-1,u+1:U];
                    imag(received(u,u)) == 0;
                    real(received(u,u)) >= 0;
                    norm([transpose(received(u,interference_users)); ...
                        1]) <= ...
                        real(received(u,u))/sqrt(Gamma);
                end

                % Per-AP power balls.
                for m = 1:M
                    index_m = block_indices(m,Nt);
                    norm(reshape(W(index_m,:),Nt*U,1)) <= sqrt(Pmax);
                end
        cvx_end
    catch solver_exception
        solve_elapsed = toc(solve_timer);
        last_status = "exception";
        message = ['CVX exception: ' solver_exception.message];
        if opts.verbose
            fprintf('[MRFP-SCA] iteration %d failed: %s\n', ...
                iteration,solver_exception.message);
        end
        solve_time_history(iteration+1) = solve_elapsed;
        break;
    end
    solve_elapsed = toc(solve_timer);
    last_status = string(cvx_status);

    solved = contains(lower(last_status),"solved") || ...
        contains(lower(last_status),"inaccurate");
    if ~solved || any(~isfinite(W(:)))
        message = sprintf('CVX subproblem failed with status %s.',last_status);
        if opts.verbose
            fprintf('[MRFP-SCA] iteration %d status=%s; stopping.\n', ...
                iteration,last_status);
        end
        solve_time_history(iteration+1) = solve_elapsed;
        cvx_status_history(iteration+1) = last_status;
        break;
    end

    Wcandidate = full(W);
    accepted = false;
    selected_backtracking = 0;
    candidate_metrics = metrics_ref;
    Waccepted = Wref;
    for backtracking_index = 0:opts.max_backtracking
        step_fraction = 0.5^backtracking_index;
        Wtrial = Wref+step_fraction*(Wcandidate-Wref);
        trial_wcell = matrix_to_wcell(cfg,Wtrial);
        trial_metrics = evaluate_beamforming_solution( ...
            cfg,channels,smat,trial_wcell);

        feasible_trial = metrics_are_feasible( ...
            trial_metrics,cfg,opts.feasibility_tol);
        monotone_trial = isfinite(trial_metrics.crlb) && ...
            trial_metrics.crlb <= metrics_ref.crlb* ...
            (1+opts.accept_rel_tol);
        if feasible_trial && monotone_trial
            accepted = true;
            selected_backtracking = backtracking_index;
            Waccepted = Wtrial;
            candidate_metrics = trial_metrics;
            break;
        end
    end

    if ~accepted
        message = ['convex trial could not pass the physical feasibility ' ...
            'and monotone-CRLB safeguard'];
        if opts.verbose
            fprintf('[MRFP-SCA] iteration %d rejected after %d halvings.\n', ...
                iteration,opts.max_backtracking);
        end
        solve_time_history(iteration+1) = solve_elapsed;
        cvx_status_history(iteration+1) = last_status;
        break;
    end

    Wk = Waccepted;
    metrics_k = candidate_metrics;
    relative_decrease = max( ...
        (metrics_ref.crlb-metrics_k.crlb)/max(metrics_ref.crlb,realmin),0);
    relative_step = norm(Wk-Wref,'fro')/max(norm(Wref,'fro'),realmin);

    num_completed = iteration;
    history_index = iteration+1;
    crlb_history(history_index) = metrics_k.crlb;
    ratio_history(history_index) = metrics_k.crlb_ratio;
    min_sinr_history(history_index) = min(metrics_k.sinr);
    max_power_history(history_index) = max(metrics_k.ap_power);
    rel_decrease_history(history_index) = relative_decrease;
    step_norm_history(history_index) = relative_step;
    surrogate_history(history_index) = ...
        max(cvx_optval-opts.prox_weight* ...
        norm(Wcandidate-Wref,'fro')^2/power_scale,0)/fim_scale;
    solve_time_history(history_index) = solve_elapsed;
    cumulative_time_history(history_index) = toc(timer_total);
    backtracking_history(history_index) = selected_backtracking;
    accepted_history(history_index) = true;
    cvx_status_history(history_index) = last_status;
    z_snapshot_history(:,history_index) = complex_matrix_to_real_vector(Wk);

    if opts.verbose && (iteration == 1 || ...
            mod(iteration,opts.print_every) == 0)
        fprintf(['[MRFP-SCA] iter %02d: CRLB=%.6g, rel=%.3e, ' ...
            'minSINR=%.3f dB, step=%.3e, CVX=%s, time=%.3fs\n'], ...
            iteration,metrics_k.crlb,relative_decrease, ...
            pow2db(max(min(metrics_k.sinr),realmin)), ...
            relative_step,last_status,solve_elapsed);
    end

    if iteration >= opts.min_iter && ...
            relative_decrease <= opts.crlb_rel_tol
        message = 'relative CRLB tolerance reached';
        break;
    end
end

total_runtime = toc(timer_total);
used_points = 1:num_completed+1;
final_wcell = matrix_to_wcell(cfg,Wk);
final_feasible = metrics_are_feasible( ...
    metrics_k,cfg,opts.feasibility_tol);

hist.iter = (0:num_completed).';
hist.crlb = crlb_history(used_points);
hist.crlb_ratio = ratio_history(used_points);
hist.min_sinr = min_sinr_history(used_points);
hist.max_ap_power = max_power_history(used_points);
hist.relative_crlb_decrease = rel_decrease_history(used_points);
hist.relative_step_norm = step_norm_history(used_points);
hist.surrogate_ratio_upper_model = surrogate_history(used_points);
hist.solve_time_s = solve_time_history(used_points);
hist.cumulative_time_s = cumulative_time_history(used_points);
hist.backtracking_halvings = backtracking_history(used_points);
hist.accepted = accepted_history(used_points);
hist.cvx_status = cvx_status_history(used_points);
hist.snapshot_iter = hist.iter;
hist.z_snapshots = z_snapshot_history(:,used_points);
hist.method = 'MRFP-SCA';
hist.total_runtime = total_runtime;
hist.num_iterations = num_completed;
hist.monotone_crlb = all(diff(hist.crlb) <= ...
    opts.accept_rel_tol*max(hist.crlb(1:end-1),realmin));

sol.wcell = final_wcell;
sol.Wmat = Wk;
sol.z = complex_matrix_to_real_vector(Wk);
sol.metrics = metrics_k;
sol.success = final_feasible && isfinite(metrics_k.crlb);
sol.message = message;
sol.method = 'MRFP-SCA';
sol.full_name = ...
    'Multi-Ratio Fractional Programming with Fisher-Matrix SCA';
sol.initializer_method = 'user-supplied common feasible initializer';
sol.initial_metrics = evaluate_beamforming_solution( ...
    cfg,channels,smat,initial_wcell);
sol.num_iterations = num_completed;
sol.total_runtime = total_runtime;
sol.cvx_status = char(last_status);
sol.first_order_only = true;
sol.negative_curvature_escape = false;
sol.hist = hist;

if opts.verbose
    fprintf(['[MRFP-SCA] finished: %s; iterations=%d, CRLB=%.6g, ' ...
        'feasible=%d, total=%.3fs\n'],message,num_completed, ...
        metrics_k.crlb,final_feasible,total_runtime);
end
end

function W = phase_align_desired_signals(W,h_stack)
U = size(W,2);
for u = 1:U
    desired = h_stack(:,u)'*W(:,u);
    if abs(desired) > realmin
        W(:,u) = W(:,u)*exp(-1j*angle(desired));
    end
end
end

function feasible = metrics_are_feasible(metrics,cfg,tolerance)
feasible = all(isfinite(metrics.sinr)) && ...
    all(isfinite(metrics.ap_power)) && ...
    all(metrics.sinr >= cfg.comm.Gamma*(1-tolerance)) && ...
    all(metrics.ap_power <= cfg.power.Pmax*(1+tolerance));
end

function opts = fill_defaults(cfg,opts)
if isfield(cfg,'mrfp')
    defaults = cfg.mrfp;
else
    defaults = struct();
end
fallbacks.max_iter = 30;
fallbacks.min_iter = 2;
fallbacks.crlb_rel_tol = 1e-5;
fallbacks.fim_floor = 1e-8;
fallbacks.prox_weight = 1e-4;
fallbacks.accept_rel_tol = 1e-8;
fallbacks.max_backtracking = 12;
fallbacks.feasibility_tol = 2e-4;
fallbacks.cvx_quiet = true;
fallbacks.cvx_solver = '';
fallbacks.verbose = true;
fallbacks.print_every = 1;

names = fieldnames(fallbacks);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(defaults,name)
        defaults.(name) = fallbacks.(name);
    end
    if ~isfield(opts,name)
        opts.(name) = defaults.(name);
    end
end

opts.max_iter = max(1,round(opts.max_iter));
opts.min_iter = max(1,min(round(opts.min_iter),opts.max_iter));
opts.max_backtracking = max(0,round(opts.max_backtracking));
opts.print_every = max(1,round(opts.print_every));
if ~(isfinite(opts.fim_floor) && opts.fim_floor > 0)
    error('opts.fim_floor must be positive and finite.');
end
if ~(isfinite(opts.prox_weight) && opts.prox_weight >= 0)
    error('opts.prox_weight must be nonnegative and finite.');
end
if ~(isfinite(opts.crlb_rel_tol) && opts.crlb_rel_tol >= 0)
    error('opts.crlb_rel_tol must be nonnegative and finite.');
end
end
