function [sol,hist] = solve_sca_sdr_rank_continuation( ...
    cfg,channels,smat,init_wcell,opts)
%SOLVE_SCA_SDR_RANK_CONTINUATION Blockwise SDR-SCA with rank continuation.
%
% The lifted model is unchanged, but the block-diagonal covariance is
% represented by U independent N-by-N PSD cones. Rank one is imposed by a
% sequence of convex inner approximations
%
%   sum_u tr((I-v_{u,k}v_{u,k}^H)W_u) <= epsilon_k,
%
% where v_{u,k} is a dominant eigenvector of the previous covariance and
% epsilon_k is reduced gradually. Each subproblem minimizes the CRLB
% surrogate itself instead of a fixed CRLB-plus-rank weighted sum.

if nargin < 5
    opts = struct();
end
opts = fill_options(cfg,opts);

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;
Pmax = cfg.power.Pmax;
beta2 = cfg.radio.effective_bandwidth_sq(:);
power_scale = max(M*Pmax,realmin);
fisher_scale = sensing_upper_scale(cfg,smat,beta2);

if nargin < 4 || isempty(init_wcell)
    error('Rank-continuation SDR-SCA requires a feasible rank-one entry.');
end

Wmat0 = vectorize_wcell_to_matrix(cfg,init_wcell);
Winitial = cell(U,1);
for u = 1:U
    Winitial{u} = hermitianize(Wmat0(:,u)*Wmat0(:,u)');
end
entry_lifted = evaluate_lifted_solution(cfg,channels,smat,Winitial);
entry_recovered = evaluate_recovered(cfg,channels,smat,Winitial);
if ~lifted_feasible(entry_lifted,cfg,opts.feasibility_tol) || ...
        ~vector_feasible(entry_recovered,cfg,opts.feasibility_tol)
    error('Rank-continuation SDR-SCA entry is not feasible.');
end

Hscaled = cell(U,1);
for u = 1:U
    hu = channels.h_stack(:,u);
    Hscaled{u} = hermitianize(hu*hu')/cfg.comm.noise_power;
end

hist = initialize_history(cfg,entry_lifted,entry_recovered, ...
    fisher_scale,power_scale);
sol = initialize_solution(Winitial,entry_lifted,entry_recovered, ...
    fisher_scale,power_scale);
solver_timer = tic;

% Phase 0: solve the unpenalized SDR once. This gives the CRLB lower-bound
% reference and a sensing-oriented covariance from which rank continuation
% starts.
Dzero = repmat({zeros(N,N)},U,1);
[relaxed_ok,Wrelaxed,relaxed_status,relaxed_obj] = ...
    solve_convex_subproblem(cfg,channels,smat,Hscaled,beta2, ...
    fisher_scale,power_scale,Winitial,Dzero,inf,false,opts);

if ~relaxed_ok
    sol.cvx_status = char(relaxed_status);
    sol.message = sprintf('Unpenalized SDR subproblem failed with status %s.', ...
        relaxed_status);
    sol.runtime_s = toc(solver_timer);
    sol.hist = hist;
    return;
end

Wcurrent = Wrelaxed;
current_lifted = evaluate_lifted_solution(cfg,channels,smat,Wcurrent);
current_recovered = evaluate_recovered(cfg,channels,smat,Wcurrent);
hist = append_history(hist,0,"relaxation",relaxed_status,relaxed_obj, ...
    current_lifted,current_recovered,inf,toc(solver_timer), ...
    cfg,opts,power_scale);
sol.num_solved_iterations = 1;
sol.cvx_status = char(relaxed_status);
sol.relaxed_Wcell = Wrelaxed;
sol.relaxed_metrics = current_lifted;

best = make_incumbent(Winitial,entry_lifted,entry_recovered,0,"entry");
best = update_incumbent(best,Wcurrent,current_lifted,current_recovered, ...
    0,"relaxation",cfg,opts,power_scale);

final_refine_count = 0;
previous_rank_gap = current_lifted.rank_gap_sum/power_scale;
previous_crlb = current_lifted.crlb;
termination_message = 'maximum rank-continuation iteration budget reached';

for iteration = 1:opts.max_continuation_iter
    [Dcell,current_linear_gap] = rank_linearization(Wcurrent,N,U);
    current_linear_gap = max(current_linear_gap/power_scale,0);

    if current_linear_gap <= opts.normalized_rank_tolerance
        desired_budget = opts.normalized_rank_tolerance;
        phase_name = "rank_refinement";
        final_refine_count = final_refine_count+1;
    else
        desired_budget = max(opts.normalized_rank_tolerance, ...
            opts.rank_budget_shrink*current_linear_gap);
        phase_name = "rank_continuation";
    end

    accepted = false;
    attempted_budget = desired_budget;
    status = "not_run";
    Wcandidate = {};
    candidate_obj = NaN;
    for backtrack = 0:opts.max_budget_backtracking
        [subproblem_ok,Wtrial,status,trial_obj] = ...
            solve_convex_subproblem(cfg,channels,smat,Hscaled,beta2, ...
            fisher_scale,power_scale,Wcurrent,Dcell, ...
            attempted_budget,true,opts);
        if subproblem_ok
            trial_lifted = evaluate_lifted_solution( ...
                cfg,channels,smat,Wtrial);
            if lifted_feasible(trial_lifted,cfg,opts.feasibility_tol)
                accepted = true;
                Wcandidate = Wtrial;
                candidate_obj = trial_obj;
                break;
            end
        end

        % Relax only the requested contraction. The previous iterate is
        % feasible at current_linear_gap, so this approaches a feasible
        % subproblem without abandoning rank progress.
        attempted_budget = attempted_budget + ...
            opts.budget_relax_fraction* ...
            (current_linear_gap-attempted_budget);
    end

    if ~accepted
        termination_message = sprintf([ ...
            'rank-continuation subproblem failed at iteration %d ' ...
            '(last status %s)'],iteration,status);
        sol.terminated_by_failure = true;
        break;
    end

    candidate_lifted = evaluate_lifted_solution( ...
        cfg,channels,smat,Wcandidate);
    candidate_recovered = evaluate_recovered( ...
        cfg,channels,smat,Wcandidate);
    hist = append_history(hist,iteration,phase_name,status,candidate_obj, ...
        candidate_lifted,candidate_recovered,attempted_budget, ...
        toc(solver_timer),cfg,opts,power_scale);

    sol.num_solved_iterations = sol.num_solved_iterations+1;
    sol.cvx_status = char(status);
    Wcurrent = Wcandidate;
    current_lifted = candidate_lifted;
    current_recovered = candidate_recovered;
    best = update_incumbent(best,Wcurrent,current_lifted, ...
        current_recovered,iteration,phase_name,cfg,opts,power_scale);

    normalized_rank_gap = current_lifted.rank_gap_sum/power_scale;
    crlb_relative_change = abs(current_lifted.crlb-previous_crlb)/ ...
        max(abs(previous_crlb),realmin);
    rank_relative_change = abs(normalized_rank_gap-previous_rank_gap)/ ...
        max(abs(previous_rank_gap),opts.normalized_rank_tolerance);

    if opts.verbose
        fprintf(['[SDR-RC] iter=%02d phase=%s status=%s CRLB=%.8g ' ...
            'rankGap=%.3e budget=%.3e recFeas=%d\n'],iteration, ...
            phase_name,status,current_lifted.crlb,normalized_rank_gap, ...
            attempted_budget,vector_feasible(current_recovered,cfg, ...
            opts.feasibility_tol));
    end

    if normalized_rank_gap <= opts.normalized_rank_tolerance && ...
            final_refine_count >= opts.min_rank_refine_iter && ...
            crlb_relative_change <= opts.relative_crlb_tolerance && ...
            rank_relative_change <= opts.relative_rank_tolerance
        termination_message = ...
            'rank tolerance and CRLB refinement tolerance reached';
        break;
    end
    if normalized_rank_gap <= opts.normalized_rank_tolerance && ...
            final_refine_count >= opts.max_rank_refine_iter
        termination_message = 'maximum rank-preserving refinement reached';
        break;
    end

    previous_rank_gap = normalized_rank_gap;
    previous_crlb = current_lifted.crlb;
end

sol.current_Wcell = Wcurrent;
sol.current_metrics_lifted = current_lifted;
sol.current_metrics_recovered = current_recovered;
sol.best_iteration = best.iteration;
sol.best_phase = char(best.phase);
sol.Wcell = best.Wcell;
sol.metrics_lifted = best.lifted;
sol.metrics_recovered = best.recovered;
sol.metrics = best.recovered;
sol.normalized_rank_gap = best.lifted.rank_gap_sum/power_scale;
sol.method = 'SDR-SCA rank continuation';
sol.algorithm_version = opts.algorithm_version;
sol.success = best.phase ~= "entry" && ...
    vector_feasible(best.recovered,cfg,opts.feasibility_tol) && ...
    best.lifted.rank_gap_sum/power_scale <= ...
    opts.incumbent_rank_tolerance;
sol.message = termination_message;
if sol.success
    sol.message = sprintf('%s; selected feasible rank-one incumbent at %s %d', ...
        termination_message,best.phase,best.iteration);
end
sol.runtime_s = toc(solver_timer);
hist.total_runtime_s = sol.runtime_s;
hist.termination_message = sol.message;
hist.best_iteration = best.iteration;
hist.best_phase = char(best.phase);
sol.hist = hist;
end

function [ok,Wcell,status,objective_value] = solve_convex_subproblem( ...
    cfg,channels,smat,Hscaled,beta2,fisher_scale,power_scale, ...
    ~,Dcell,rank_budget,enforce_rank,opts)

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;
Gamma = cfg.comm.Gamma;
Pmax = cfg.power.Pmax;
ok = false;
Wcell = {};
status = "not_run";
objective_value = NaN;

try
    cvx_begin sdp
        if opts.cvx_quiet
            cvx_quiet true
        end
        if ~isempty(opts.cvx_solver)
            cvx_solver(opts.cvx_solver)
        end

        variable Wvar(N,N,U) hermitian semidefinite
        variable ta nonnegative
        variable tb nonnegative
        variable kappa_a nonnegative
        variable kappa_b nonnegative
        variable kappa_ab nonnegative
        expression sens_a
        expression sens_b
        expression sens_ab
        expression rank_leakage_normalized

        sens_a = 0;
        sens_b = 0;
        sens_ab = 0;
        rank_leakage_normalized = 0;
        for u = 1:U
            Wu = Wvar(:,:,u);
            sens_a = sens_a + beta2(u)* ...
                real(trace(smat.Aa*Wu))/fisher_scale;
            sens_b = sens_b + beta2(u)* ...
                real(trace(smat.Ab*Wu))/fisher_scale;
            sens_ab = sens_ab + beta2(u)* ...
                real(trace(smat.Aab*Wu))/fisher_scale;
            rank_leakage_normalized = rank_leakage_normalized + ...
                real(trace(Dcell{u}*Wu))/power_scale;
        end
        minimize(inv_pos(ta)+inv_pos(tb))
        subject to
            ta >= opts.normalized_positive_floor;
            tb >= opts.normalized_positive_floor;
            kappa_a >= opts.normalized_positive_floor;
            kappa_b >= opts.normalized_positive_floor;
            kappa_ab >= 0;

            norm([2*kappa_ab;kappa_b-kappa_a+ta]) <= ...
                kappa_b+kappa_a-ta;
            norm([2*kappa_ab;kappa_a-kappa_b+tb]) <= ...
                kappa_a+kappa_b-tb;
            kappa_a <= sens_a;
            kappa_b <= sens_b;
            -kappa_ab <= sens_ab;
            sens_ab <= kappa_ab;

            for u = 1:U
                signal_scaled = real(trace(Hscaled{u}*Wvar(:,:,u)));
                interference_scaled = 0;
                for iuser = 1:U
                    if iuser ~= u
                        interference_scaled = interference_scaled + ...
                            real(trace(Hscaled{u}*Wvar(:,:,iuser)));
                    end
                end
                interference_scaled+1-signal_scaled/Gamma <= 0;
            end
            for m = 1:M
                normalized_power = 0;
                for u = 1:U
                    normalized_power = normalized_power + ...
                        real(trace(smat.Bm{m}*Wvar(:,:,u)))/Pmax;
                end
                normalized_power <= 1;
            end
            if enforce_rank
                rank_leakage_normalized <= rank_budget;
            end
    cvx_end
catch
    status = "failed_exception";
    return;
end

status = string(cvx_status);
if ~(contains(lower(status),"solved") || ...
        contains(lower(status),"inaccurate")) || ...
        any(~isfinite(Wvar(:)))
    return;
end

Wcell = cell(U,1);
for u = 1:U
    Wcell{u} = full(hermitianize(Wvar(:,:,u)));
end
objective_value = cvx_optval;
ok = true;
end

function [Dcell,linear_gap] = rank_linearization(Wcell,N,U)
Dcell = cell(U,1);
linear_gap = 0;
for u = 1:U
    Wu = hermitianize(Wcell{u});
    [V,D] = eig(Wu);
    [~,index] = max(real(diag(D)));
    principal_vector = V(:,index);
    Dcell{u} = hermitianize(eye(N)-principal_vector*principal_vector');
    linear_gap = linear_gap+max(real(trace(Dcell{u}*Wu)),0);
end
end

function hist = initialize_history(cfg,lifted,recovered, ...
    fisher_scale,power_scale)
hist.iter = 0;
hist.phase = "entry";
hist.time_s = 0;
hist.cvx_status = "initial";
hist.surrogate_obj = NaN;
hist.crlb = lifted.crlb;
hist.lifted_crlb = lifted.crlb;
hist.true_ratio = lifted.crlb_ratio;
hist.recovered_crlb = recovered.crlb;
hist.rank_gap_sum = lifted.rank_gap_sum;
hist.normalized_rank_gap = lifted.rank_gap_sum/power_scale;
hist.rank_budget = NaN;
hist.recovered_min_sinr = min(recovered.sinr);
hist.recovered_max_power = max(recovered.ap_power);
hist.lifted_feasible = true;
hist.recovered_feasible = true;
hist.fisher_scale = fisher_scale;
hist.power_scale = power_scale;
hist.snapshot_iter = 0;
hist.z_snapshots = complex_matrix_to_real_vector( ...
    vectorize_wcell_to_matrix(cfg,recovered.wcell));
end

function hist = append_history(hist,iteration,phase,status,obj,lifted, ...
    recovered,rank_budget,time_s,cfg,opts,power_scale)
hist.iter(end+1) = iteration;
hist.phase(end+1) = string(phase);
hist.time_s(end+1) = time_s;
hist.cvx_status(end+1) = string(status);
hist.surrogate_obj(end+1) = obj;
hist.crlb(end+1) = lifted.crlb;
hist.lifted_crlb(end+1) = lifted.crlb;
hist.true_ratio(end+1) = lifted.crlb_ratio;
hist.recovered_crlb(end+1) = recovered.crlb;
hist.rank_gap_sum(end+1) = lifted.rank_gap_sum;
hist.normalized_rank_gap(end+1) = lifted.rank_gap_sum/power_scale;
hist.rank_budget(end+1) = rank_budget;
hist.recovered_min_sinr(end+1) = min(recovered.sinr);
hist.recovered_max_power(end+1) = max(recovered.ap_power);
hist.lifted_feasible(end+1) = ...
    lifted_feasible(lifted,cfg,opts.feasibility_tol);
hist.recovered_feasible(end+1) = ...
    vector_feasible(recovered,cfg,opts.feasibility_tol);
hist.snapshot_iter(end+1) = hist.snapshot_iter(end)+1;
hist.z_snapshots(:,end+1) = complex_matrix_to_real_vector( ...
    vectorize_wcell_to_matrix(cfg,recovered.wcell));
end

function sol = initialize_solution(Wcell,lifted,recovered, ...
    fisher_scale,power_scale)
sol.success = false;
sol.cvx_status = 'not_run';
sol.message = 'no SDR subproblem solved';
sol.Wcell = Wcell;
sol.metrics_lifted = lifted;
sol.metrics_recovered = recovered;
sol.num_solved_iterations = 0;
sol.terminated_by_failure = false;
sol.fisher_scale = fisher_scale;
sol.power_scale = power_scale;
end

function incumbent = make_incumbent(Wcell,lifted,recovered,iteration,phase)
incumbent.Wcell = Wcell;
incumbent.lifted = lifted;
incumbent.recovered = recovered;
incumbent.iteration = iteration;
incumbent.phase = string(phase);
end

function incumbent = update_incumbent(incumbent,Wcell,lifted,recovered, ...
    iteration,phase,cfg,opts,power_scale)
strict_rank_one = lifted.rank_gap_sum/power_scale <= ...
    opts.incumbent_rank_tolerance;
candidate_ok = strict_rank_one && ...
    vector_feasible(recovered,cfg,opts.feasibility_tol);
if candidate_ok && recovered.crlb < incumbent.recovered.crlb
    incumbent = make_incumbent(Wcell,lifted,recovered,iteration,phase);
end
end

function metrics = evaluate_recovered(cfg,channels,smat,Wcell)
wcell = recover_beamformers_from_W(cfg,Wcell);
metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
metrics.wcell = wcell;
end

function tf = lifted_feasible(metrics,cfg,tolerance)
tf = all(isfinite(metrics.sinr)) && all(isfinite(metrics.ap_power)) && ...
    isfinite(metrics.crlb) && metrics.crlb > 0 && ...
    all(metrics.sinr >= cfg.comm.Gamma*(1-tolerance)) && ...
    all(metrics.ap_power <= cfg.power.Pmax*(1+tolerance));
end

function tf = vector_feasible(metrics,cfg,tolerance)
tf = all(isfinite(metrics.sinr)) && all(isfinite(metrics.ap_power)) && ...
    isfinite(metrics.crlb) && metrics.crlb > 0 && ...
    all(metrics.sinr >= cfg.comm.Gamma*(1-tolerance)) && ...
    all(metrics.ap_power <= cfg.power.Pmax*(1+tolerance));
end

function scale = sensing_upper_scale(cfg,smat,beta2)
M = cfg.sys.M;
Nt = cfg.sys.Nt;
Pmax = cfg.power.Pmax;
upper_a = 0;
upper_b = 0;
upper_ab = 0;
for m = 1:M
    indices = block_indices(m,Nt);
    upper_a = upper_a + ...
        max(real(eig(hermitianize(smat.Aa(indices,indices)))));
    upper_b = upper_b + ...
        max(real(eig(hermitianize(smat.Ab(indices,indices)))));
    upper_ab = upper_ab + norm(smat.Aab(indices,indices),2);
end
scale = max(beta2)*Pmax*max([upper_a,upper_b,upper_ab,realmin]);
end

function opts = fill_options(cfg,opts)
defaults.max_continuation_iter = 20;
defaults.rank_budget_shrink = 0.25;
defaults.max_budget_backtracking = 4;
defaults.budget_relax_fraction = 0.5;
defaults.normalized_rank_tolerance = 1e-8;
defaults.incumbent_rank_tolerance = 1e-6;
defaults.min_rank_refine_iter = 2;
defaults.max_rank_refine_iter = 5;
defaults.relative_crlb_tolerance = cfg.opt.tol;
defaults.relative_rank_tolerance = 1e-3;
defaults.normalized_positive_floor = 1e-9;
defaults.feasibility_tol = 2e-4;
defaults.cvx_quiet = cfg.opt.cvx_quiet;
defaults.cvx_solver = cfg.opt.cvx_solver;
defaults.verbose = cfg.opt.verbose;
defaults.algorithm_version = 'blockwise_rank_budget_continuation_v1';
names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts,name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
end
