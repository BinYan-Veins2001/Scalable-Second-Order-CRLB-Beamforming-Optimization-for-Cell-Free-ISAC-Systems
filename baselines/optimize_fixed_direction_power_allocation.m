function [wcell,info] = optimize_fixed_direction_power_allocation( ...
    cfg,channels,smat,directions,opts)
%OPTIMIZE_FIXED_DIRECTION_POWER_ALLOCATION Power allocation for fixed beams.
%
% Let d_{m,u} be a unit-norm fixed beam direction and let
%
%                    w_{m,u} = F_{m,u} d_{m,u},  F_{m,u} >= 0,
%
% where F_{m,u}=sqrt(P_{m,u}).  This function implements the power-variable
% structure used by (35), (36), and Appendix C of the reference paper.
%
% The communication SINR and per-AP power constraints are imposed as SOC
% constraints.  A feasible minimum-power point is found first.  With the
% directions fixed, a proximal SCA/SOCP step then minimizes a first-order
% model of the project's physical CRLB while preserving those constraints.
%
% If the requested common SINR is not achievable by the supplied fixed
% directions, a bisection finds the largest achievable common SINR.  The
% returned point is then deliberately marked target-infeasible in INFO.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(cfg,smat,opts);

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
Pmax = cfg.power.Pmax;

if ~iscell(directions) || ~isequal(size(directions),[M,U])
    error('directions must be an M-by-U cell array.');
end

% Normalize every direction and rotate its irrelevant scalar phase so that
% the intended local channel gain is real and nonnegative.  This makes the
% desired coherent sum an affine real quantity in the SOC constraints.
for m = 1:M
    for u = 1:U
        d = directions{m,u}(:);
        if numel(d) ~= Nt || norm(d) <= realmin
            error('Invalid fixed direction at AP %d, user %d.',m,u);
        end
        d = d/norm(d);
        desired_gain = channels.h{m,u}'*d;
        if abs(desired_gain) > realmin
            d = d*exp(-1i*angle(desired_gain));
        end
        directions{m,u} = d;
    end
end

% Normalized effective channels.  Dividing by sqrt(noise) avoids the
% 1e-8--1e-4 coefficient range that otherwise causes unnecessary CVX
% numerical warnings in high-SINR tests.
noise_sqrt = sqrt(max(cfg.comm.noise_power,realmin));
effective_channel = complex(zeros(U,U,M));
for rx_user = 1:U
    for stream = 1:U
        for m = 1:M
            effective_channel(rx_user,stream,m) = ...
                channels.h{m,rx_user}'*directions{m,stream}/noise_sqrt;
        end
    end
end

gamma_requested = cfg.comm.Gamma;
gamma_soc = gamma_requested*(1+opts.sinr_margin);
[target_soc_feasible,F,status] = ...
    solve_fixed_direction_min_power_socp( ...
    effective_channel,gamma_soc,cfg,opts);

used_gamma = gamma_requested;
target_direction_feasible = target_soc_feasible;
if ~target_soc_feasible
    % The direction family itself may be unable to meet the requested SINR
    % (especially SO-BF).  Find its maximum common-SINR operating point so
    % the returned baseline remains meaningful and reproducible.
    [F,used_gamma,status] = ...
        solve_fixed_direction_max_common_sinr( ...
        effective_channel,gamma_requested,cfg,opts);
end

if isempty(F) || any(~isfinite(F(:)))
    F = sqrt(Pmax/U)*ones(M,U);
    status = [char(status) '; fallback_equal_power'];
    used_gamma = 0;
end

initial_wcell = assemble_wcell(F);
[current_obj,current_crlb,current_grad] = objective_and_gradient(F);
crlb_history = current_crlb;
objective_history = current_obj;
accepted_history = false(0,1);
radius_history = zeros(0,1);
cvx_status_history = string(status);

if opts.optimize_sensing && ~isempty(smat) && used_gamma > 0
    radius = opts.initial_radius*sqrt(Pmax);
    for sca_iter = 1:opts.max_sca_iter
        grad_norm = norm(current_grad,'fro');
        if ~isfinite(grad_norm) || grad_norm <= opts.gradient_tol
            break;
        end
        normalized_grad = current_grad/grad_norm;

        accepted = false;
        local_radius = radius;
        trial_status = "not_run";
        for backtrack_iter = 1:opts.max_backtracking
            [trial_ok,Ftrial,trial_status] = ...
                solve_fixed_direction_proximal_socp( ...
                effective_channel,used_gamma,normalized_grad,F, ...
                local_radius,cfg,opts);
            if ~trial_ok
                local_radius = local_radius*opts.radius_shrink;
                continue;
            end

            [trial_obj,trial_crlb,trial_grad] = ...
                objective_and_gradient(Ftrial);
            sufficient_decrease = trial_obj <= current_obj - ...
                opts.objective_tol*max(1,abs(current_obj));
            if isfinite(trial_obj) && sufficient_decrease
                F = Ftrial;
                current_obj = trial_obj;
                current_crlb = trial_crlb;
                current_grad = trial_grad;
                radius = min(opts.radius_expand*local_radius, ...
                    opts.max_radius*sqrt(Pmax));
                accepted = true;
                break;
            end
            local_radius = local_radius*opts.radius_shrink;
        end

        crlb_history(end+1,1) = current_crlb; %#ok<AGROW>
        objective_history(end+1,1) = current_obj; %#ok<AGROW>
        accepted_history(end+1,1) = accepted; %#ok<AGROW>
        radius_history(end+1,1) = local_radius; %#ok<AGROW>
        cvx_status_history(end+1,1) = string(trial_status); %#ok<AGROW>

        if ~accepted || local_radius <= opts.min_radius*sqrt(Pmax)
            break;
        end
    end
end

wcell = assemble_wcell(F);
if isempty(smat)
    metrics = evaluate_communication_metrics(cfg,channels,wcell);
else
    metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
end
target_feasible = all(metrics.sinr >= gamma_requested*(1-opts.check_tol)) && ...
    all(metrics.ap_power <= Pmax*(1+opts.check_tol));

info.method = opts.method_name;
info.power_variable = F;
info.allocated_power = F.^2;
info.requested_gamma = gamma_requested;
info.requested_gamma_dB = pow2db(max(gamma_requested,realmin));
info.soc_gamma = used_gamma;
info.soc_gamma_dB = pow2db(max(used_gamma,realmin));
info.direction_target_feasible = target_direction_feasible;
info.target_feasible = target_feasible;
info.cvx_status = char(status);
info.initial_wcell = initial_wcell;
info.initial_crlb = crlb_history(1);
info.final_crlb = metrics.crlb;
info.final_min_sinr = min(metrics.sinr);
info.final_min_sinr_dB = pow2db(max(info.final_min_sinr,realmin));
info.final_max_ap_power = max(metrics.ap_power);
info.crlb_history = crlb_history;
info.objective_history = objective_history;
info.accepted_history = accepted_history;
info.radius_history = radius_history;
info.cvx_status_history = cvx_status_history;
info.num_sca_iterations = numel(accepted_history);
info.num_accepted_sca_steps = nnz(accepted_history);
info.options = opts;

    function [obj,crlb,gradF] = objective_and_gradient(Fvalue)
        wvalue = assemble_wcell(Fvalue);
        Wvalue = vectorize_wcell_to_matrix(cfg,wvalue);
        if isempty(smat)
            obj = NaN;
            crlb = NaN;
            gradF = zeros(M,U);
            return;
        end
        cfg_objective = cfg;
        if ~isfield(cfg_objective,'alm')
            cfg_objective.alm = struct();
        end
        cfg_objective.alm.objective_type = 'log_crlb_ratio';
        [obj,gradW,gout] = ...
            compute_crlb_objective_grad(cfg_objective,smat,Wvalue);
        crlb = gout.crlb;
        gradF = zeros(M,U);
        for u1 = 1:U
            for m1 = 1:M
                rows = block_indices(m1,Nt);
                gradF(m1,u1) = real( ...
                    directions{m1,u1}'*gradW(rows,u1));
            end
        end
    end

    function wvalue = assemble_wcell(Fvalue)
        wvalue = cell(M,U);
        for m1 = 1:M
            for u1 = 1:U
                wvalue{m1,u1} = Fvalue(m1,u1)*directions{m1,u1};
            end
        end
    end
end

function [ok,Fsol,solver_status] = ...
        solve_fixed_direction_min_power_socp( ...
        effective_channel,gamma,cfg,opts)
M = cfg.sys.M;
U = cfg.sys.U;
Pmax = cfg.power.Pmax;
ok = false;
Fsol = [];
solver_status = 'not_run';
gamma = max(real(gamma),realmin);
power_radius = sqrt(Pmax)*(1-opts.power_margin);
try
    cvx_begin
        if opts.cvx_quiet
            cvx_quiet true
        end
        if ~isempty(opts.cvx_solver)
            cvx_solver(opts.cvx_solver)
        end
        variable Fcvx(M,U) nonnegative
        minimize(sum_square(Fcvx(:)))
        subject to
            for m = 1:M
                norm(Fcvx(m,:),2) <= power_radius;
            end
            for u = 1:U
                coefficient = reshape(effective_channel(u,:,:),[U,M]).';
                response_real = sum(real(coefficient).*Fcvx,1);
                response_imag = sum(imag(coefficient).*Fcvx,1);
                if U == 1
                    1 <= response_real(u)/sqrt(gamma);
                else
                    other_streams = [1:u-1,u+1:U];
                    norm([transpose(response_real(other_streams)); ...
                        transpose(response_imag(other_streams));1],2) ...
                        <= response_real(u)/sqrt(gamma);
                end
            end
    cvx_end
    solver_status = cvx_status;
    Fcandidate = full(Fcvx);
    ok = cvx_status_is_solved(solver_status) && ...
        all(isfinite(Fcandidate(:)));
    if ok
        Fsol = max(real(Fcandidate),0);
    end
catch ME
    cvx_clear
    solver_status = ['CVX error: ' getReport(ME,'extended','hyperlinks','off')];
end
end

function [Fbest,gamma_best,solver_status] = ...
        solve_fixed_direction_max_common_sinr( ...
        effective_channel,gamma_requested,cfg,opts)
requested_dB = pow2db(max(gamma_requested,realmin));
low_dB = min(opts.minimum_gamma_dB, ...
    requested_dB-opts.gamma_search_span_dB);
[low_ok,Fbest,solver_status] = ...
    solve_fixed_direction_min_power_socp( ...
    effective_channel,db2pow(low_dB),cfg,opts);
if ~low_ok
    Fbest = [];
    gamma_best = 0;
    return;
end

high_dB = requested_dB;
for bisect_iter = 1:opts.gamma_bisection_iters
    mid_dB = 0.5*(low_dB+high_dB);
    [mid_ok,Fmid,mid_status] = ...
        solve_fixed_direction_min_power_socp( ...
        effective_channel,db2pow(mid_dB),cfg,opts);
    solver_status = mid_status;
    if mid_ok
        low_dB = mid_dB;
        Fbest = Fmid;
    else
        high_dB = mid_dB;
    end
end
gamma_best = db2pow(low_dB)*(1-opts.infeasible_gamma_backoff);
[backoff_ok,Fbackoff,backoff_status] = ...
    solve_fixed_direction_min_power_socp( ...
    effective_channel,gamma_best,cfg,opts);
solver_status = backoff_status;
if backoff_ok
    Fbest = Fbackoff;
end
end

function [ok,Fsol,solver_status] = ...
        solve_fixed_direction_proximal_socp( ...
        effective_channel,gamma,normalized_grad,Fcenter,radius,cfg,opts)
M = cfg.sys.M;
U = cfg.sys.U;
Pmax = cfg.power.Pmax;
ok = false;
Fsol = [];
solver_status = 'not_run';
gamma = max(real(gamma),realmin);
radius = max(real(radius),realmin);
power_radius = sqrt(Pmax)*(1-opts.power_margin);
try
    cvx_begin
        if opts.cvx_quiet
            cvx_quiet true
        end
        if ~isempty(opts.cvx_solver)
            cvx_solver(opts.cvx_solver)
        end
        variable Fcvx(M,U) nonnegative
        minimize(normalized_grad(:)'*(Fcvx(:)-Fcenter(:)) + ...
            sum_square(Fcvx(:)-Fcenter(:))/(2*radius))
        subject to
            for m = 1:M
                norm(Fcvx(m,:),2) <= power_radius;
            end
            for u = 1:U
                coefficient = reshape(effective_channel(u,:,:),[U,M]).';
                response_real = sum(real(coefficient).*Fcvx,1);
                response_imag = sum(imag(coefficient).*Fcvx,1);
                if U == 1
                    1 <= response_real(u)/sqrt(gamma);
                else
                    other_streams = [1:u-1,u+1:U];
                    norm([transpose(response_real(other_streams)); ...
                        transpose(response_imag(other_streams));1],2) ...
                        <= response_real(u)/sqrt(gamma);
                end
            end
    cvx_end
    solver_status = cvx_status;
    Fcandidate = full(Fcvx);
    ok = cvx_status_is_solved(solver_status) && ...
        all(isfinite(Fcandidate(:)));
    if ok
        Fsol = max(real(Fcandidate),0);
    end
catch ME
    cvx_clear
    solver_status = ['CVX error: ' getReport(ME,'extended','hyperlinks','off')];
end
end

function tf = cvx_status_is_solved(status)
status = lower(strtrim(char(status)));
tf = contains(status,'solved') && ~contains(status,'failed');
end

function metrics = evaluate_communication_metrics(cfg,channels,wcell)
M = cfg.sys.M;
U = cfg.sys.U;
sinr = zeros(U,1);
for u = 1:U
    desired = 0;
    for m = 1:M
        desired = desired+channels.h{m,u}'*wcell{m,u};
    end
    interference = 0;
    for stream = 1:U
        if stream == u
            continue;
        end
        response = 0;
        for m = 1:M
            response = response+channels.h{m,u}'*wcell{m,stream};
        end
        interference = interference+abs(response)^2;
    end
    sinr(u) = abs(desired)^2/ ...
        max(interference+cfg.comm.noise_power,realmin);
end
ap_power = zeros(M,1);
for m = 1:M
    for u = 1:U
        ap_power(m) = ap_power(m)+norm(wcell{m,u})^2;
    end
end
metrics.sinr = sinr;
metrics.ap_power = ap_power;
metrics.crlb = NaN;
end

function opts = fill_defaults(cfg,smat,opts)
if isfield(cfg,'baseline') && isfield(cfg.baseline,'power_alloc')
    defaults = cfg.baseline.power_alloc;
else
    defaults = struct();
end
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(opts,fields{k})
        opts.(fields{k}) = defaults.(fields{k});
    end
end

if ~isfield(opts,'method_name'); opts.method_name='Fixed-direction BF'; end
if ~isfield(opts,'optimize_sensing'); opts.optimize_sensing=~isempty(smat); end
if ~isfield(opts,'max_sca_iter'); opts.max_sca_iter=8; end
if ~isfield(opts,'max_backtracking'); opts.max_backtracking=6; end
if ~isfield(opts,'initial_radius'); opts.initial_radius=0.25; end
if ~isfield(opts,'max_radius'); opts.max_radius=1.0; end
if ~isfield(opts,'min_radius'); opts.min_radius=1e-5; end
if ~isfield(opts,'radius_shrink'); opts.radius_shrink=0.5; end
if ~isfield(opts,'radius_expand'); opts.radius_expand=1.5; end
if ~isfield(opts,'gradient_tol'); opts.gradient_tol=1e-8; end
if ~isfield(opts,'objective_tol'); opts.objective_tol=1e-8; end
if ~isfield(opts,'sinr_margin'); opts.sinr_margin=1e-3; end
if ~isfield(opts,'power_margin'); opts.power_margin=1e-6; end
if ~isfield(opts,'check_tol'); opts.check_tol=1e-4; end
if ~isfield(opts,'minimum_gamma_dB'); opts.minimum_gamma_dB=-40; end
if ~isfield(opts,'gamma_search_span_dB'); opts.gamma_search_span_dB=60; end
if ~isfield(opts,'gamma_bisection_iters'); opts.gamma_bisection_iters=10; end
if ~isfield(opts,'infeasible_gamma_backoff')
    opts.infeasible_gamma_backoff=2e-3;
end
if ~isfield(opts,'cvx_quiet')
    opts.cvx_quiet = isfield(cfg,'opt') && ...
        isfield(cfg.opt,'cvx_quiet') && cfg.opt.cvx_quiet;
end
if ~isfield(opts,'cvx_solver')
    if isfield(cfg,'opt') && isfield(cfg.opt,'cvx_solver')
        opts.cvx_solver = cfg.opt.cvx_solver;
    else
        opts.cvx_solver = '';
    end
end
end
