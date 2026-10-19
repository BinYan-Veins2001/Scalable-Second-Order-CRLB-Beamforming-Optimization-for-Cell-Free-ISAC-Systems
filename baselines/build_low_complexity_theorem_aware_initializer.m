function [wcell,info] = build_low_complexity_theorem_aware_initializer( ...
    cfg,channels,smat,opts)
%BUILD_LOW_COMPLEXITY_THEOREM_AWARE_INITIALIZER No-CVX COSTAR warm start.
%
% This independent prototype is intended for the small-system regime
% Nt >= U, where the local equation-(35) ZF directions exist.  It uses:
%
%  1) exact per-AP local-ZF directions;
%  2) a closed-form SINR power floor;
%  3) strict SINR and AP-power margins motivated by the Phase-II
%     initialization certificate; and
%  4) a few projected CRLB-gradient redistribution steps that never reduce
%     any power coefficient below the closed-form feasible floor.
%
% No CVX problem is solved.  The construction cost is dominated by M
% inversions of U-by-U Gram matrices and O(max_refine_iter*M*U) scalar power
% updates, in addition to CRLB gradient evaluations.
%
% Important theoretical scope:
%   - the returned point is directly checkable for SINR/AP-power feasibility;
%   - EVALUATE_AS_RTR_INITIALIZATION_CERTIFICATE can certify any required
%     local normal restoration;
%   - feasibility alone does not prove membership in the attraction basin of
%     a particular second-order KKT phase orbit.

if nargin < 4
    opts = struct();
end
opts = fill_defaults(opts);

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
Pmax = cfg.power.Pmax;

if Nt < U
    error(['The theorem-aware closed-form local-ZF initializer requires ' ...
        'Nt >= U. Use a separate centralized/RZF extension when Nt < U.']);
end

timer_total = tic;
directions = cell(M,U);
gain = zeros(M,U);

timer_direction = tic;
for m = 1:M
    Hm = zeros(U,Nt);
    for u = 1:U
        Hm(u,:) = channels.h{m,u}';
    end
    gram = hermitianize(Hm*Hm');
    gram_scale = max(real(trace(gram))/max(U,1),realmin);
    gram_inverse = pinv(gram,1e-12*gram_scale);
    raw_direction = Hm'*gram_inverse;
    for u = 1:U
        d = raw_direction(:,u);
        if norm(d) <= realmin
            error('Degenerate local-ZF direction at AP %d, user %d.',m,u);
        end
        d = d/norm(d);
        local_gain = channels.h{m,u}'*d;
        d = d*exp(-1i*angle(local_gain));
        directions{m,u} = d;
        gain(m,u) = max(real(channels.h{m,u}'*d),0);
    end
end
direction_runtime = toc(timer_direction);

% With exact local ZF, all multiuser interference is zero and
%
%   SINR_u = (gain(:,u)'*F(:,u))^2 / sigma^2.
%
% The minimum Euclidean-norm nonnegative vector attaining the strict target
% is therefore available in closed form.
strict_gamma = cfg.comm.Gamma*(1+opts.sinr_margin_fraction);
required_amplitude = sqrt(strict_gamma*cfg.comm.noise_power);
Ffloor = zeros(M,U);
for u = 1:U
    gain_norm_sq = sum(gain(:,u).^2);
    if gain_norm_sq <= realmin
        error('User %d has zero aggregate local-ZF gain.',u);
    end
    Ffloor(:,u) = required_amplitude*gain(:,u)/gain_norm_sq;
end

power_cap = (1-opts.power_margin_fraction)*Pmax;
floor_ap_power = sum(Ffloor.^2,2);
closed_form_power_feasible = all(floor_ap_power <= ...
    power_cap*(1+opts.feasibility_tolerance));

if closed_form_power_feasible
    % Start at the closed-form minimum-power point.  This keeps the SINR
    % constraints near active, which is empirically important for entering
    % the same AS-RTR basin as the SOCP power allocator.
    F = Ffloor;
else
    % Return a power-bounded diagnostic point.  INFO marks it infeasible so
    % the caller can select a fallback without silently entering AS RTR.
    F = Ffloor;
    for m = 1:M
        row_norm = norm(F(m,:));
        if row_norm > sqrt(power_cap)
            F(m,:) = F(m,:)*sqrt(power_cap)/row_norm;
        end
    end
end

timer_refine = tic;
[current_obj,current_crlb,current_grad] = objective_and_gradient(F);
objective_history = current_obj;
crlb_history = current_crlb;
accepted_history = false(0,1);
radius = opts.initial_radius*sqrt(Pmax);

if closed_form_power_feasible && opts.max_refine_iter > 0 && ~isempty(smat)
    for refine_iter = 1:opts.max_refine_iter
        grad_norm = norm(current_grad,'fro');
        if ~isfinite(grad_norm) || grad_norm <= opts.gradient_tolerance
            break;
        end
        normalized_grad = current_grad/grad_norm;
        local_radius = radius;
        accepted = false;
        for backtrack_iter = 1:opts.max_backtracking
            trial = F-local_radius*normalized_grad;
            trial = project_zf_power_feasible_dykstra( ...
                trial,gain,required_amplitude,sqrt(power_cap), ...
                opts.projection_cycles);
            [trial_obj,trial_crlb,trial_grad] = ...
                objective_and_gradient(trial);
            if isfinite(trial_obj) && trial_obj <= current_obj - ...
                    opts.objective_tolerance*max(1,abs(current_obj))
                F = trial;
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
        accepted_history(end+1,1) = accepted; %#ok<AGROW>
        objective_history(end+1,1) = current_obj; %#ok<AGROW>
        crlb_history(end+1,1) = current_crlb; %#ok<AGROW>
        if ~accepted || local_radius <= opts.min_radius*sqrt(Pmax)
            break;
        end
    end
end
refine_runtime = toc(timer_refine);

wcell = assemble_wcell(F);
metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
entry_feasible = all(metrics.sinr >= ...
    cfg.comm.Gamma*(1-opts.feasibility_tolerance)) && ...
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

info.method = 'LC theorem-aware ZF initializer';
info.wcell = wcell;
info.directions = directions;
info.effective_gain = gain;
info.power_floor = Ffloor;
info.power_variable = F;
info.allocated_power = F.^2;
info.floor_ap_power = floor_ap_power;
info.closed_form_power_feasible = closed_form_power_feasible;
info.entry_feasible = entry_feasible;
info.strict_gamma = strict_gamma;
info.power_cap = power_cap;
info.metrics = metrics;
info.objective_history = objective_history;
info.crlb_history = crlb_history;
info.accepted_history = accepted_history;
info.num_accepted_refine_steps = nnz(accepted_history);
info.normal_certificate = normal_certificate;
if opts.evaluate_normal_certificate
    % If no SINR constraint is near active, no normal restoration is
    % required. Otherwise use the finite-step Kantorovich certificate.
    info.feasibility_entry_certificate = entry_feasible && ...
        (normal_certificate.active_sinr_count == 0 || ...
        normal_certificate.kantorovich_finite_step_certificate);
else
    info.feasibility_entry_certificate = entry_feasible;
end
info.direction_runtime = direction_runtime;
info.refine_runtime = refine_runtime;
info.total_runtime = toc(timer_total);
info.options = opts;

    function [obj,crlb,gradF] = objective_and_gradient(Fvalue)
        wvalue = assemble_wcell(Fvalue);
        Wvalue = vectorize_wcell_to_matrix(cfg,wvalue);
        cfg_objective = cfg;
        cfg_objective.alm.objective_type = 'log_crlb_ratio';
        [obj,gradW,gout] = compute_crlb_objective_grad( ...
            cfg_objective,smat,Wvalue);
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
                wvalue{m1,u1} = ...
                    Fvalue(m1,u1)*directions{m1,u1};
            end
        end
    end
end

function F = project_zf_power_feasible_dykstra( ...
        V,gain,required_amplitude,power_radius,num_cycles)
% Dykstra projection onto nonnegativity, AP balls, and ZF SINR halfspaces.
[M,U] = size(V);
num_sets = 1+M+U;
correction = zeros(M,U,num_sets);
F = real(V);
for cycle = 1:num_cycles
    set_index = 1;
    Y = F+correction(:,:,set_index);
    Fnew = max(Y,0);
    correction(:,:,set_index) = Y-Fnew;
    F = Fnew;

    for m = 1:M
        set_index = set_index+1;
        Y = F+correction(:,:,set_index);
        Fnew = Y;
        row_norm = norm(Y(m,:));
        if row_norm > power_radius
            Fnew(m,:) = Y(m,:)*power_radius/row_norm;
        end
        correction(:,:,set_index) = Y-Fnew;
        F = Fnew;
    end

    for u = 1:U
        set_index = set_index+1;
        Y = F+correction(:,:,set_index);
        Fnew = Y;
        deficit = required_amplitude-gain(:,u)'*Y(:,u);
        if deficit > 0
            Fnew(:,u) = Y(:,u) + ...
                deficit*gain(:,u)/max(sum(gain(:,u).^2),realmin);
        end
        correction(:,:,set_index) = Y-Fnew;
        F = Fnew;
    end
end
F = max(real(F),0);
end

function opts = fill_defaults(opts)
if ~isfield(opts,'sinr_margin_fraction')
    opts.sinr_margin_fraction = 1e-3;
end
if ~isfield(opts,'power_margin_fraction')
    opts.power_margin_fraction = 1e-3;
end
if ~isfield(opts,'max_refine_iter'); opts.max_refine_iter=20; end
if ~isfield(opts,'max_backtracking'); opts.max_backtracking=8; end
if ~isfield(opts,'initial_radius'); opts.initial_radius=0.25; end
if ~isfield(opts,'max_radius'); opts.max_radius=1; end
if ~isfield(opts,'min_radius'); opts.min_radius=1e-6; end
if ~isfield(opts,'radius_shrink'); opts.radius_shrink=0.5; end
if ~isfield(opts,'radius_expand'); opts.radius_expand=1.5; end
if ~isfield(opts,'gradient_tolerance'); opts.gradient_tolerance=1e-8; end
if ~isfield(opts,'objective_tolerance'); opts.objective_tolerance=1e-9; end
if ~isfield(opts,'projection_cycles'); opts.projection_cycles=60; end
if ~isfield(opts,'feasibility_tolerance')
    opts.feasibility_tolerance=1e-4;
end
if ~isfield(opts,'evaluate_normal_certificate')
    opts.evaluate_normal_certificate=true;
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
