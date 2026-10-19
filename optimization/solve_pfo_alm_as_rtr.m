function [sol,hist] = solve_pfo_alm_as_rtr( ...
    cfg,channels,smat,init_wcell,rtr_opts)
%SOLVE_PFO_ALM_AS_RTR PFO-ALM followed by active-set RTR refinement.
%
% Phase I runs the unperturbed projected first-order ALM. Its returned
% beamformer is passed directly to the same Phase-II AS RTR routine used by
% COSTAR. No feasibility restoration or alternative initialization is
% inserted between the two phases.

if nargin < 5
    rtr_opts = struct();
end
if ~isfield(rtr_opts,'max_iter')
    rtr_opts.max_iter = cfg.costar.rtr_max_iter;
end
if ~isfield(rtr_opts,'seed')
    rtr_opts.seed = cfg.seed+cfg.costar.rtr_seed_offset;
end

% Guard against a caller carrying PPFO settings into this pure PFO phase.
cfg_pfo = cfg;
cfg_pfo.ppfo.enabled = false;

timer_total = tic;
timer_phase1 = tic;
[sol_phase1,hist_phase1] = solve_alm_first_order_projected_power( ...
    cfg_pfo,channels,smat,init_wcell);
time_phase1 = toc(timer_phase1);

timer_rtr = tic;
[sol_rtr,hist_rtr] = polish_beamformer_active_rtr( ...
    cfg,channels,smat,sol_phase1.wcell,rtr_opts);
time_rtr = toc(timer_rtr);

sol = sol_rtr;
sol.phase1_solution = sol_phase1;
sol.total_runtime = toc(timer_total);
sol.phase1_runtime = time_phase1;
sol.rtr_runtime = time_rtr;
sol.method = 'PFO-ALM AS RTR';
sol.phase1_method = 'PFO-ALM';
sol.phase2_method = 'AS RTR';

hist.phase1 = hist_phase1;
hist.rtr = hist_rtr;
hist.total_recorded_points = numel(hist_phase1.crlb)+numel(hist_rtr.crlb);
hist.total_runtime = sol.total_runtime;
hist.phase1_runtime = time_phase1;
hist.rtr_runtime = time_rtr;
hist.method = sol.method;
sol.hist = hist;
end
