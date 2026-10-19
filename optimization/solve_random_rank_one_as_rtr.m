function [sol,hist] = solve_random_rank_one_as_rtr( ...
    cfg,channels,smat,init_wcell,rtr_opts)
%SOLVE_RANDOM_RANK_ONE_AS_RTR Random rank-one point followed by AS RTR.
%
% The supplied random rank-one beamformer is passed directly to AS RTR. No
% Phase-I optimization or feasibility restoration is applied.

if nargin < 5
    rtr_opts = struct();
end
if ~isfield(rtr_opts,'max_iter')
    rtr_opts.max_iter = cfg.costar.rtr_max_iter;
end
if ~isfield(rtr_opts,'seed')
    rtr_opts.seed = cfg.seed+cfg.costar.rtr_seed_offset;
end

timer_total = tic;
timer_phase1 = tic;
sol_phase1 = make_phase1_solution( ...
    cfg,channels,smat,init_wcell,'Random rank-one');
time_phase1 = toc(timer_phase1);

timer_rtr = tic;
[sol_rtr,hist_rtr] = solve_active_rtr_phase2( ...
    cfg,channels,smat,sol_phase1.wcell,rtr_opts);
time_rtr = toc(timer_rtr);

sol = sol_rtr;
sol.phase1_solution = sol_phase1;
sol.total_runtime = toc(timer_total);
sol.phase1_runtime = time_phase1;
sol.rtr_runtime = time_rtr;
sol.method = 'Random rank-one AS RTR';
sol.phase1_method = 'Random rank-one';
sol.phase2_method = sol_rtr.phase2_method;

hist.phase1 = struct('method','Random rank-one','crlb', ...
    sol_phase1.metrics.crlb,'time_s',time_phase1);
hist.rtr = hist_rtr;
hist.total_recorded_points = 1+numel(hist_rtr.crlb);
hist.total_runtime = sol.total_runtime;
hist.phase1_runtime = time_phase1;
hist.rtr_runtime = time_rtr;
hist.method = sol.method;
hist.phase2_method = sol.phase2_method;
sol.hist = hist;
end

function sol = make_phase1_solution(cfg,channels,smat,wcell,method)
Wmat = vectorize_wcell_to_matrix(cfg,wcell);
sol.wcell = wcell;
sol.Wmat = Wmat;
sol.z = complex_matrix_to_real_vector(Wmat);
sol.metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
sol.success = true;
sol.message = [method ' point supplied directly'];
sol.method = method;
end
