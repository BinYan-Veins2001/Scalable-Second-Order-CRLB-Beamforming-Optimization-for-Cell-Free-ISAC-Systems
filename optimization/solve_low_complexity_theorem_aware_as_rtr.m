function [sol,hist] = solve_low_complexity_theorem_aware_as_rtr( ...
    cfg,channels,smat,rtr_opts)
%SOLVE_LOW_COMPLEXITY_THEOREM_AWARE_AS_RTR No-CVX initializer plus AS RTR.

if nargin < 4
    rtr_opts = struct();
end
if ~isfield(rtr_opts,'max_iter')
    rtr_opts.max_iter = cfg.costar.rtr_max_iter;
end
if ~isfield(rtr_opts,'seed')
    rtr_opts.seed = cfg.seed+cfg.costar.rtr_seed_offset;
end
if isfield(cfg.costar,'low_complexity_init')
    initializer_opts = cfg.costar.low_complexity_init;
else
    initializer_opts = struct();
end

timer_total = tic;
timer_phase1 = tic;
[phase1_wcell,initializer_info] = ...
    build_low_complexity_theorem_aware_initializer( ...
    cfg,channels,smat,initializer_opts);
sol_phase1 = make_phase1_solution( ...
    cfg,channels,smat,phase1_wcell,initializer_info.method);
sol_phase1.initializer_info = initializer_info;
time_phase1 = toc(timer_phase1);

timer_rtr = tic;
[sol_rtr,hist_rtr] = solve_active_rtr_phase2( ...
    cfg,channels,smat,sol_phase1.wcell,rtr_opts);
time_rtr = toc(timer_rtr);

sol = sol_rtr;
sol.phase1_solution = sol_phase1;
sol.initializer_info = initializer_info;
sol.initializer_entry_feasible = initializer_info.entry_feasible;
sol.initializer_certificate_passed = ...
    initializer_info.feasibility_entry_certificate;
sol.total_runtime = toc(timer_total);
sol.phase1_runtime = time_phase1;
sol.rtr_runtime = time_rtr;
sol.method = 'LC theorem-aware ZF AS RTR';
sol.phase1_method = initializer_info.method;
sol.phase2_method = sol_rtr.phase2_method;

hist.phase1 = struct( ...
    'method',initializer_info.method, ...
    'crlb',sol_phase1.metrics.crlb, ...
    'time_s',time_phase1, ...
    'entry_feasible',initializer_info.entry_feasible, ...
    'certificate_passed', ...
        initializer_info.feasibility_entry_certificate, ...
    'initializer_info',initializer_info);
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
sol.message = [method ' construction completed'];
sol.method = method;
end
