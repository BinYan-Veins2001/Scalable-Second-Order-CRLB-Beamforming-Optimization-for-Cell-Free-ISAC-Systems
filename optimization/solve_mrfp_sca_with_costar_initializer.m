function [sol,hist,initializer_info] = ...
    solve_mrfp_sca_with_costar_initializer( ...
    cfg,channels,smat,shared_init_wcell,opts)
%SOLVE_MRFP_SCA_WITH_COSTAR_INITIALIZER Paper-facing MRFP-SCA wrapper.
%
% MRFP-SCA requires a feasible rank-one entry point. This wrapper constructs
% the entry point selected by cfg.costar.initializer_mode, so every figure
% uses the same initialization policy as COSTAR:
%   - global_aggregate_zf: centralized exact ZF plus AP-power projection;
%   - low_complexity_theorem_aware: no-CVX local-ZF initializer;
%   - zf_bf: equation-(35) ZF/RZF plus fixed-direction power allocation;
%   - random_rank_one: the supplied shared point (diagnostic only).
%
% Initialization time is included in sol.total_runtime and in the
% end-to-end time history used by runtime/convergence plots.

if nargin < 4
    shared_init_wcell = [];
end
if nargin < 5
    opts = struct();
end

if isfield(cfg,'costar') && isfield(cfg.costar,'initializer_mode')
    initializer_meta = describe_costar_initializer( ...
        cfg.costar.initializer_mode,cfg);
else
    initializer_meta = describe_costar_initializer('auto',cfg);
end

initializer_timer = tic;
switch string(initializer_meta.mode)
    case "global_aggregate_zf"
        if isfield(cfg.costar,'global_zf_init')
            initializer_options = cfg.costar.global_zf_init;
        else
            initializer_options = struct();
        end
        [init_wcell,initializer_info] = ...
            build_global_aggregate_zf_initializer( ...
            cfg,channels,smat,initializer_options);
    case "low_complexity_theorem_aware"
        if isfield(cfg.costar,'low_complexity_init')
            initializer_options = cfg.costar.low_complexity_init;
        else
            initializer_options = struct();
        end
        [init_wcell,initializer_info] = ...
            build_low_complexity_theorem_aware_initializer( ...
            cfg,channels,smat,initializer_options);
    case "zf_bf"
        [init_wcell,initializer_info] = ...
            build_zf_beamformer(cfg,channels,smat);
    case "random_rank_one"
        if isempty(shared_init_wcell)
            error(['MRFP-SCA with COSTAR random_rank_one mode requires ' ...
                'the same shared rank-one point to be supplied explicitly.']);
        end
        init_wcell = shared_init_wcell;
        initializer_info = struct();
        initializer_info.method = 'shared random rank-one initializer';
        initializer_info.target_feasible = false;
    otherwise
        error('Unsupported COSTAR initializer mode for MRFP-SCA: %s', ...
            initializer_meta.mode);
end
initializer_runtime = toc(initializer_timer);

entry_metrics = evaluate_beamforming_solution( ...
    cfg,channels,smat,init_wcell);
if isfield(cfg,'mrfp') && isfield(cfg.mrfp,'feasibility_tol')
    feasibility_tolerance = cfg.mrfp.feasibility_tol;
else
    feasibility_tolerance = 2e-4;
end
entry_feasible = isfinite(entry_metrics.crlb) && entry_metrics.crlb > 0 && ...
    all(entry_metrics.sinr >= ...
    cfg.comm.Gamma*(1-feasibility_tolerance)) && ...
    all(entry_metrics.ap_power <= ...
    cfg.power.Pmax*(1+feasibility_tolerance));
if ~entry_feasible
    error(['The COSTAR-selected MRFP-SCA initializer is infeasible: ' ...
        'min SINR %.3f dB (target %.3f dB), max AP power %.6g W.'], ...
        pow2db(max(min(entry_metrics.sinr),realmin)), ...
        cfg.comm.Gamma_dB,max(entry_metrics.ap_power));
end

solver_timer = tic;
[sol,hist] = solve_mrfp_sca( ...
    cfg,channels,smat,init_wcell,opts);
solver_runtime = toc(solver_timer);

phase1_solution = struct();
phase1_solution.wcell = init_wcell;
phase1_solution.Wmat = vectorize_wcell_to_matrix(cfg,init_wcell);
phase1_solution.z = complex_matrix_to_real_vector(phase1_solution.Wmat);
phase1_solution.metrics = entry_metrics;
phase1_solution.success = entry_feasible;
phase1_solution.message = initializer_info.method;

initializer_info.configured_mode = initializer_meta.configured_mode;
initializer_info.mode = initializer_meta.mode;
initializer_info.description = initializer_meta.description;
initializer_info.version = initializer_meta.version;
initializer_info.runtime = initializer_runtime;
initializer_info.entry_metrics = entry_metrics;
initializer_info.entry_feasible = entry_feasible;

sol.method = 'MRFP-SCA';
sol.phase1_method = initializer_info.method;
sol.phase2_method = 'MRFP-SCA';
sol.phase1_solution = phase1_solution;
sol.initializer_info = initializer_info;
sol.phase1_runtime = initializer_runtime;
sol.phase2_runtime = solver_runtime;
sol.total_runtime = initializer_runtime+solver_runtime;
sol.message = sprintf('%s; initializer=%s', ...
    sol.message,initializer_info.method);

hist.initializer = struct( ...
    'method',initializer_info.method, ...
    'time_s',initializer_runtime, ...
    'crlb',entry_metrics.crlb, ...
    'entry_feasible',entry_feasible);
if isfield(hist,'cumulative_time_s')
    hist.end_to_end_time_s = initializer_runtime+ ...
        hist.cumulative_time_s(:);
elseif isfield(hist,'time_s')
    hist.end_to_end_time_s = initializer_runtime+hist.time_s(:);
else
    hist.end_to_end_time_s = [initializer_runtime;sol.total_runtime];
end
hist.total_runtime = sol.total_runtime;
hist.method = 'MRFP-SCA';
sol.hist = hist;
end
