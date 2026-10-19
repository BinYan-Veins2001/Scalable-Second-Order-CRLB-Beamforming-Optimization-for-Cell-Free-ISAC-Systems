function [sol,hist] = solve_tr_alm(cfg,channels,smat,init_wcell)
%SOLVE_TR_ALM Projected trust-region ALM (PTR-ALM) baseline.
%
% The default paper mode starts from the same theorem-aware global-ZF point
% as PFO-ALM, uses PFO-ALM as a globalization stage, and then applies a
% bounded PTR refinement. The best feasible stage endpoint is retained, so
% the refinement cannot numerically degrade a feasible PFO incumbent. Set
% cfg.tr.pfo_warmstart_enabled=false to recover the direct legacy PTR path.

if nargin < 4; init_wcell = []; end
use_safeguard = isfield(cfg,'tr') && ...
    isfield(cfg.tr,'pfo_warmstart_enabled') && ...
    logical(cfg.tr.pfo_warmstart_enabled);

if use_safeguard
    cfg_pfo = cfg;
    if isfield(cfg_pfo,'ppfo'); cfg_pfo.ppfo.enabled = false; end
    stage_timer = tic;
    [pfo_sol,pfo_hist] = solve_alm_first_order_projected_power( ...
        cfg_pfo,channels,smat,init_wcell);
    pfo_runtime = toc(stage_timer);

    cfg_tr = configure_costar_tr_stage(cfg);
    cfg_tr.tr.max_outer = min(cfg.tr.refine_max_outer,cfg.alm.max_outer);
    cfg_tr.tr.max_inner = min(cfg.tr.refine_max_inner,cfg.alm.max_inner);
    cfg_tr.tr.delta0 = cfg.tr.refine_delta0;
    cfg_tr.tr.delta_max = cfg.tr.refine_delta_max;
    cfg_tr.tr.gamma_expand = cfg.tr.refine_gamma_expand;
    cfg_tr.tr.gamma_pg_success = cfg.tr.refine_gamma_pg_success;
    stage_timer = tic;
    [ptr_sol,ptr_hist] = solve_alm_trust_region( ...
        cfg_tr,channels,smat,pfo_sol.wcell);
    ptr_runtime = toc(stage_timer);

    [sol,selected_stage] = select_feasible_incumbent( ...
        cfg,channels,smat,pfo_sol,ptr_sol);
    hist = combine_stage_histories(pfo_hist,ptr_hist,pfo_runtime, ...
        ptr_runtime,sol.metrics,selected_stage,cfg);
    sol.success = logical(pfo_sol.success || ptr_sol.success);
    sol.converged_pfo = logical(pfo_sol.success);
    sol.converged_ptr = logical(ptr_sol.success);
    sol.pfo_stage_runtime_s = pfo_runtime;
    sol.ptr_stage_runtime_s = ptr_runtime;
    sol.selected_stage = selected_stage;
    sol.pfo_stage_crlb = pfo_sol.metrics.crlb;
    sol.ptr_stage_crlb = ptr_sol.metrics.crlb;
    sol.safeguard_enabled = true;
    sol.message = sprintf('PTR-ALM; selected %s',selected_stage);
else
    cfg_tr = configure_costar_tr_stage(cfg);
    [sol,hist] = solve_alm_trust_region( ...
        cfg_tr,channels,smat,init_wcell);
    sol.safeguard_enabled = false;
end

if ~isempty(init_wcell)
    initial_Wmat = project_per_ap_power_matrix( ...
        cfg,vectorize_wcell_to_matrix(cfg,init_wcell));
    caller_supplied_initializer = true;
else
    initial_Wmat = sol.initial_Wmat;
    caller_supplied_initializer = false;
end

sol.method = 'PTR-ALM';
sol.max_outer = cfg_tr.tr.max_outer;
sol.max_inner = cfg_tr.tr.max_inner;
sol.initial_Wmat = initial_Wmat;
sol.initial_wcell = matrix_to_wcell(cfg,initial_Wmat);
sol.initial_z = complex_matrix_to_real_vector(initial_Wmat);
sol.caller_supplied_initializer = caller_supplied_initializer;
hist.method = sol.method;
hist.max_outer = cfg_tr.tr.max_outer;
hist.max_inner = cfg_tr.tr.max_inner;
hist.initial_Wmat = initial_Wmat;
hist.caller_supplied_initializer = caller_supplied_initializer;
sol.hist = hist;
end

function value = ternary(condition,a,b)
if condition; value=a; else; value=b; end
end

function [sol,stage] = select_feasible_incumbent( ...
    cfg,channels,smat,pfo_sol,ptr_sol)
pfo_metrics = evaluate_beamforming_solution( ...
    cfg,channels,smat,pfo_sol.wcell);
ptr_metrics = evaluate_beamforming_solution( ...
    cfg,channels,smat,ptr_sol.wcell);
tol = cfg.tr.incumbent_feas_tol;
pfo_feasible = endpoint_feasible(cfg,pfo_metrics,tol);
ptr_feasible = endpoint_feasible(cfg,ptr_metrics,tol);
if pfo_feasible && (~ptr_feasible || pfo_metrics.crlb <= ptr_metrics.crlb)
    sol = pfo_sol; sol.metrics = pfo_metrics; stage = 'PFO safeguard';
else
    sol = ptr_sol; sol.metrics = ptr_metrics; stage = 'PTR refinement';
end
sol.pfo_endpoint_feasible = pfo_feasible;
sol.ptr_endpoint_feasible = ptr_feasible;
end

function yes = endpoint_feasible(cfg,metrics,tol)
yes = all(metrics.sinr >= cfg.comm.Gamma*(1-tol)) && ...
    all(metrics.ap_power <= cfg.power.Pmax*(1+tol)) && ...
    isfinite(metrics.crlb);
end

function hist = combine_stage_histories(pfo,ptr,pfo_runtime,ptr_runtime, ...
    selected_metrics,selected_stage,cfg)
hist = struct();
hist.time_s = [pfo.time_s(:); pfo_runtime+ptr.time_s(:)];
hist.crlb = [pfo.crlb(:); ptr.crlb(:)];
hist.max_sinr_violation = [pfo.max_sinr_violation(:); ...
    ptr.max_sinr_violation(:)];
hist.max_power_violation = [pfo.max_power_violation(:); ...
    ptr.max_power_violation(:)];
hist.phase = [repmat({ 'PFO globalization' },numel(pfo.crlb),1); ...
    repmat({ 'PTR refinement' },numel(ptr.crlb),1)];
hist.time_s(end+1,1) = pfo_runtime+ptr_runtime;
hist.crlb(end+1,1) = selected_metrics.crlb;
hist.max_sinr_violation(end+1,1) = max(0,max( ...
    cfg.comm.Gamma(:)-selected_metrics.sinr(:)) ./ ...
    max(cfg.comm.Gamma(:),realmin));
hist.max_power_violation(end+1,1) = max(0,max( ...
    selected_metrics.ap_power(:)/max(cfg.power.Pmax,realmin)-1));
hist.phase{end+1,1} = selected_stage;
hist.stage1_pfo = pfo;
hist.stage2_ptr = ptr;
hist.method = 'PTR-ALM';
hist.total_runtime_s = pfo_runtime+ptr_runtime;
hist.ptr_start_time_s = pfo_runtime;
end
