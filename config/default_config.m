function cfg = default_config()
%DEFAULT_CONFIG Centralized configuration for the reproduction package.

cfg.seed = 2026;

%% System size
cfg.sys.M  = 16;        % number of APs/BSs
cfg.sys.U  = 8;        % number of users
cfg.sys.Nt = 4;        % Tx antennas per AP; use 8 for paper-like setting
cfg.sys.Nr = 4;        % Rx antennas per AP; use 8 for paper-like setting

%% Geometry. Units: meter.
cfg.area.ap_x = [-70, 70];
cfg.area.ap_y = [-100, 100];
cfg.area.user_x = [-100, 100];
cfg.area.user_y = [-150, 150];
cfg.area.target_x = [-20, 20];
cfg.area.target_y = [-20, 20];
cfg.area.min_ap_distance = 10;
cfg.area.min_user_ap_distance = 10;
cfg.area.min_target_ap_distance = 10;

% If true and M=4, uses a paper-like AP geometry. Users and target remain random.
cfg.area.use_paper_like_ap_layout = false;

%% Carrier and array
cfg.radio.fc = 4.2e9;
cfg.radio.c = 3e8;
cfg.radio.lambda = cfg.radio.c / cfg.radio.fc;
cfg.radio.antenna_spacing = cfg.radio.lambda / 2;
cfg.radio.bandwidth = 100e6;
cfg.radio.frame_duration = 1e-3;  % T in CRLB expression

% Squared effective bandwidth beta_u^2.
% For a flat spectrum over bandwidth B, E[f^2] = B^2/12.
cfg.radio.effective_bandwidth_sq = (cfg.radio.bandwidth^2)/12 * ones(cfg.sys.U,1);

%% Communication channel
cfg.comm.noise_dBm = -80;
cfg.comm.noise_power = dbm2w(cfg.comm.noise_dBm);
cfg.comm.pathloss_ref_dB = -30;  % large-scale power gain at d_ref
cfg.comm.pathloss_ref_dist = 1;
cfg.comm.pathloss_exp = 3;
cfg.comm.rician_K_dB = 3;
cfg.comm.rician_K = db2pow(cfg.comm.rician_K_dB);
cfg.comm.Gamma_dB = 10;
cfg.comm.Gamma = db2pow(cfg.comm.Gamma_dB);

%% Sensing channel / CRLB
cfg.sensing.noise_dBm = -80;
cfg.sensing.N0 = dbm2w(cfg.sensing.noise_dBm);  % used as N0-like scalar in implementation
cfg.sensing.kappa_ref_dB = -30;
cfg.sensing.kappa_ref_amp = 10^(cfg.sensing.kappa_ref_dB/20);
cfg.sensing.d_ref = 1;
cfg.sensing.rcs_amp = 1;
cfg.sensing.random_rcs_phase = true;

%% Power constraint
cfg.power.Pmax_dBm = 30;
cfg.power.Pmax = dbm2w(cfg.power.Pmax_dBm);

%% SCA / SDR solver
cfg.opt.max_iter = 50;
cfg.opt.tol = 1e-5;
cfg.opt.penalty_eta = 1e-2;
% Avoid arbitrary dominant-eigenvector bias at W_u^(0)=alpha I.
% In the first SCA step, the rank-one penalty is disabled; afterwards it is
% applied using the dominant eigenvector of the previous iterate.
cfg.opt.first_iter_no_rank_penalty = true;
cfg.opt.tie_aware_rank_penalty = true;
cfg.opt.eig_tie_tol = 1e-8;
cfg.opt.eps_pos = 1e-8;
cfg.opt.cvx_quiet = true;
cfg.opt.cvx_solver = 'sdpt3';        % '', 'sdpt3', or 'mosek'
cfg.opt.init_power_scale = 0.8;  % initial W trace uses this fraction of total AP power
cfg.opt.verbose = true;

% Rank-continuation SDR-SCA. The first SDP computes the unpenalized lower
% bound. Subsequent equivalent block-cone SDPs reduce a linearized rank
% leakage budget while continuing to minimize the CRLB surrogate itself.
cfg.sdr.max_continuation_iter = 20;
cfg.sdr.rank_budget_shrink = 0.25;
cfg.sdr.max_budget_backtracking = 4;
cfg.sdr.budget_relax_fraction = 0.5;
cfg.sdr.normalized_rank_tolerance = 1e-8;
cfg.sdr.incumbent_rank_tolerance = 1e-6;
cfg.sdr.min_rank_refine_iter = 2;
cfg.sdr.max_rank_refine_iter = 5;
cfg.sdr.relative_crlb_tolerance = cfg.opt.tol;
cfg.sdr.relative_rank_tolerance = 1e-3;
cfg.sdr.normalized_positive_floor = 1e-9;
cfg.sdr.feasibility_tol = 2e-4;
cfg.sdr.cvx_quiet = cfg.opt.cvx_quiet;
cfg.sdr.cvx_solver = cfg.opt.cvx_solver;
cfg.sdr.verbose = cfg.opt.verbose;
cfg.sdr.algorithm_version = 'blockwise_rank_budget_continuation_v1';

%% Fixed-direction ZF-BF / SO-BF power allocation
% Equations (35)-(36) use F_{m,u}=sqrt(P_{m,u}) as optimization variables.
% The implementation first solves the SINR/per-AP-power SOCP and then uses
% feasible proximal SCA/SOCP steps to reduce the physical sensing CRLB.
cfg.baseline.power_alloc.max_sca_iter = 8;
cfg.baseline.power_alloc.max_backtracking = 6;
cfg.baseline.power_alloc.initial_radius = 0.25;
cfg.baseline.power_alloc.max_radius = 1.0;
cfg.baseline.power_alloc.min_radius = 1e-5;
cfg.baseline.power_alloc.radius_shrink = 0.5;
cfg.baseline.power_alloc.radius_expand = 1.5;
cfg.baseline.power_alloc.gradient_tol = 1e-8;
cfg.baseline.power_alloc.objective_tol = 1e-8;
cfg.baseline.power_alloc.sinr_margin = 1e-3;
cfg.baseline.power_alloc.power_margin = 1e-6;
cfg.baseline.power_alloc.check_tol = 1e-4;
cfg.baseline.power_alloc.minimum_gamma_dB = -40;
cfg.baseline.power_alloc.gamma_search_span_dB = 60;
cfg.baseline.power_alloc.gamma_bisection_iters = 10;
cfg.baseline.power_alloc.infeasible_gamma_backoff = 2e-3;

%% MRFP-SCA: multi-ratio / matrix-fractional SCA baseline
% The CRLB ratio is written as
%   tr(G^{-1}) = 1/(g_a-g_ab^2/g_b) + 1/(g_b-g_ab^2/g_a).
% Each iteration constructs a Loewner-affine lower model of the 2-by-2
% Fisher matrix and solves the resulting vector-domain SOCP.  No lifted
% N-by-N transmit covariance variable is introduced.
cfg.mrfp.max_iter = 30;
cfg.mrfp.min_iter = 2;
cfg.mrfp.crlb_rel_tol = 1e-5;
cfg.mrfp.fim_floor = 1e-8;
cfg.mrfp.prox_weight = 1e-4;
cfg.mrfp.accept_rel_tol = 1e-8;
cfg.mrfp.max_backtracking = 12;
cfg.mrfp.feasibility_tol = 2e-4;
cfg.mrfp.cvx_quiet = true;
cfg.mrfp.cvx_solver = cfg.opt.cvx_solver;
cfg.mrfp.verbose = true;
cfg.mrfp.print_every = 1;



%% First-order ALM baseline without saddle escaping
cfg.alm.max_outer = 50;
cfg.alm.max_inner = 100;
cfg.alm.min_outer = 3;
cfg.alm.stationarity_tol = 1e-3;
cfg.alm.outer_crlb_rel_tol = 1e-5;
cfg.alm.stop_on_feasible_only = false;
cfg.alm.rho0 = 1;
cfg.alm.rho_growth = 3;
cfg.alm.rho_max = 1e10;
cfg.alm.rho_update_threshold = 0.9;
cfg.alm.feas_tol = 1e-4;
cfg.alm.inner_grad_tol = 1e-6;
cfg.alm.inner_stepsize = 1e-3;
cfg.alm.min_stepsize = 1e-12;
cfg.alm.max_backtracking = 25;
cfg.alm.backtracking_beta = 0.5;
cfg.alm.armijo_c = 1e-4;
cfg.alm.grad_clip = 1e3;
% The logarithmic ratio has the same positive-domain minimizers and avoids
% network-size-dependent premature stopping under fixed KKT tolerances.
cfg.alm.objective_type = 'log_crlb_ratio';
cfg.alm.objective_scale_mode = 'none';
cfg.alm.init_mode = 'zf';              % shared global aggregate-ZF entry
cfg.alm.identity_power_scale = cfg.opt.init_power_scale;
cfg.alm.real_domain = true;              % FO-ALM descent is carried out on [Re(w); Im(w)]
cfg.alm.det_floor = 1e-30;
cfg.alm.use_power_safeguard = false;
cfg.alm.power_safeguard_factor = 5;
cfg.alm.verbose = true;
cfg.alm.print_every = 20;

%% Perturbed projected first-order ALM (PPFO-ALM)
% Practical constrained adaptation of Jin et al., ICML 2017:
% https://proceedings.mlr.press/v70/jin17a.html
% The base PFO solver leaves this disabled; the PPFO wrapper enables it.
cfg.ppfo.enabled = false;
cfg.ppfo.seed = cfg.seed+5000;
cfg.ppfo.pg_trigger_tol = 0.88;
cfg.ppfo.perturb_radius_relative = 2e-3;
cfg.ppfo.min_perturb_radius = 1e-8;
cfg.ppfo.power_active_tol = 1e-6;
cfg.ppfo.escape_window_iters = 20;
cfg.ppfo.escape_decrease_relative = 1e-5;
cfg.ppfo.escape_requires_crlb_improvement = true;
cfg.ppfo.escape_crlb_decrease_relative = 0;
cfg.ppfo.max_perturbations = 4;
cfg.ppfo.min_outer = 3;
cfg.ppfo.max_sinr_violation = 1e-2;
cfg.ppfo.outer_patience = 4;
cfg.ppfo.incumbent_feas_tol = 1e-4;


%% Trust-region ALM baseline with projected AP power constraints
% Safeguarded PTR-ALM. Only TR settings are changed here; the SCA
% and first-order ALM settings above are intentionally left untouched.
% Use the same maximum iteration budget as PFO-/PPFO-ALM. Paper-facing
% PTR-ALM synchronizes from cfg.alm in configure_costar_tr_stage, so cfg.alm
% is the single ALM-baseline budget source. COSTAR does not use this budget.
cfg.tr.max_outer = cfg.alm.max_outer;
cfg.tr.max_inner = cfg.alm.max_inner;
cfg.tr.min_outer = 3;
cfg.tr.min_inner = 3;
cfg.tr.rho0 = 1;
cfg.tr.rho_growth = 3;
cfg.tr.rho_max = 1e10;
cfg.tr.rho_update_threshold = 0.9;
cfg.tr.feas_tol = 1e-4;

% Conservative trust-region radii. Large radii are problematic because the
% per-AP power projection can strongly distort a tCG step at the active
% power boundary.
cfg.tr.delta0 = 1e-3;
cfg.tr.delta_max = 5e-2;
cfg.tr.delta_min = 1e-8;
cfg.tr.eta_accept = 0.10;
cfg.tr.eta_shrink = 0.25;
cfg.tr.eta_expand = 0.75;
cfg.tr.gamma_shrink = 0.5;
cfg.tr.gamma_expand = 1.5;
cfg.tr.gamma_pg_success = 1.0;
cfg.tr.boundary_tol = 0.9;
cfg.tr.pred_tol = 1e-14;
cfg.tr.min_actual_decrease = 1e-12;
cfg.tr.accept_actual_decrease_without_ratio = true;

% tCG and finite-difference Hessian-vector products.
cfg.tr.cg_max_iter = 60;
cfg.tr.cg_tol = 1e-3;
cfg.tr.cg_nc_tol = 1e-8;       % true negative curvature only if curvature <= -tol
cfg.tr.cg_flat_tol = 1e-12;
cfg.tr.hessian_mode = 'structured';  % 'structured' or 'finite_difference'
cfg.tr.hv_eps = 1e-5;

% Joint tangent-space model for active per-AP power balls and nearly active
% SINR inequalities. This prevents oversized tCG steps from crossing the
% shifted-ALM active-set kink in first order.
cfg.tr.use_tangent_projection = true;
cfg.tr.power_active_tol = 1e-5;
cfg.tr.tangent_grad_floor = 1e-12;
cfg.tr.use_sinr_tangent_projection = true;
cfg.tr.sinr_tangent_active_tol = 1e-3;
cfg.tr.tangent_rank_tol = 1e-10;
cfg.tr.return_best_feasible = true;
cfg.tr.incumbent_feas_tol = 1e-4;

% Projected-gradient fallback candidate. This safeguard prevents the TR
% method from performing much worse than Proj-FO-ALM when the FD Hessian
% model is inaccurate.
cfg.tr.use_pg_fallback = true;
cfg.tr.pg_step = 1e-3;
cfg.tr.pg_alpha0 = 1e-3;
cfg.tr.pg_respect_delta = false;
cfg.tr.pg_max_step_fraction = 1.0;
cfg.tr.pg_max_backtracking = 30;
cfg.tr.pg_backtracking_beta = 0.5;
cfg.tr.pg_armijo_c = 1e-4;
cfg.tr.pg_accept_best_decrease = true;
cfg.tr.grad_clip = 1e3;

% Delay second-order tCG steps until the SINR ALM subproblem is close to
% feasible. Early infeasible iterates are dominated by feasibility
% restoration, where projected-gradient steps are more reliable.
cfg.tr.tcg_start_outer = 3;
cfg.tr.tcg_start_sinr_violation = cfg.tr.feas_tol;

cfg.tr.inner_pg_tol = 1e-4;
cfg.tr.inner_plateau_pg_tol = 1e-2;
cfg.tr.outer_pg_tol = 1e-4;
cfg.tr.inner_crlb_rel_tol = 1e-6;
cfg.tr.outer_crlb_rel_tol = 1e-5;
cfg.tr.inner_patience = 6;
cfg.tr.outer_patience = 2;
cfg.tr.step_tol = 1e-10;
cfg.tr.neg_curv_tol = 1e-5;
cfg.tr.enable_curvature_diagnostic = false;
cfg.tr.lanczos_max_iter = 20;
cfg.tr.lanczos_tol = 1e-8;
cfg.tr.require_stationarity_for_stop = false;
cfg.tr.stop_inner_on_feasible = false;
cfg.tr.stop_on_feasible_only = false;
cfg.tr.clip_feasibility_crossing = false;
cfg.tr.feasibility_crossing_bisection_iters = 12;
cfg.tr.objective_type = 'log_crlb_ratio';  % 'crlb_ratio' or 'log_crlb_ratio'
cfg.tr.objective_scale_mode = 'none';
cfg.tr.det_floor = 1e-30;
cfg.tr.init_mode = 'zf';               % same entry as PFO-/PPFO-ALM
cfg.tr.identity_power_scale = cfg.opt.init_power_scale;
cfg.tr.verbose = true;
cfg.tr.print_every = 5;

% Safeguarded PTR baseline: PFO globalization followed by a bounded PTR
% refinement. Both stages are included in the reported wall-clock time.
cfg.tr.pfo_warmstart_enabled = true;
cfg.tr.refine_max_outer = 8;
cfg.tr.refine_max_inner = 50;
cfg.tr.refine_delta0 = 1e-2;
cfg.tr.refine_delta_max = 1e-1;
cfg.tr.refine_gamma_expand = 2;
cfg.tr.refine_gamma_pg_success = 1.25;

%% COSTAR: shared feasible initialization and a safeguarded exact-HVP loop.
% All five iterative comparators receive one shared stacked-ZF point.
% initialize_shared_stacked_zf performs no CRLB descent and no CVX solve.
% initializer_mode is compatibility metadata, not a solve_costar dispatch.
% Legacy initializer parameters below are not used by the figure programs.
cfg.costar.initializer_mode = 'global_aggregate_zf';
cfg.costar.print_initializer = false;
cfg.costar.phase1_tcg_start_outer = 2;
cfg.costar.phase1_verbose = false;
cfg.costar.rtr_max_iter = 200;
cfg.costar.rtr_seed_offset = 4000;
cfg.costar.phase2_variant = 'core_kkt_triggered';
cfg.costar.phase2_version = 'shared_zf_accuracy_controlled_sosp_v1';
cfg.costar.sosp = struct('max_iter',120,'epsilon_KKT',1e-5, ...
    'epsilon_H',1e-5,'failure_probability',1e-3,'meo_max_iter',Inf, ...
    'use_structured_psd',true,'curvature_space','ambient', ...
    'hessian_model','exact_lagrangian','enable_tcg_negative_curvature',true, ...
    'sinr_active_tol',1e-4,'power_active_tol',1e-5);

% Core AS-RTR keeps the constrained second-order model while removing
% repeated HVP/trial evaluations. Lanczos is triggered only at a genuine
% first-order stopping candidate; tCG detects negative curvature otherwise.
cfg.costar.core_rtr.curvature_policy = 'terminal';
cfg.costar.core_rtr.reuse_tcg_model = true;
cfg.costar.core_rtr.fuse_trial_evaluation = true;
cfg.costar.core_rtr.enable_sinr_correction = true;
cfg.costar.core_rtr.multiplier_mode = 'nnls';
cfg.costar.core_rtr.cg_max_iter = 30;
cfg.costar.core_rtr.sinr_retraction_max_iter = 4;
cfg.costar.core_rtr.lanczos_max_iter = 20;
cfg.costar.core_rtr.hessian_model = 'exact_lagrangian';
cfg.costar.core_rtr.enable_tcg_negative_curvature = true;
% The explicit Ritz trial was accepted zero times in the 30-instance audit.
% The fast-stage finite Lanczos test is diagnostic only. The subsequent
% accuracy-controlled stage supplies the conditional curvature enclosure.
cfg.costar.core_rtr.enable_ritz_escape = false;
% The central-difference settings below are retained only for the
% hessian_model='lagrangian_fd' ablation.
cfg.costar.core_rtr.fd_trust_scaling = true;
cfg.costar.core_rtr.fd_trust_fraction = 0.1;
cfg.costar.core_rtr.fd_min_displacement = 1e-8;

% Low-complexity theorem-aware initializer parameters.
cfg.costar.low_complexity_init.sinr_margin_fraction = 1e-3;
cfg.costar.low_complexity_init.power_margin_fraction = 1e-3;
cfg.costar.low_complexity_init.max_refine_iter = 20;
cfg.costar.low_complexity_init.max_backtracking = 8;
cfg.costar.low_complexity_init.initial_radius = 0.25;
cfg.costar.low_complexity_init.max_radius = 1;
cfg.costar.low_complexity_init.min_radius = 1e-6;
cfg.costar.low_complexity_init.radius_shrink = 0.5;
cfg.costar.low_complexity_init.radius_expand = 1.5;
cfg.costar.low_complexity_init.gradient_tolerance = 1e-8;
cfg.costar.low_complexity_init.objective_tolerance = 1e-9;
cfg.costar.low_complexity_init.projection_cycles = 60;
cfg.costar.low_complexity_init.feasibility_tolerance = 1e-4;
cfg.costar.low_complexity_init.evaluate_normal_certificate = true;
cfg.costar.low_complexity_init.certificate_delta0 = 5e-2;
cfg.costar.low_complexity_init.certificate_feas_tol = 1e-6;
cfg.costar.low_complexity_init.certificate_sinr_retraction_margin = 1e-5;
cfg.costar.low_complexity_init.certificate_max_normal_steps = 4;

% Global aggregate theorem-aware ZF initializer parameters.
cfg.costar.global_zf_init.sinr_margin_fraction = 1e-3;
cfg.costar.global_zf_init.power_margin_fraction = 1e-3;
cfg.costar.global_zf_init.rank_tolerance = 1e-11;
cfg.costar.global_zf_init.max_projection_cycles = 300;
cfg.costar.global_zf_init.projection_tolerance = 1e-9;
cfg.costar.global_zf_init.max_refine_iter = 0;
cfg.costar.global_zf_init.max_backtracking = 8;
cfg.costar.global_zf_init.initial_radius = 0.05;
cfg.costar.global_zf_init.max_radius = 0.5;
cfg.costar.global_zf_init.min_radius = 1e-7;
cfg.costar.global_zf_init.radius_shrink = 0.5;
cfg.costar.global_zf_init.radius_expand = 1.5;
cfg.costar.global_zf_init.gradient_tolerance = 1e-8;
cfg.costar.global_zf_init.objective_tolerance = 1e-9;
cfg.costar.global_zf_init.feasibility_tolerance = 1e-4;
cfg.costar.global_zf_init.evaluate_normal_certificate = false;
cfg.costar.global_zf_init.certificate_delta0 = 5e-2;
cfg.costar.global_zf_init.certificate_feas_tol = 1e-6;
cfg.costar.global_zf_init.certificate_sinr_retraction_margin = 1e-5;
cfg.costar.global_zf_init.certificate_max_normal_steps = 4;

%% SINR sweep
cfg.sweep.Gamma_dB_list = 0:5:20;

%% I/O
cfg.io.results_dir = fullfile(pwd, 'results');

%% Paper-aligned COSTAR: one protected loop, deterministic tangent inertia.
% The old core settings above remain available to independent ablations;
% solve_costar does not dispatch the old fast/terminal two-stage pipeline.
cfg.costar.paper = costar_paper_options(struct());
cfg.costar.sosp = cfg.costar.paper; % simulation/plot compatibility alias
cfg.costar.rtr_max_iter = cfg.costar.paper.max_iter;
cfg.costar.phase2_variant = 'single_safeguarded_paper_loop';
cfg.costar.phase2_version = 'safeguard_adaptive_20261019';
end
