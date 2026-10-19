function s = simtwc_config()
%SIMTWC_CONFIG Frozen specification for the TWC simulation figures.

s.root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
s.data_dir = fullfile(s.root,'data');
s.figure_dir = fullfile(s.root,'figures');
s.schema_version = 'costar_twc_shared_zf_sosp_20261003_v2';
s.code_version = 'SharedZF_PTR_ExactHVP_AccuracyControlled_COSTAR';
s.smoke = strcmpi(strtrim(getenv('COSTAR_SMOKE')),'1');
s.quick = strcmpi(strtrim(getenv('COSTAR_QUICK')),'1') || s.smoke;
s.resume = ~strcmpi(strtrim(getenv('COSTAR_NO_RESUME')),'1');
s.feasibility_tol = 1e-4;
s.kkt_tol = 1e-3;
s.kkt_active_tol = 1e-3;
s.kkt_sinr_active_tol = 1e-4;
s.kkt_power_active_tol = 1e-5;
s.initialization_policy = 'external_shared_stacked_zf_without_crlb_descent';
s.runtime_definition = 'shared initialization + optimizer; diagnostics excluded';
defaults=default_config();
s.costar_options=defaults.costar.sosp;
s.timing_repetitions = 1;
s.seeds = 2026:2035;
s.methods = {'zf_bf','so_bf','pfo_alm','ptr_alm', ...
    'mrfp_sca','costar','sdr_sca'};
s.iterative_methods = {'pfo_alm','ptr_alm','mrfp_sca','costar','sdr_sca'};
s.default = struct('M',16,'U',8,'Nt',4,'Nr',4, ...
    'gamma_dB',10,'min_ap_distance',10);

s.fig01 = s.default;
s.fig01.seed = 2026;
s.fig01.grid_axis = linspace(-20,20,101);
s.fig01.map_methods = {'sdr_sca','costar'};

s.fig02 = s.default;
s.fig02.seed = 2026;

s.fig03 = s.default;
s.fig03.seeds = s.seeds;

s.fig04 = s.fig03;
s.fig04.crlb_factor = 1.01;
s.fig04.kkt_tol = 1e-3;

s.fig05 = s.default;
s.fig05.seed = 2026;

s.fig06.seeds = s.seeds;
s.fig06.regimes = struct( ...
    'name',{'Default','High SINR','More APs'}, ...
    'M',{16,16,32},'U',{8,8,8},'Nt',{4,4,4},'Nr',{4,4,4}, ...
    'gamma_dB',{10,20,10},'min_ap_distance',{10,10,10});

s.fig07 = s.default;
s.fig07.seeds = s.seeds;
s.fig07.gamma_dB = [5 10 15 20 25];

s.fig08.seeds_compact = s.seeds;
s.fig08.compact = struct('M_values',16:4:40,'U',8,'Nt',4,'Nr',4, ...
    'gamma_dB',10,'min_ap_distance',10);
% SDR-SCA is evaluated only on the sparse grid below and only while
% M <= sdr_sca_max_M.  This is an inclusive, user-editable memory/runtime
% safeguard: use 16 for only M=16, 28 for M={16,28}, or 40 to restore all
% three original SDR-SCA points.  Other algorithms are unaffected.
s.fig08.sdr_sca_max_M = 28;
s.fig08.compact_sdr_M = [16 28 40];
s.fig08.seeds_large = 2001:2010;
s.fig08.large = struct('M_values',20:5:40,'U',16,'Nt',8,'Nr',8, ...
    'gamma_dB',10,'min_ap_distance',10,'area',[-500 500]);

s.fig09 = s.default;
s.fig09.seeds = s.seeds;
s.fig09.power_dBm = [20 25 30 35 40];
s.fig09.methods = {'costar','mrfp_sca','sdr_sca'};

s.fig10 = s.default;
s.fig10.seeds = s.seeds;
s.fig10.tau = [1e-4 1e-3 1e-2 1e-1];

% Fig. 12: common endpoint second-order audit.  The beamformers are reused
% from Fig. 3; only the quotient-space finite-Lanczos diagnostics are new.
s.fig12 = s.default;
s.fig12.source_figure = 3;
s.fig12.seeds = s.seeds;
s.fig12.methods = s.methods;
s.fig12.lanczos_max_iter = 40;
s.fig12.lanczos_tol = 1e-9;
s.fig12.sosp_tol = 1e-3;
s.fig12.sinr_active_tol = s.kkt_sinr_active_tol;
s.fig12.power_active_tol = s.kkt_power_active_tol;

% Fig. 13: controlled local one-step test for the conditional quadratic
% convergence statement.  A high-accuracy COSTAR trajectory supplies
% feasible points near its final phase orbit; every selected point is then
% advanced by one accurately solved, nominally interior AS-RTR step.
s.fig13 = s.default;
s.fig13.seed = 2026;
s.fig13.reference_max_iter = 300;
s.fig13.reference_stationarity_tol = 1e-8;
s.fig13.reference_curvature_tol = 1e-7;
s.fig13.reference_cg_tol = 1e-10;
s.fig13.reference_cg_max_iter = 100;
s.fig13.reference_lanczos_max_iter = 40;
s.fig13.local_points = 10;
s.fig13.local_error_max = 2e-1;
s.fig13.local_error_floor = 1e-10;
s.fig13.one_step_delta = 5e-1;
s.fig13.sinr_active_tol = 1e-3;
s.fig13.power_active_tol = 1e-5;

% Fig. 14: separate the variable-dimension effects of the number of users
% and antennas.  SDR-SCA is intentionally excluded from this scaling test.
% Edit these three fields to reduce or enlarge the publication experiment.
s.fig14 = struct();
s.fig14.base = s.default;
s.fig14.seeds = s.seeds;
s.fig14.methods = {'zf_bf','so_bf','pfo_alm','ptr_alm', ...
    'mrfp_sca','costar'};
s.fig14.U_values = [4 8 12 16];
s.fig14.Nt_values = [2 4 6 8];

% Fig. 15: unified 2-by-3 network-dimension scaling figure.  The top row
% reports CRLB and the bottom row reports end-to-end runtime.  Each column
% changes only one declared dimension; all remaining parameters are fixed.
s.fig15 = struct();
s.fig15.seeds = s.seeds;       % Monte Carlo seeds; edit to change trial count.
s.fig15.methods = s.methods;   % All seven paper algorithms.
s.fig15.gamma_dB = 10;         % Common target SINR for all three sweeps.
s.fig15.min_ap_distance = 10;  % Common minimum AP separation (m).

% Column 1: AP-number sweep.  SDR-SCA is launched only for M values not
% exceeding the inclusive, user-editable cap sdr_sca_max_M.  For example,
% use 20 for M={16,20}, or Inf to remove the resource safeguard.
s.fig15.ap.M_values = 16:4:40;
s.fig15.ap.U = 8;
s.fig15.ap.Nt = 4;
s.fig15.ap.Nr = 4;
s.fig15.sdr_sca_max_M = 20;

% Column 2: transmit-antenna sweep.  M=16 and U=8 are fixed.  By default
% only Nt changes and Nr stays at 4.  Set couple_receive_antennas=true if
% the desired experiment is Nt=Nr at every point.
s.fig15.antenna.Nt_values = [4 8 12 16];
s.fig15.antenna.M = 16;
s.fig15.antenna.U = 8;
s.fig15.antenna.Nr = 4;
s.fig15.antenna.couple_receive_antennas = false;

% Column 3: user-number sweep.  M=16 and Nt=Nr=4 are fixed.
s.fig15.user.U_values = [8 12 14 16];
s.fig15.user.M = 16;
s.fig15.user.Nt = 4;
s.fig15.user.Nr = 4;

% Fig. 16: one-realization COSTAR convergence and phase-reduced curvature.
% Curvature is audited from every stored RTR iterate with the Fig. 12 rule.
s.fig16 = s.default;
s.fig16.seed = 2026;

% Fig. 17: endpoint multiplier/constraint audit.  Increase mc_trials to
% aggregate consecutive seeds beginning at base_seed; the default is five.
s.fig17 = s.default;
s.fig17.base_seed = 2026;
s.fig17.mc_trials = 5;
s.fig17.seeds = s.fig17.base_seed + (0:s.fig17.mc_trials-1);

% Fig. 18: numerical diagnostic for the conditional local quadratic-rate
% theorem. A high-accuracy endpoint is reused from Fig. 13, while local
% feasible-side perturbations test the one-step relation e_{k+1}=O(e_k^2).
s.fig18 = s.fig13;
s.fig18.seeds = s.fig18.seed;
fig18_seed_text = strtrim(getenv('COSTAR_FIG18_SEEDS'));
if ~isempty(fig18_seed_text)
    parsed_seeds = sscanf(strrep(fig18_seed_text,',',' '),'%f').';
    if ~isempty(parsed_seeds); s.fig18.seeds = parsed_seeds; end
end
% The window is below the observed local-regime boundary and above the
% double-precision/reference-solve error floor for the default realization.
s.fig18.local_test_errors = ...
    logspace(log10(5e-6),log10(2.5e-5),11);
s.fig18.reference_max_iter = 300;
s.fig18.reference_stationarity_tol = 1e-8;
s.fig18.reference_cg_tol = 1e-10;
s.fig18.reference_cg_max_iter = 100;
s.fig18.local_cg_tol = 1e-10;
s.fig18.local_cg_max_iter = 100;
s.fig18.local_forcing_constant = 0.1;
s.fig18.local_cg_absolute_floor = 1e-14;
s.fig18.local_manifold_margin = 1e-10;
s.fig18.local_manifold_corrections = 12;
s.fig18.local_manifold_tolerance = 1e-11;
s.fig18.local_manifold_regularization = 1e-12;
s.fig18.local_correction_cap_fraction = 0.5;
s.fig18.local_correction_backtracks = 12;
s.fig18.hvp_fd_steps = logspace(-2,-7,11);

% Fig. 19: paired default-system endpoint curvature comparison. Reuse
% the Fig. 12 diagnostics (which in turn reuse Fig. 3 beamformers).
s.fig19 = s.default;
s.fig19.source_figure = 12;
s.fig19.seeds = s.fig12.seeds;
s.fig19.methods = s.methods;

if s.quick
    s.data_dir = fullfile(s.root,'data','quick');
    s.costar_options.max_iter = 8;
    s.costar_options.meo_max_iter = 24;
    s.figure_dir = fullfile(s.root,'figures','quick');
    s.seeds = s.seeds(1);
    s.fig01.grid_axis = linspace(-20,20,21);
    s.fig03.seeds = s.fig03.seeds(1);
    s.fig04.seeds = s.fig04.seeds(1);
    s.fig06.seeds = s.fig06.seeds(1);
    s.fig07.seeds = s.fig07.seeds(1);
    s.fig07.gamma_dB = s.fig07.gamma_dB([1 end]);
    s.fig08.seeds_compact = s.fig08.seeds_compact(1);
    s.fig08.compact.M_values = s.fig08.compact.M_values([1 4 end]);
    s.fig08.seeds_large = s.fig08.seeds_large(1);
    s.fig08.large.M_values = s.fig08.large.M_values([1 end]);
    s.fig09.seeds = s.fig09.seeds(1);
    s.fig09.power_dBm = s.fig09.power_dBm([1 end]);
    s.fig10.seeds = s.fig10.seeds(1);
    s.fig10.tau = s.fig10.tau([1 end]);
    s.fig12.seeds = s.fig12.seeds(1);
    s.fig12.lanczos_max_iter = 6;
    s.fig13.reference_max_iter = 8;
    s.fig13.reference_cg_max_iter = 12;
    s.fig13.reference_lanczos_max_iter = 6;
    s.fig13.local_points = 3;
    s.fig14.seeds = s.fig14.seeds(1);
    s.fig14.U_values = s.fig14.U_values([1 end]);
    s.fig14.Nt_values = s.fig14.Nt_values([1 end]);
    s.fig15.seeds = s.fig15.seeds(1);
    % Retain the two SDR-SCA AP points and one deliberately capped point.
    s.fig15.ap.M_values = [16 20 40];
    % Use two modest points in the remaining quick-test sweeps so the
    % interface test does not itself become a large SDR-SCA experiment.
    s.fig15.antenna.Nt_values = s.fig15.antenna.Nt_values(1:2);
    s.fig15.user.U_values = s.fig15.user.U_values(1:2);
    s.fig17.mc_trials = 1;
    s.fig17.seeds = s.fig17.base_seed;
    s.fig18 = s.fig13;
    s.fig18.seeds = s.seeds(1);
    s.fig18.local_test_errors = logspace(log10(7e-6),log10(1.75e-5),3);
    s.fig18.local_cg_tol = 1e-10;
    s.fig18.local_cg_max_iter = 40;
    s.fig18.local_forcing_constant = 0.1;
    s.fig18.local_cg_absolute_floor = 1e-14;
    s.fig18.local_manifold_margin = 1e-10;
    s.fig18.local_manifold_corrections = 12;
    s.fig18.local_manifold_tolerance = 1e-11;
    s.fig18.local_manifold_regularization = 1e-12;
    s.fig18.local_correction_cap_fraction = 0.5;
    s.fig18.local_correction_backtracks = 12;
    s.fig18.hvp_fd_steps = logspace(-2,-7,7);
    s.fig19.seeds = s.fig12.seeds;
end
if s.smoke
    % Interface tests only. These results never enter publication archives.
    s.data_dir=fullfile(s.root,'data','smoke');
    s.figure_dir=fullfile(s.root,'figures','smoke');
    s=rewrite_smoke_dimensions(s);
    s.fig08.compact.M_values=[4 6]; s.fig08.compact_sdr_M=[4 6];
    s.fig08.sdr_sca_max_M=6; s.fig08.large.M_values=[4 6];
    s.fig14.U_values=[2 3]; s.fig14.Nt_values=[2 3];
    s.fig15.ap.M_values=[4 6]; s.fig15.sdr_sca_max_M=4;
    s.fig15.antenna.Nt_values=[2 3]; s.fig15.user.U_values=[2 3];
    s.fig06.regimes(3).M=6;
end
end

function v=rewrite_smoke_dimensions(v)
if ~isstruct(v); return; end
for k=1:numel(v)
    names=fieldnames(v(k));
    for j=1:numel(names)
        name=names{j}; value=v(k).(name);
        if strcmp(name,'M'); value=4;
        elseif any(strcmp(name,{'U','Nt','Nr'})); value=2;
        elseif isstruct(value); value=rewrite_smoke_dimensions(value);
        end
        v(k).(name)=value;
    end
end
end
