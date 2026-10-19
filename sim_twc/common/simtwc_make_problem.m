function p = simtwc_make_problem(spec,seed)
%SIMTWC_MAKE_PROBLEM Deterministically construct one simulation realization.

cfg = default_config();
cfg.seed = seed;
cfg.sys.M=spec.M; cfg.sys.U=spec.U;
cfg.sys.Nt=spec.Nt; cfg.sys.Nr=spec.Nr;
cfg.radio.effective_bandwidth_sq = ...
    (cfg.radio.bandwidth^2)/12*ones(cfg.sys.U,1);
cfg.comm.Gamma_dB=spec.gamma_dB;
cfg.comm.Gamma=db2pow(cfg.comm.Gamma_dB);
cfg.area.min_ap_distance=spec.min_ap_distance;
if isfield(spec,'area')
    cfg.area.ap_x=spec.area; cfg.area.ap_y=spec.area;
    cfg.area.user_x=spec.area; cfg.area.user_y=spec.area;
    cfg.area.target_x=0.2*spec.area; cfg.area.target_y=0.2*spec.area;
end
if isfield(spec,'power_dBm')
    cfg.power.Pmax_dBm=spec.power_dBm;
    cfg.power.Pmax=dbm2w(spec.power_dBm);
end
cfg.costar.initializer_mode='global_aggregate_zf';
cfg.costar.global_zf_init.max_refine_iter=0;
cfg.costar.global_zf_init.evaluate_normal_certificate=false;
cfg.costar.core_rtr.hessian_model='exact_lagrangian';
cfg.costar.print_initializer=false;
cfg.opt.verbose=false; cfg.opt.cvx_quiet=true;
cfg.mrfp.verbose=false; cfg.alm.verbose=false; cfg.tr.verbose=false;
cfg.diagnostics.store_iterates=true;
cfg.diagnostics.store_stride=1;

s=simtwc_config();
cfg.costar.sosp=s.costar_options;
if s.quick
    cfg.baseline.power_alloc.max_sca_iter=1;
    cfg.opt.max_iter=1;
    cfg.sdr.max_continuation_iter=0;
    cfg.sdr.normalized_rank_tolerance=1;
    cfg.sdr.incumbent_rank_tolerance=1;
    cfg.alm.max_outer=1; cfg.alm.max_inner=3;
    cfg.tr.max_outer=1; cfg.tr.max_inner=3;
    cfg.tr.refine_max_outer=1; cfg.tr.refine_max_inner=3;
    cfg.mrfp.max_iter=1;
    cfg.costar.rtr_max_iter=3;
end

rng(seed,'twister');
scenario=generate_scenario(cfg);
channels=generate_channels(cfg,scenario);
smat=build_sensing_matrices(cfg,scenario);
p=struct('cfg',cfg,'scenario',scenario,'channels',channels, ...
    'smat',smat,'seed',seed,'spec',spec);
p.initialization=initialize_shared_stacked_zf(cfg,channels,smat);
end
