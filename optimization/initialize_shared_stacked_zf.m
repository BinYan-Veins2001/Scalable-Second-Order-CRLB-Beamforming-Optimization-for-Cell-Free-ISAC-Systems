function entry=initialize_shared_stacked_zf(cfg,channels,smat)
%INITIALIZE_SHARED_STACKED_ZF Independent common entry; no CRLB descent.
options=cfg.costar.global_zf_init;
options.max_refine_iter=0;
options.evaluate_normal_certificate=false;
entry=struct('success',false,'wcell',[],'Wmat',[],'z',[], ...
    'runtime_s',NaN,'info',struct(),'message','', ...
    'method','Stacked global ZF','version','shared_stacked_zf_no_descent_v1');
timer=tic;
try
    [entry.wcell,entry.info]=build_global_aggregate_zf_initializer( ...
        cfg,channels,smat,options);
    entry.Wmat=vectorize_wcell_to_matrix(cfg,entry.wcell);
    entry.z=complex_matrix_to_real_vector(entry.Wmat);
    entry.metrics=evaluate_beamforming_solution(cfg,channels,smat,entry.wcell);
    entry.success=logical(entry.info.entry_feasible) && ...
        isfinite(entry.metrics.crlb) && entry.metrics.crlb>0;
    if ~entry.success
        entry.message='Stacked-ZF feasibility or positive-FIM check failed.';
    end
catch exception
    entry.message=exception.message;
    entry.error_identifier=exception.identifier;
end
entry.runtime_s=toc(timer);
end
