function run = simtwc_run_method(key,p,opts)
%SIMTWC_RUN_METHOD Uniform end-to-end adapter for all seven methods.

if nargin<3; opts=struct(); end
if ~isfield(opts,'capture_history'); opts.capture_history=true; end
if ~isfield(opts,'rtr_opts'); opts.rtr_opts=struct(); end
if ~isfield(opts.rtr_opts,'max_iter')
    opts.rtr_opts.max_iter=p.cfg.costar.rtr_max_iter;
end
if ~isfield(opts.rtr_opts,'seed'); opts.rtr_opts.seed=p.seed+4000; end

sol=struct(); hist=struct(); entry_w=[]; entry_time=0;
iterative=any(strcmp(key,{'pfo_alm','ptr_alm','mrfp_sca','costar','sdr_sca'}));
if iterative
    if ~isfield(p,'initialization') || ~p.initialization.success
        error('Simulation:SharedZFFailed', ...
            'No valid shared stacked-ZF entry for seed %d.',p.seed);
    end
    entry_w=p.initialization.wcell;
    entry_time=p.initialization.runtime_s;
end
rng(p.seed,'twister');
timer=tic;
switch key
    case 'zf_bf'
        [wcell,info]=build_zf_beamformer(p.cfg,p.channels,p.smat);
        sol.wcell=wcell; sol.info=info; sol.success=true;
    case 'so_bf'
        [wcell,info]=build_sensing_oriented_beamformer( ...
            p.cfg,p.scenario,p.channels,p.smat);
        sol.wcell=wcell; sol.info=info; sol.success=true;
    case {'pfo_alm','ptr_alm'}
        if strcmp(key,'pfo_alm')
            [sol,hist]=solve_alm_first_order_projected_power( ...
                p.cfg,p.channels,p.smat,entry_w);
        else
            [sol,hist]=solve_tr_alm(p.cfg,p.channels,p.smat,entry_w);
        end
    case 'mrfp_sca'
        mopts=struct('max_iter',p.cfg.mrfp.max_iter);
        [sol,hist]=solve_mrfp_sca( ...
            p.cfg,p.channels,p.smat,entry_w,mopts);
    case 'costar'
        [sol,hist]=solve_costar( ...
            p.cfg,p.channels,p.smat,entry_w,opts.rtr_opts);
    case 'sdr_sca'
        [lifted,hist]=solve_sca_sdr( ...
            p.cfg,p.channels,p.smat,entry_w);
        if isempty(lifted.Wcell) || ~logical(lifted.success)
            error('SDR-SCA returned no certified solution (%s).',lifted.cvx_status);
        end
        sol=lifted;
        sol.wcell=recover_beamformers_from_W(p.cfg,lifted.Wcell);
    otherwise
        error('Unknown method key: %s',key);
end
optimization_time=toc(timer);
elapsed=entry_time+optimization_time;

if ~isfield(sol,'wcell') || isempty(sol.wcell)
    if isfield(sol,'Wmat') && ~isempty(sol.Wmat)
        sol.wcell=matrix_to_wcell(p.cfg,sol.Wmat);
    else
        error('%s did not return a rank-one beamformer.',key);
    end
end
Wmat=vectorize_wcell_to_matrix(p.cfg,sol.wcell);
metrics=evaluate_beamforming_solution(p.cfg,p.channels,p.smat,sol.wcell);
kkt=simtwc_unified_kkt(p,Wmat);
s=simtwc_config();
feasible=isfinite(metrics.crlb) && metrics.crlb>0 && ...
    all(metrics.sinr>=p.cfg.comm.Gamma*(1-s.feasibility_tol)) && ...
    all(metrics.ap_power<=p.cfg.power.Pmax*(1+s.feasibility_tol));

run=struct();
run.key=key; run.seed=p.seed; run.elapsed_s=elapsed;
run.initialization_runtime_s=entry_time;
run.optimization_runtime_s=optimization_time;
run.runtime_definition='shared stacked-ZF initialization + optimization';
run.algorithm_version=field_or(sol,'algorithm_version','baseline_20261003');
if iterative
    run.initial_Wmat=p.initialization.Wmat;
    run.initializer_info=p.initialization.info;
    run.initializer_version=p.initialization.version;
    run.initialization_success=p.initialization.success;
end
run.conditional_sosp_reached=field_or(sol,'conditional_sosp_reached',false);
run.metrics=metrics; run.feasible=feasible; run.kkt=kkt;
run.Wmat=Wmat; run.wcell=sol.wcell;
run.solver_success=field_or(sol,'success',true);
run.message=string(field_or(sol,'message','completed'));
if opts.capture_history
    run.history=simtwc_canonical_history( ...
        key,p,sol,hist,Wmat,elapsed,entry_w,entry_time);
else
    run.history=struct();
end
run.solution=strip_solution(sol);
end

function sol=strip_solution(sol)
for name={'wcell','Wmat','Wcell','z','hist','phase1_solution'}
    if isfield(sol,name{1}); sol=rmfield(sol,name{1}); end
end
end

function value=field_or(s,name,fallback)
if isfield(s,name); value=s.(name); else; value=fallback; end
end
