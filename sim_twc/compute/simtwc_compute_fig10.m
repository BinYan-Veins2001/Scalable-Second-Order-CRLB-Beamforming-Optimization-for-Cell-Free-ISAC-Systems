function archive = simtwc_compute_fig10()
%SIMTWC_COMPUTE_FIG10 Active-set threshold sensitivity.

s=simtwc_setup(); spec=s.fig10;
data=simtwc_initialize_data(10,spec,{'costar'}); path=simtwc_archive_path(10);
nQ=numel(spec.tau); nT=numel(spec.seeds);
if ~isfield(data,'runs') || ~isequal(size(data.runs),[nQ,nT])
    data.runs=cell(nQ,nT);
end
for q=1:nQ
    for t=1:nT
        if ok(data.runs{q,t}); continue; end
        p=simtwc_make_problem(spec,spec.seeds(t));
        ropts.rtr_opts=struct('sinr_active_tol',spec.tau(q), ...
            'power_active_tol',spec.tau(q));
        fprintf('[Fig10] tau=%.1e, seed=%d\n',spec.tau(q),spec.seeds(t));
        run=simtwc_try_run('costar',p,true,ropts);
        if ok(run) && isfield(run.history,'raw')
            rtr=run.history.raw.rtr;
            run.active_jitter=sum(rtr.active_change_count);
            run.rtr_iterations=run.solution.work.rtr_iterations;
        else
            run.active_jitter=NaN; run.rtr_iterations=NaN;
        end
        data.runs{q,t}=run;
        simtwc_save_atomic(path,data);
    end
end
data.tau=spec.tau; data.seeds=spec.seeds;
[data.final_crlb,data.kkt_residual,data.active_jitter, ...
    data.rtr_iterations]=extract_results(data.runs);
data.complete=all(cellfun(@ok,data.runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function yes=ok(x)
yes=isstruct(x)&&isfield(x,'status')&&strcmp(x.status,'ok');
end

function [C,K,J,I]=extract_results(runs)
sz=size(runs); C=nan(sz); K=nan(sz); J=nan(sz); I=nan(sz);
for k=1:numel(runs)
    q=runs{k}; if ~ok(q); continue; end
    C(k)=q.metrics.crlb; K(k)=q.kkt.r_kkt;
    if isfield(q,'active_jitter'); J(k)=q.active_jitter; end
    if isfield(q,'rtr_iterations'); I(k)=q.rtr_iterations; end
end
end
