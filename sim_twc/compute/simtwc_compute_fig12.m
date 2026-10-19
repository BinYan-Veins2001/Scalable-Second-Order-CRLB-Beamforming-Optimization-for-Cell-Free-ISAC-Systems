function archive = simtwc_compute_fig12()
%SIMTWC_COMPUTE_FIG12 Common quotient-space second-order endpoint audit.
%
% Final beamformers are reused from Fig. 3.  No baseline solver is rerun;
% this routine only applies the same active-set, phase-quotient and
% finite-Lanczos diagnostic to every available endpoint.

s=simtwc_setup(); spec=s.fig12;
source=simtwc_load_data(spec.source_figure);
methods=cellstr(string(spec.methods));
data=simtwc_initialize_data(12,spec,methods); path=simtwc_archive_path(12);
nA=numel(methods); nT=numel(spec.seeds);
if ~isfield(data,'runs') || ~isequal(size(data.runs),[nA nT])
    data.runs=cell(nA,nT);
end

for t=1:nT
    seed=spec.seeds(t);
    source_t=find(source.seeds==seed,1);
    if isempty(source_t)
        error('Fig. 3 does not contain seed %d required by Fig. 12.',seed);
    end
    p=simtwc_make_problem(spec,seed);
    for a=1:nA
        if ok(data.runs{a,t}); continue; end
        source_a=find(strcmp(source.methods,methods{a}),1);
        if isempty(source_a)
            data.runs{a,t}=failed(methods{a},seed, ...
                'Method is absent from the Fig. 3 archive.');
            simtwc_save_atomic(path,data);
            continue;
        end
        q=source.runs{source_a,source_t};
        if ~isstruct(q) || ~isfield(q,'status') || ...
                ~strcmp(q.status,'ok') || ~isfield(q,'Wmat') || isempty(q.Wmat)
            data.runs{a,t}=failed(methods{a},seed, ...
                'No usable endpoint is present in the Fig. 3 archive.');
            simtwc_save_atomic(path,data);
            continue;
        end
        fprintf('[Fig12] seed=%d, endpoint SOSP audit: %s\n',seed,methods{a});
        try
            z=complex_matrix_to_real_vector(q.Wmat);
            opts=struct('hessian_model','exact_lagrangian', ...
                'lanczos_max_iter',spec.lanczos_max_iter, ...
                'lanczos_tol',spec.lanczos_tol,'seed',seed+7000+a, ...
                'sinr_active_tol',spec.sinr_active_tol, ...
                'power_active_tol',spec.power_active_tol);
            diagnostic=evaluate_phase2_quotient_second_order_stationarity( ...
                p.cfg,p.channels,p.smat,z,opts);
            record=struct('status','ok','key',methods{a},'seed',seed, ...
                'source_feasible',logical(q.feasible), ...
                'source_crlb',q.metrics.crlb,'diagnostic',diagnostic);
            data.runs{a,t}=record;
        catch ME
            data.runs{a,t}=failed(methods{a},seed,ME.message);
        end
        simtwc_save_atomic(path,data);
    end
end

data.seeds=spec.seeds; data.methods=methods;
[data.kkt_residual,data.normalized_negative_curvature, ...
    data.sosp_residual,data.lambda_min_horizontal,data.feasible]= ...
    extract(data.runs);
data.complete=all(cellfun(@ok,data.runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function [K,N,S,L,F]=extract(runs)
K=nan(size(runs)); N=nan(size(runs)); S=nan(size(runs));
L=nan(size(runs)); F=false(size(runs));
for k=1:numel(runs)
    q=runs{k}; if ~ok(q); continue; end
    d=q.diagnostic;
    K(k)=d.kkt_residual;
    N(k)=d.normalized_negative_curvature;
    S(k)=d.sosp_residual;
    L(k)=d.lambda_min_horizontal;
    F(k)=q.source_feasible;
end
end

function q=failed(key,seed,message)
q=struct('status','failed','key',key,'seed',seed,'message',message);
fprintf('    FAILED %s: %s\n',key,message);
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
