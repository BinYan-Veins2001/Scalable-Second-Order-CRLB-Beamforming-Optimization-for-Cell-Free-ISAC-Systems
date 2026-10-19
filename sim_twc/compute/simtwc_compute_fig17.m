function archive = simtwc_compute_fig17()
%SIMTWC_COMPUTE_FIG17 Endpoint KKT multipliers and normalized constraints.

s=simtwc_setup(); spec=s.fig17;
if ~isscalar(spec.mc_trials) || spec.mc_trials<1 || spec.mc_trials~=fix(spec.mc_trials)
    error('s.fig17.mc_trials must be a positive integer.');
end
expected_seeds=spec.base_seed+(0:spec.mc_trials-1);
if ~isequal(spec.seeds(:).',expected_seeds)
    error('s.fig17.seeds is inconsistent with base_seed and mc_trials.');
end

data=simtwc_initialize_data(17,spec,{'costar'});
path=simtwc_archive_path(17);
nT=numel(spec.seeds); J=spec.U+spec.M;
if ~isfield(data,'runs') || numel(data.runs)~=nT
    data.runs=cell(1,nT);
end
for t=1:nT
    if ok(data.runs{t}); continue; end
    seed=spec.seeds(t);
    fprintf('[Fig17] seed=%d, COSTAR endpoint multiplier audit\n',seed);
    p=simtwc_make_problem(spec,seed);
    data.runs{t}=simtwc_try_run('costar',p,false,struct());
    simtwc_save_atomic(path,data);
end

lambda=nan(J,nT); constraints=nan(J,nT); active=false(J,nT);
feasible=false(1,nT); crlb=nan(1,nT); kkt=nan(1,nT);
sinr_margin_dB=nan(spec.U,nT);
power_margin_percent=nan(spec.M,nT);
for t=1:nT
    q=data.runs{t}; if ~ok(q); continue; end
    lambda(:,t)=q.kkt.lambda(:);
    constraints(:,t)=q.kkt.constraints(:);
    active(:,t)=q.kkt.active(:);
    feasible(t)=logical(q.feasible);
    crlb(t)=q.metrics.crlb;
    kkt(t)=q.kkt.r_kkt;
    sinr_margin_dB(:,t)=10*log10(max(q.metrics.sinr(:),realmin))- ...
        spec.gamma_dB;
    power_margin_percent(:,t)=-100*q.kkt.constraints(spec.U+1:end);
end
data.seeds=spec.seeds;
data.lambda=lambda;
data.constraints=constraints;
data.active=active;
data.complementarity=abs(lambda.*constraints);
data.feasible=feasible;
data.crlb=crlb;
data.kkt_residual=kkt;
data.sinr_margin_dB=sinr_margin_dB;
data.power_margin_percent=power_margin_percent;
data.sinr_indices=1:spec.U;
data.power_indices=spec.U+(1:spec.M);
data.complete=all(cellfun(@ok,data.runs));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data);
archive=path;
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
