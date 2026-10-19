function archive = simtwc_compute_fig16()
%SIMTWC_COMPUTE_FIG16 COSTAR iteration and horizontal-curvature audit.

s=simtwc_setup(); spec=s.fig16;
data=simtwc_initialize_data(16,spec,{'costar'});
path=simtwc_archive_path(16);
if ~isfield(data,'run') || ~ok(data.run)
    p=simtwc_make_problem(spec,spec.seed);
    run=simtwc_try_run('costar',p,true,struct());
    if ~ok(run)
        error('Fig. 16 COSTAR run failed: %s',run.message);
    end
    if ~isfield(run.history,'raw') || ~isfield(run.history.raw,'rtr')
        error('Fig. 16 requires the raw COSTAR RTR history.');
    end
    rtr=run.history.raw.rtr;
    required={'iter','crlb','tangent_grad_norm'};
    for j=1:numel(required)
        if ~isfield(rtr,required{j}) || isempty(rtr.(required{j}))
            error('Fig. 16 history field %s is missing.',required{j});
        end
    end
    n=min([numel(rtr.iter),numel(rtr.crlb),numel(rtr.tangent_grad_norm)]);
    data.seed=spec.seed;
    data.run=run;
    data.iter=real(rtr.iter(1:n)).';
    data.crlb=real(rtr.crlb(1:n)).';
    data.tangent_grad_norm=real(rtr.tangent_grad_norm(1:n)).';
    data.accepted=logical(field_prefix(rtr,'accepted',n,false));
    data.complete=false;
    simtwc_save_atomic(path,data);
end

rtr=data.run.history.raw.rtr;
n=numel(data.iter);
if ~isfield(rtr,'z_snapshots') || size(rtr.z_snapshots,2)<n
    error('Fig. 16 needs one saved beamformer for each RTR iteration.');
end
opts=struct('hessian_model','exact_lagrangian', ...
    'lanczos_max_iter',s.fig12.lanczos_max_iter, ...
    'lanczos_tol',s.fig12.lanczos_tol, ...
    'seed',spec.seed+7000+6, ...
    'sinr_active_tol',s.fig12.sinr_active_tol, ...
    'power_active_tol',s.fig12.power_active_tol);
valid_cache=isfield(data,'curvature_opts') && isequaln(data.curvature_opts,opts) && ...
    isfield(data,'curvature_done') && numel(data.curvature_done)==n && ...
    isfield(data,'lambda_min_horizontal') && ...
    numel(data.lambda_min_horizontal)==n;
if ~valid_cache
    data.lambda_min_horizontal=nan(1,n);
    data.curvature_kkt=nan(1,n);
    data.curvature_lanczos_iter=nan(1,n);
    data.curvature_horizontal_dim=nan(1,n);
    data.curvature_done=false(1,n);
    data.curvature_opts=opts;
    data.complete=false;
end
if all(data.curvature_done)
    archive=path;
    return;
end

p=simtwc_make_problem(spec,spec.seed);
for k=find(~data.curvature_done)
    fprintf('[Fig16] horizontal curvature at RTR iteration %d/%d\n',k,n);
    z=real(rtr.z_snapshots(:,k));
    diagnostic=evaluate_phase2_quotient_second_order_stationarity( ...
        p.cfg,p.channels,p.smat,z,opts);
    data.lambda_min_horizontal(k)=diagnostic.lambda_min_horizontal;
    data.curvature_kkt(k)=diagnostic.kkt_residual;
    data.curvature_lanczos_iter(k)=diagnostic.lanczos.iter;
    data.curvature_horizontal_dim(k)=diagnostic.horizontal_dim;
    data.curvature_done(k)=true;
    if mod(k,10)==0 || k==n
        simtwc_save_atomic(path,data);
    end
end
data.complete=all(data.curvature_done);
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data);
archive=path;
end

function v=field_prefix(s,name,n,fallback)
if isfield(s,name) && numel(s.(name))>=n
    v=s.(name)(1:n);
else
    v=repmat(fallback,n,1);
end
v=v(:);
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
