function archive=simtwc_compute_fig19()
%SIMTWC_COMPUTE_FIG19 Compare common phase-reduced endpoint curvatures.

s=simtwc_setup(); spec=s.fig19;
source_path=simtwc_compute_fig12();
source=load(source_path);
if ~source.complete
    error('Fig. 19 needs a complete Fig. 12 endpoint-curvature archive.');
end
methods=cellstr(string(spec.methods));
data=simtwc_initialize_data(19,spec,methods);
path=simtwc_archive_path(19);
if isfield(data,'complete') && data.complete && ...
        isfield(data,'source_completed_at') && ...
        strcmp(data.source_completed_at,source.completed_at)
    archive=path;
    return;
end

nA=numel(methods); nT=numel(spec.seeds);
data.lambda_min_horizontal=nan(nA,nT);
data.kkt_residual=nan(nA,nT);
data.feasible=false(nA,nT);
data.lanczos_iter=nan(nA,nT);
for a=1:nA
    source_a=find(strcmp(source.methods,methods{a}),1);
    if isempty(source_a)
        error('Method %s is missing from Fig. 12.',methods{a});
    end
    for t=1:nT
        source_t=find(source.seeds==spec.seeds(t),1);
        if isempty(source_t)
            error('Seed %d is missing from Fig. 12.',spec.seeds(t));
        end
        q=source.runs{source_a,source_t};
        if ~isstruct(q) || ~isfield(q,'status') || ~strcmp(q.status,'ok')
            error('Fig. 12 audit failed for %s, seed %d.', ...
                methods{a},spec.seeds(t));
        end
        data.lambda_min_horizontal(a,t)=q.diagnostic.lambda_min_horizontal;
        data.kkt_residual(a,t)=q.diagnostic.kkt_residual;
        data.feasible(a,t)=q.source_feasible;
        data.lanczos_iter(a,t)=q.diagnostic.lanczos.iter;
    end
end
data.seeds=spec.seeds;
data.methods=methods;
data.source_completed_at=source.completed_at;
data.complete=true;
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data);
archive=path;
end
