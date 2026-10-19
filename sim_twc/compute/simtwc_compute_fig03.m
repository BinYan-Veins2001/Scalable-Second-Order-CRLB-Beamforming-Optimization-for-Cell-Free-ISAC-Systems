function archive = simtwc_compute_fig03()
%SIMTWC_COMPUTE_FIG03 Seven-method paired endpoint comparison.

s=simtwc_setup(); spec=s.fig03; methods=s.methods; seeds=spec.seeds;
data=simtwc_initialize_data(3,spec,methods); path=simtwc_archive_path(3);
if ~isfield(data,'runs') || ~isequal(size(data.runs),[numel(methods),numel(seeds)])
    data.runs=cell(numel(methods),numel(seeds));
end
for t=1:numel(seeds)
    p=simtwc_make_problem(spec,seeds(t));
    for a=1:numel(methods)
        if ok(data.runs{a,t}); continue; end
        fprintf('[Fig03] seed=%d, %s\n',seeds(t),methods{a});
        data.runs{a,t}=simtwc_try_run(methods{a},p,true,struct());
        simtwc_save_atomic(path,data);
    end
end
data.seeds=seeds;
[data.final_crlb,data.runtime_s,data.feasible,data.kkt_residual]= ...
    endpoint_arrays(data.runs);
data.complete=all(cellfun(@ok,data.runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end

function [C,T,F,K]=endpoint_arrays(runs)
[nA,nT]=size(runs); C=nan(nA,nT); T=nan(nA,nT);
F=false(nA,nT); K=nan(nA,nT);
for a=1:nA
    for t=1:nT
        r=runs{a,t}; if ~ok(r); continue; end
        C(a,t)=r.metrics.crlb; T(a,t)=r.runtime_median_s;
        F(a,t)=r.feasible; K(a,t)=r.kkt.r_kkt;
    end
end
end
