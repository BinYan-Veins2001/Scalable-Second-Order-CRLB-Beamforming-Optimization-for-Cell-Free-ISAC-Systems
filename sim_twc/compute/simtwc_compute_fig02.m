function archive = simtwc_compute_fig02()
%SIMTWC_COMPUTE_FIG02 Complete wall-clock CRLB histories for five methods.

s=simtwc_setup(); spec=s.fig02; methods=s.iterative_methods;
data=simtwc_initialize_data(2,spec,methods); path=simtwc_archive_path(2);
if ~isfield(data,'runs'); data.runs=cell(1,numel(methods)); end
p=simtwc_make_problem(spec,spec.seed);
for a=1:numel(methods)
    if ok(data.runs{a}); continue; end
    fprintf('[Fig02] %s\n',methods{a});
    data.runs{a}=simtwc_try_run(methods{a},p,true,struct());
    simtwc_save_atomic(path,data);
end
data.scenario=p.scenario;
data.complete=all(cellfun(@ok,data.runs));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
