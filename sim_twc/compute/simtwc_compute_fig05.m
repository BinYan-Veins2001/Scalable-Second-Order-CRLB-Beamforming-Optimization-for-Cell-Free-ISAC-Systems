function archive = simtwc_compute_fig05()
%SIMTWC_COMPUTE_FIG05 COSTAR trust-region and active-set diagnostics.

s=simtwc_setup(); spec=s.fig05;
data=simtwc_initialize_data(5,spec,{'costar'}); path=simtwc_archive_path(5);
if data.complete; archive=path; return; end
p=simtwc_make_problem(spec,spec.seed);
run=simtwc_try_run('costar',p,true,struct());
if ~ok(run); error('Fig. 5 COSTAR run failed: %s',run.message); end
rtr=run.history.raw.rtr;
data.run=run;
data.iter=rtr.iter(:);
data.tangent_grad_norm=rtr.tangent_grad_norm(:);
data.delta=rtr.delta(:);
data.active_count=rtr.active_count(:);
data.active_membership=rtr.active_membership;
data.active_change_count=rtr.active_change_count(:);
data.complete=true; data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
