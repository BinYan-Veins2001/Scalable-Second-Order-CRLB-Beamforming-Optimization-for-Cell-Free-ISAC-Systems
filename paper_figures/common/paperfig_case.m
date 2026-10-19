function [p,run]=paperfig_case(regime,spec,seed,method,o)
%PAPERFIG_CASE Share fresh endpoints across figures without legacy archives.
s=paperfig_setup(); signature=paperfig_signature(); physical_signature=paperfig_audit_signature(0);
folder=fullfile(s.case_dir,sprintf('%s_seed%d',regime,seed));
if ~isfolder(folder); mkdir(folder); end
problem_path=fullfile(folder,'problem.mat');
if isfile(problem_path)
    old=load(problem_path);
else
    old=struct();
end
if isfield(old,'signature') && strcmp(old.signature,signature) && isequaln(old.spec,spec) && ...
        isfield(old,'physical_signature') && strcmp(old.physical_signature,physical_signature)
    p=old.p;
else
    p=simtwc_make_problem(spec,seed);
    save(problem_path,'p','spec','signature','physical_signature','-v7.3');
end
path=fullfile(folder,[method '.mat']);
if o.resume && isfile(path)
    old=load(path);
    if isfield(old,'signature') && strcmp(old.signature,signature) && ...
            isequaln(old.solver_options,o.solver_options) && isfield(old,'problem_spec') && ...
            isequaln(old.problem_spec,spec) && isfield(old,'seed') && old.seed==seed && isfield(old,'run') && ...
            isfield(old,'physical_signature') && strcmp(old.physical_signature,physical_signature)
        run=old.run; return;
    end
end
fprintf('[Paper case] %s, seed %d, %s\n',regime,seed,method);
if strcmp(method,'costar')
    assert(p.initialization.success,p.initialization.message);
    opts=o.solver_options; opts.seed=seed+4000; opts.record_snapshots=true;
    timer=tic;
    [sol,hist]=solve_costar(p.cfg,p.channels,p.smat,p.initialization.wcell,opts);
    time=toc(timer);
    run=struct('key',method,'seed',seed,'status','ok','solution',sol, ...
        'history',hist,'Wmat',sol.Wmat,'metrics',sol.metrics, ...
        'runtime_median_s',time+p.initialization.runtime_s, ...
        'optimization_runtime_s',time,'initialization_runtime_s',p.initialization.runtime_s, ...
        'feasible',sol.primal_residual<=sol.options.feas_tol, ...
        'message',sol.message,'solver_success',sol.success);
else
    run=simtwc_try_run(method,p,false,struct());
end
solver_options=o.solver_options;
problem_spec=spec;
save(path,'run','signature','solver_options','problem_spec','seed','physical_signature','-v7.3');
end
