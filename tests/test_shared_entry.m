function report=test_shared_entry()
%TEST_SHARED_ENTRY Actual small-system execution of all seven methods.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
previous=getenv('COSTAR_SMOKE'); cleanup=onCleanup(@() setenv('COSTAR_SMOKE',previous)); %#ok<NASGU>
setenv('COSTAR_SMOKE','1'); s=simtwc_setup();
p=simtwc_make_problem(s.default,2026);
assert(p.initialization.success,p.initialization.message);
assert(p.initialization.info.options.max_refine_iter==0);
assert(p.initialization.info.num_accepted_refine_steps==0);
before=p.initialization.Wmat; methods=s.methods;
report=struct('methods',{methods},'runs',{cell(size(methods))}, ...
    'shared_entry_pass',false,'created_at',char(datetime('now')));
for j=1:numel(methods)
    q=simtwc_run_method(methods{j},p);
    q.status='ok'; report.runs{j}=q;
    if any(strcmp(methods{j},s.iterative_methods))
        assert(isequal(q.initial_Wmat,before),'Different algorithm entry detected.');
        assert(abs(q.elapsed_s-q.initialization_runtime_s-q.optimization_runtime_s)<1e-10);
        assert(isequal(q.history.z_snapshots(:,1),p.initialization.z));
        if isfield(q.solution,'initial_Wmat')
            assert(norm(q.solution.initial_Wmat-before,'fro')<1e-12*max(1,norm(before,'fro')));
        elseif strcmp(methods{j},'costar')
            assert(isequal(q.history.raw.fast.z_snapshots(:,1),p.initialization.z));
        elseif any(strcmp(methods{j},{'mrfp_sca','sdr_sca'}))
            % SOCP phase normalization and covariance recovery may rotate
            % each column, but must preserve the supplied rank-one input.
            actual=real_vector_to_complex_matrix(q.history.z_snapshots(:,2), ...
                p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
            for u=1:p.cfg.sys.U
                target=before(:,u)*before(:,u)';
                assert(norm(actual(:,u)*actual(:,u)'-target,'fro')<1e-10*max(1,norm(target,'fro')));
            end
        end
    end
    assert(all(diff(q.history.time_s)>=0));
    fprintf('TEST %s: CRLB %.8g, KKT %.3g, feasible %d, %.3fs\n', ...
        methods{j},q.metrics.crlb,q.kkt.r_kkt,q.feasible,q.elapsed_s);
end
assert(isequal(p.initialization.Wmat,before),'Shared point was overwritten.');
q=report.runs{strcmp(methods,'costar')};
assert(strcmp(q.algorithm_version,'COSTAR_SOSP_ExternalSharedZF_20261003'));
assert(isfield(q.history.raw,'terminal') && isfield(q.history.raw,'fast'));
assert(q.solution.work.fd_hvp_calls==0);
assert(numel(q.history.raw.rtr.iter)==size(q.history.raw.rtr.z_snapshots,2));
try
    solve_costar(p.cfg,p.channels,p.smat,[]);
    error('Test:MissingError','COSTAR unexpectedly performed initialization.');
catch exception
    assert(strcmp(exception.identifier,'COSTAR:ExternalInitializationRequired'));
end
report.shared_entry_pass=true;
save(fullfile(root,'tests','shared_entry_report.mat'),'report','-v7.3');
fprintf('SHARED ENTRY TEST PASSED\n');
end
