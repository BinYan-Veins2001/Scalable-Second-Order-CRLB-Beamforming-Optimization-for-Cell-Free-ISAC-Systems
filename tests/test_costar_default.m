function report=test_costar_default()
%TEST_COSTAR_DEFAULT Full optimizer settings at the default network size.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
quick=getenv('COSTAR_QUICK'); smoke=getenv('COSTAR_SMOKE');
cleanup=onCleanup(@() restore(quick,smoke)); %#ok<NASGU>
setenv('COSTAR_QUICK',''); setenv('COSTAR_SMOKE','');
s=simtwc_setup(); p=simtwc_make_problem(s.default,2026);
assert(p.initialization.success,p.initialization.message);
q=simtwc_run_method('costar',p);
assert(q.feasible,'Default COSTAR endpoint is infeasible.');
assert(q.solution.work.fd_hvp_calls==0);
assert(q.solution.terminal_options.epsilon_KKT==s.costar_options.epsilon_KKT);
assert(isequal(q.initial_Wmat,p.initialization.Wmat));
assert(isequal(q.history.raw.rtr.z_snapshots(:,end), ...
    complex_matrix_to_real_vector(q.Wmat)));
report=q;
save(fullfile(root,'tests','costar_default_report.mat'),'report','-v7.3');
fprintf('DEFAULT COSTAR: M=%d U=%d Nt=%d, CRLB %.9g, KKT %.3g, runtime %.3fs\n', ...
    p.cfg.sys.M,p.cfg.sys.U,p.cfg.sys.Nt,q.metrics.crlb,q.kkt.r_kkt,q.elapsed_s);
fprintf('Conditional SOSP flag %d; %s\n',q.conditional_sosp_reached,q.message);
end

function restore(quick,smoke)
setenv('COSTAR_QUICK',quick); setenv('COSTAR_SMOKE',smoke);
end
