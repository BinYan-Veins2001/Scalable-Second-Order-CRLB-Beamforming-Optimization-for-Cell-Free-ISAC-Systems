function report = test_theorem3_terminal_random()
%TEST_THEOREM3_TERMINAL_RANDOM Exercise the actual Gaussian tangent oracle.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
old_quick=getenv('COSTAR_QUICK'); old_smoke=getenv('COSTAR_SMOKE');
cleanup=onCleanup(@() restore(old_quick,old_smoke)); %#ok<NASGU>
setenv('COSTAR_QUICK',''); setenv('COSTAR_SMOKE','');
s=simtwc_setup(); spec=s.default; spec.gamma_dB=20;
p=simtwc_make_problem(spec,2026); assert(p.initialization.success);
term=struct('curvature_space','active_tangent','use_structured_psd',false, ...
    'meo_record_history',true,'meo_overlap_safety_factor',1,'max_iter',200);
[sol,hist]=solve_costar(p.cfg,p.channels,p.smat,p.initialization.wcell, ...
    struct('terminal_options',term));
assert(sol.work_terminal.meo_calls>=1,'Random terminal branch was not exercised.');
assert(sol.conditional_sosp_reached,sol.message);
assert(~sol.certificate.oracle.ambient);
assert(~sol.certificate.floating_point_validated);
assert(~isempty(sol.certificate.oracle.trace.basis_dim));
o=simtwc_theorem3_options(); p.cfg.costar.sosp=sol.terminal_options;
a=simtwc_theorem3_audit(p,sol.z,o,true);
assert(a.joint_ratio<=1 && a.positive_multipliers && a.licq);
report=struct('solution',sol,'history',hist,'audit',a);
save(fullfile(root,'tests','theorem3_terminal_random_report.mat'),'report','-v7.3');
fprintf('REAL TANGENT ORACLE PASSED: KKT %.3g, lambda_T %.3g, search %d/%d\n', ...
    a.kkt,a.dense_lambda_min,sol.certificate.oracle.iterations,a.tangent_dimension);
end

function restore(quick,smoke)
setenv('COSTAR_QUICK',quick); setenv('COSTAR_SMOKE',smoke);
end
