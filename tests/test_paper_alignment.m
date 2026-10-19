function report = test_paper_alignment()
%TEST_PAPER_ALIGNMENT Small physical problem with an independent dense audit.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root),'-begin');
setenv('COSTAR_QUICK',''); setenv('COSTAR_SMOKE','');
assert(startsWith(which('solve_costar'),root),'Wrong project on path.');
spec=struct('M',4,'U',2,'Nt',2,'Nr',2,'gamma_dB',10,'min_ap_distance',10);
p=simtwc_make_problem(spec,2026);
assert(p.initialization.success,p.initialization.message);
initial=p.initialization.Wmat;
initial_metrics=evaluate_beamforming_solution(p.cfg,p.channels,p.smat, ...
    p.initialization.wcell);
[sol,hist]=solve_costar(p.cfg,p.channels,p.smat,p.initialization.wcell);
assert(isequal(initial,p.initialization.Wmat),'Initial beamformer was modified.');
assert(isfield(sol,'lambda') && isfield(sol,'active'), ...
    'The paper-aligned endpoint must expose its own multipliers/working set.');
audit=benchmark_endpoint_audit(p,sol.z,struct('lambda',sol.lambda, ...
    'active',sol.active,'dense_max_dim',2048));
assert(isfinite(audit.crlb) && audit.crlb>0,'Invalid physical CRLB.');
assert(audit.crlb<=initial_metrics.crlb*(1+1e-10),'CRLB increased.');
assert(audit.primal<=1e-10,'Endpoint violates normalized feasibility.');
assert(abs(audit.kkt-sol.kkt_residual)<=1e-8*max(1,audit.kkt), ...
    'Reported KKT residual disagrees with independent calculation.');
assert(isfield(hist,'rtr') && isfield(hist.rtr,'crlb'), ...
    'Canonical committed-iterate history is missing.');
crlb=hist.rtr.crlb(:); crlb=crlb(isfinite(crlb));
assert(all(diff(crlb)<=1e-10*max(1,initial_metrics.crlb)), ...
    'Committed objective trajectory is not monotone.');
assert(isfield(sol.work,'safeguard_calls'),'First-order safeguard was not instrumented.');
assert(isfield(sol.work,'structured_checks'),'SPSC was not instrumented.');
if isfield(sol.work,'meo_calls')
    assert(sol.work.meo_calls==0,'Randomized Krylov is still used.');
end
if sol.success
    assert(audit.kkt<=1e-5*(1+1e-6),'False first-order certificate.');
    assert(audit.lambda_min_tangent>=-1e-5-1e-8, ...
        'False tangent-curvature certificate.');
    assert(sol.curvature_lower_bound>=-1e-5,'Invalid reported lower bound.');
    assert(sol.critical_cone_tolerance<=1e-5 || ...
        (isfield(sol.certificate,'oracle') && sol.certificate.oracle.ambient), ...
        'Success lacks the required cone/full-space curvature check.');
end
report=struct('passed',true,'solver_success',sol.success,'audit',audit, ...
    'message',sol.message,'work',sol.work);
save(fullfile(root,'tests','paper_alignment_report.mat'),'report','-v7.3');
fprintf('PAPER ALIGNMENT small: success=%d CRLB=%.9g KKT=%.3g eigT=%.3g; %s\n', ...
    sol.success,audit.crlb,audit.kkt,audit.lambda_min_tangent,sol.message);
end
