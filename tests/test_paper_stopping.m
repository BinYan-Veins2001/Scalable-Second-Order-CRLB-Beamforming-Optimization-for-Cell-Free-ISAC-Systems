function report=test_paper_stopping()
%TEST_PAPER_STOPPING Analytic endpoints and one-step NC restoration gate.
% Both examples use the actual objective, constraint, geometry, SPSC and
% solver entry points. No oracle or trial point is mocked.

% With a=b=p, f=log(2)-log(p). The power endpoint p=1 is KKT,
% with lambda_power=1 and zero phase curvature.
[cfg,channels]=small_configuration(1);
smat=struct('Aa',1,'Ab',1,'Aab',0);
options=struct('max_iter',4,'epsilon_KKT',1e-6,'epsilon_H',1e-5, ...
    'nu0',1e-6,'verbose',false);
[endpoint,endpoint_hist]=solve_costar(cfg,channels,smat,{1},options);
assert(endpoint.success && endpoint.numerical_joint_pass);
assert(endpoint.kkt_residual<=options.epsilon_KKT);
assert(endpoint.primal_residual<=0);
assert(endpoint.curvature_lower_bound>=-options.epsilon_H);
assert(endpoint.critical_cone_tolerance<=options.epsilon_H);
assert(endpoint.certificate.oracle_is_at_final_point);
assert(endpoint.work.meo_calls==0 && endpoint.work.meo_hvp_calls==0);
assert(endpoint.work.structured_checks>=1);
assert(abs(endpoint.lambda(2)-1)<1e-10 && endpoint.lambda(1)==0);
assert(abs(endpoint.metrics.crlb-2)<1e-12);
assert(all(diff(endpoint_hist.rtr.objective)<=1e-12));

% For a=b=w'*diag(1,2)*w, w=e1 is a power-constrained saddle.
% The Lagrangian Hessian in [Re(w);Im(w)] coordinates is
% diag(4,-2,0,-2), with strictly positive power multiplier.
[cfg,channels]=small_configuration(2);
smat=struct('Aa',diag([1,2]),'Ab',diag([1,2]),'Aab',zeros(2));
initial=[1;0]; initial_objective=log(2);
options=struct('max_iter',1,'restoration_max_iter',0, ...
    'epsilon_KKT',1e-6,'epsilon_H',1e-5,'nu0',1e-6, ...
    'delta0',0.05,'normal_inward_margin',1,'feas_tol',0,'verbose',false);
[escaped,escape_hist]=solve_costar(cfg,channels,smat,{initial},options);
final_objective=costar_paper_value(cfg,smat,escaped.Wmat);
[constraints,~]=costar_sosp_constraints(cfg,channels,escaped.Wmat);
assert(final_objective<initial_objective-1e-6);
assert(max(constraints)<=0 && escaped.primal_residual<=0);
assert(escaped.work.normal_newton_accepted>=1);
assert(escaped.work.negative_curvature_steps==1);
assert(escaped.work.accepted_steps==1);
assert(escaped.work.meo_calls==0 && escaped.work.meo_hvp_calls==0);
assert(any(strcmp(escape_hist.rtr.step_type,'spsc_negative_curvature')));
assert(all(diff(escape_hist.rtr.objective)<=1e-12));
% The one-step budget stops at the new point, whose old certificate cannot
% be reused. An accepted negative-curvature step is not a final SOSP flag.
assert(~escaped.success && ~escaped.numerical_joint_pass);
assert(~escaped.certificate.oracle_is_at_final_point);
assert(~escaped.certificate.deterministic_numerical);
assert(isinf(escaped.curvature_lower_bound) && escaped.curvature_lower_bound<0);

report=struct('passed',true,'analytic_kkt_certified',endpoint.success, ...
    'analytic_kkt_residual',endpoint.kkt_residual, ...
    'analytic_curvature_lower_bound',endpoint.curvature_lower_bound, ...
    'nc_objective_decrease',initial_objective-final_objective, ...
    'nc_primal_residual',max([0;constraints]), ...
    'normal_newton_accepted',escaped.work.normal_newton_accepted, ...
    'stale_certificate_rejected',~escaped.success, ...
    'randomized_krylov_calls',endpoint.work.meo_calls+escaped.work.meo_calls);
disp(report);
end

function [cfg,channels]=small_configuration(Nt)
cfg=default_config();
cfg.sys.M=1; cfg.sys.Nt=Nt; cfg.sys.Nr=1; cfg.sys.U=1;
cfg.power.Pmax=1; cfg.comm.noise_power=1;
cfg.comm.Gamma=0.1; cfg.comm.Gamma_dB=-10;
cfg.radio.effective_bandwidth_sq=1;
cfg.radio.bandwidth=sqrt(12); cfg.radio.c=1; cfg.radio.frame_duration=1;
cfg.sensing.N0=8*pi^2; cfg.alm.objective_type='log_crlb_ratio';
cfg.alm.det_floor=1e-30; cfg.costar.paper=struct();
h=zeros(Nt,1); h(1)=1;
channels=struct('h_stack',h,'h',{{h}});
end
