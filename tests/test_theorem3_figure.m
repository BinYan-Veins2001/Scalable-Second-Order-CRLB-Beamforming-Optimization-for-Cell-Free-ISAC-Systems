function report = test_theorem3_figure()
%TEST_THEOREM3_FIGURE Algebraic oracle and dense/analytic-HVP regression tests.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
s=simtwc_setup();
spec=s.default; spec.M=4; spec.U=2; spec.Nt=2; spec.Nr=2;
p=simtwc_make_problem(spec,2026);
assert(p.initialization.success,p.initialization.message);
o=simtwc_theorem3_options(struct('cone_samples',20));
[a,op]=simtwc_theorem3_audit(p,p.initialization.z,o,true);
assert(a.hvp_relative_error<1e-10);
assert(norm(op.Z'*op.Z-eye(size(op.Z,2)),'fro')<1e-10);

H=diag([-0.8,linspace(0.2,1,19)]);
settings=struct('epsilon_H',1e-4,'norm_bound',1,'failure_probability',1e-3, ...
    'effective_dim',20,'record_history',true,'overlap_safety_factor',1,'seed',2026);
negative=costar_sosp_meo(@(v) H*v,20,@(v) v,settings);
assert(strcmp(negative.status,'negative_curvature'));
assert(negative.rayleigh<=-settings.epsilon_H/2);
assert(negative.eigenvalue_lower_bound<=-0.8+1e-10);
assert(~isempty(negative.trace.basis_dim));

H=diag(linspace(0,1,20));
positive=costar_sosp_meo(@(v) H*v,20,@(v) v,settings);
assert(positive.conditional_psd);
assert(positive.eigenvalue_lower_bound<=1e-10);
assert(positive.eigenvalue_lower_bound>=-settings.epsilon_H);
settings.max_iter=4;
capped=costar_sosp_meo(@(v) H*v,20,@(v) v,settings);
assert(capped.budget_capped && ~capped.conditional_psd);

% A weak active constraint must not be promoted to a critical-cone certificate.
[sol,~]=polish_costar_sosp(p.cfg,p.channels,p.smat,p.initialization.wcell, ...
    struct('max_iter',2,'curvature_space','active_tangent'));
assert(strcmp(sol.options.curvature_space,'active_tangent'));
assert(~sol.certificate.floating_point_validated);
report=struct('dense_hvp_error',a.hvp_relative_error, ...
    'negative_oracle',negative,'positive_oracle',positive,'capped_oracle',capped);
save(fullfile(root,'tests','theorem3_figure_report.mat'),'report');
fprintf('THEOREM 3 FIGURE TEST PASSED: dense HVP error %.3g\n',a.hvp_relative_error);
end
