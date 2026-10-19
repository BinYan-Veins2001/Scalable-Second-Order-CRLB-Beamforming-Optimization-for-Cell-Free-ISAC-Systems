function report=test_curvature_oracle()
%TEST_CURVATURE_ORACLE Hidden curvature, PSD enclosure, and capped-budget gate.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
before=rng; cleanup=onCleanup(@() rng(before)); %#ok<NASGU>
rng(913,'twister'); n=30; [Q,~]=qr(randn(n));
e=linspace(0.01,2,n); e(1)=-0.05; A=Q*diag(e)*Q';
options=struct('epsilon_H',1e-3,'norm_bound',2, ...
    'failure_probability',1e-3,'seed',1926,'max_iter',Inf);
negative=costar_sosp_meo(@(v) A*v,n,@(v) v,options);
assert(strcmp(negative.status,'negative_curvature'));
assert(negative.negative_direction'*A*negative.negative_direction<=-options.epsilon_H/2);
positive=costar_sosp_meo(@(v) (A+0.1*eye(n))*v,n,@(v) v,options);
assert(positive.conditional_psd && positive.eigenvalue_lower_bound>=-options.epsilon_H);
options.max_iter=1;
capped=costar_sosp_meo(@(v) A*v,n,@(v) v,options);
assert(~capped.conditional_psd && capped.budget_capped);
H=diag([-1 2]); g=[1;0];
[step_on,cg_on]=trust_region_tcg(g,@(v) H*v,1,struct('detect_negative_curvature',true));
[step_off,cg_off]=trust_region_tcg(g,@(v) H*v,1,struct('detect_negative_curvature',false));
assert(cg_on.negative_curvature && ~cg_off.negative_curvature);
assert(strcmp(cg_off.reason,'negative_curvature_disabled') && norm(step_off)<norm(step_on));
report=struct('negative',negative,'positive',positive,'capped',capped, ...
    'tcg_on',cg_on,'tcg_off',cg_off,'passed',true);
save(fullfile(root,'tests','curvature_oracle_report.mat'),'report');
fprintf('CURVATURE ORACLE TEST PASSED (negative / PSD / capped-budget)\n');
end
