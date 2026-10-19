function results = test_paper_geometry()
%TEST_PAPER_GEOMETRY Analytic geometry and safeguard identities.
results=struct(); o=struct();

% A weakly active SINR boundary must retain its feasible inward descent.
g=[-2;0]; G=[-2;0]; c=0;
s=costar_paper_geometry(g,c,G,1,1e-5,o);
assert(s.success && s.lambda==0 && norm(s.gt)<1e-12);
assert(G'*(-s.zeta)<=1e-12 && abs(g'*(-s.zeta)+norm(s.zeta)^2)<1e-12);
assert(norm(s.normal_inward-[0.5;0])<1e-10);
results.weak_active_descent=true;

% Scaling a constraint does not change its tangent space or stationarity.
g=[-2;-3]; G=diag([2,300]); c=[0;0];
s=costar_paper_geometry(g,c,G,1,1e-5,o);
assert(s.success && norm(s.zeta)<1e-11);
assert(norm(s.lambda-[1;0.01])<1e-11);
assert(abs(s.sigma_min-2)<1e-12 && s.full_column_rank);
assert(norm(s.normal_inward-[-0.5;-1/300])<1e-9);
results.unscaled_geometry=true;

% Dependent, feasible normals: strongest inequality determines the norm.
in=costar_paper_inward([1,2;0,0],o);
assert(in.success && norm(in.direction-[-1;0])<1e-8);
assert(in.complementarity_residual<1e-9);
results.dependent_inward=true;

% Incompatible near-active constraints must be released by thresholding.
o2=struct('tau_bar',1e-3,'sinr_active_tol',1e-3, ...
    'power_active_tol',1e-3,'inward_max_sweeps',50);
s=costar_paper_geometry([0;1],[-1e-4;-1e-4],[1,-1;0,0],1,1e-5,o2);
assert(s.success && ~any(s.active) && s.work.threshold_refinements>0);
results.incompatible_near_active_release=true;

% Truly incompatible active normals remain explicitly uncertified.
o3=struct('threshold_max_iter',2,'inward_max_sweeps',10);
s=costar_paper_geometry([0;1],[0;0],[1,-1;0,0],1,1e-5,o3);
assert(~s.success && contains(s.message,'uncertified'));
results.infeasible_budget_status=true;

% Empty sets and the strict c > -threshold selection rule.
o4=struct('tau_bar',1e-4,'kappa_tau',0.5);
s=costar_paper_geometry([1;0],-1e-4,[1;0],1,1e-5,o4);
assert(s.success && ~any(s.active) && isempty(s.Q));
assert(norm(s.normal_inward)==0 && s.positive_multipliers);
results.empty_and_strict_selection=true;

% A narrow feasible cone checks the dual support Newton polish.
G=[1,-1;0,0.02];
in=costar_paper_inward(G,o);
expected=-(G')\ones(2,1);
assert(in.success && norm(in.direction-expected)/norm(expected)<1e-7);
assert(max(G'*in.direction+1)<1e-9);
results.narrow_inward_cone=true;
disp(results);
end
