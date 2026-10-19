function results = test_guard_scaled_inward()
%TEST_GUARD_SCALED_INWARD Analytic tests for signed directional right sides.
results=struct(); o=struct('inward_tol',1e-9,'inward_max_sweeps',200);

[v,s]=costar_guard_scaled_inward(eye(2),[2;-3],.5,o);
check(v,s,[-2.5;0]); results.mixed_signs=true;

G=[1,0,1;0,1,1];
[v,s]=costar_guard_scaled_inward(G,[1;-2;.5],0,o);
check(v,s,[-1;0]); results.inactive_negative_rhs=true;

% A negative right side can still bind when another inequality forces it.
G=[1,-1;0,0];
[v,s]=costar_guard_scaled_inward(G,[-1;.5],0,o);
check(v,s,[.5;0]); results.opposing_feasible_inequalities=true;

[v,s]=costar_guard_scaled_inward([1,2;0,0],[1;3],0,o);
check(v,s,[-1.5;0]); results.dependent_normals=true;

G=[1,-1;0,.02]; b=[1;-.5];
[v,s]=costar_guard_scaled_inward(G,b,0,o);
check(v,s,[-1;-25]);
scale=[1e6;1e-3];
[vs,ss]=costar_guard_scaled_inward(G.*scale.',b.*scale,0,o);
check(vs,ss,[-1;-25]); results.narrow_cone_and_scaling=true;

[v,s]=costar_guard_scaled_inward(eye(3),[-2;-1;0],0,o);
check(v,s,zeros(3,1)); results.zero_optimum=true;
[v,s]=costar_guard_scaled_inward(zeros(3,0),zeros(0,1),1,o);
check(v,s,zeros(3,1)); results.empty_set=true;

% Unlike a common maximum coefficient, the AP-like second coordinate
% receives only the curvature correction it requires.
[v,s]=costar_guard_scaled_inward(eye(2),[5;2],1,o);
check(v,s,[-6;-3]);
assert(norm(v)<norm([-6;-6])); results.per_constraint_correction=true;

[~,s]=costar_guard_scaled_inward([1,-1;0,0],[1;1],0, ...
    struct('inward_max_sweeps',12));
assert(~s.success); results.infeasible_is_unresolved=true;
[~,s]=costar_guard_scaled_inward([0;0],1,0,o);
assert(~s.success); results.zero_gradient_failure=true;
[v,s]=costar_guard_scaled_inward([1,0;0,0],[1;-1],0,o);
check(v,s,[-1;0]); results.harmless_zero_gradient=true;
[~,s]=costar_guard_scaled_inward(eye(2),ones(2,1),0, ...
    struct('inward_max_sweeps',0));
assert(~s.success); results.budget_cannot_certify=true;
disp(results);
end

function check(v,s,expected)
assert(s.success,s.message);
assert(norm(v-expected)<=1e-7*max(1,norm(expected)));
assert(s.scaled_primal_residual<=1e-9);
assert(s.stationarity_residual<=1e-9 && s.complementarity_residual<=1e-9);
assert(all(s.multipliers>=0));
end
