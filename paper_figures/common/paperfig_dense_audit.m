function a=paperfig_dense_audit(p,sol,o)
%PAPERFIG_DENSE_AUDIT Independent dense tangent spectrum at the saved endpoint.
% Uses the solver's final NNLS multipliers and actual coordinated working set.
[geo,g,c,G,out]=paperfig_geometry(p,sol.z,sol.endpoint_geometry_nu,sol.options);
N=p.cfg.sys.M*p.cfg.sys.Nt; U=p.cfg.sys.U; d=numel(sol.z);
assert(d<=o.dense_max_dim,'Offline dense audit exceeds declared dimension cap.');
assert(norm(geo.lambda-sol.lambda)<=1e-9*max(1,norm(sol.lambda)) && ...
    isequal(geo.active,sol.active),'Saved endpoint geometry is stale.');
W=real_vector_to_complex_matrix(sol.z,N,U);
[~,info]=exact_lagrangian_hess_vec(p.cfg,p.channels,p.smat,W,zeros(N,U),geo.lambda);
alpha=info.objective.alpha; Hphi=info.objective.outer_hessian;
AA={p.smat.Aa,p.smat.Ab,p.smat.Aab}; beta=p.cfg.radio.effective_bandwidth_sq(:);
Ab=zeros(N); J=zeros(d,3);
for j=1:3
    Ab=Ab+alpha(j)*AA{j};
    J(:,j)=complex_matrix_to_real_vector(2*(AA{j}*W).*beta.');
end
H=zeros(d); C=p.channels.h_stack;
power=2*kron(geo.lambda(U+1:end)/p.cfg.power.Pmax,ones(p.cfg.sys.Nt,1));
for u=1:U
    weights=p.cfg.comm.Gamma*geo.lambda(1:U); weights(u)=-geo.lambda(u);
    B=2*beta(u)*Ab+2/p.cfg.comm.noise_power*(C.*weights.')*C'+diag(power);
    B=(B+B')/2; ii=(u-1)*N+(1:N); jj=N*U+ii;
    H(ii,ii)=real(B); H(ii,jj)=-imag(B); H(jj,ii)=imag(B); H(jj,jj)=real(B);
end
H=H+J*Hphi*J'; H=(H+H')/2;
if isempty(geo.Q); Z=eye(d);
else; [Q,~]=qr(geo.Q); Z=Q(:,size(geo.Q,2)+1:end); end
Ht=Z'*H*Z; Ht=(Ht+Ht')/2;
if isempty(Ht); lambda=inf; else; lambda=min(real(eig(Ht,'vector'))); end
stream=RandStream('mt19937ar','Seed',p.seed+701); v=randn(stream,d,1); v=v/norm(v);
hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec(p.cfg,p.channels,p.smat,W, ...
    real_vector_to_complex_matrix(v,N,U),geo.lambda));
error_hvp=norm(H*v-hv)/max([1,norm(hv),norm(H*v)]);
assert(error_hvp<1e-9,'Dense and analytic Hessian actions differ.');
q_dense=max(geo.kkt/o.epsilon_KKT,max(0,-lambda)/o.epsilon_H);
epsilon2=sol.critical_cone_tolerance;
q_certificate=max(geo.kkt/o.epsilon_KKT,epsilon2/o.epsilon_H);
% q_stop retains the paper's two-term independent tangent-spectrum ratio.
% The complete calibrated stopping condition is recorded separately.
a=struct('crlb',out.crlb,'kkt',geo.kkt,'primal',max([0;c]), ...
    'stationarity',geo.rg,'complementarity',geo.r_lambda, ...
    'dense_lambda_min',lambda,'joint_ratio_dense',q_dense, ...
    'epsilon2',epsilon2,'joint_ratio_certificate',q_certificate, ...
    'joint_ratio',q_dense,'joint_ratio_full',max(q_dense,q_certificate), ...
    'solver_success',sol.success,'fullspace_certificate',sol.certificate.oracle.ambient, ...
    'curvature_lower_bound',sol.curvature_lower_bound, ...
    'geometry_success',geo.success,'sigma_min',geo.sigma_min, ...
    'minimum_active_multiplier',min([inf;geo.lambda(geo.active)]), ...
    'active_count',nnz(geo.active),'rank',size(geo.Q,2), ...
    'hvp_relative_error',error_hvp,'lambda',geo.lambda,'active',geo.active, ...
    'nnls_norm',norm(geo.zeta),'gradient_norm',norm(g), ...
    'constraint_gradient_norm',norm(G,'fro'));
end
