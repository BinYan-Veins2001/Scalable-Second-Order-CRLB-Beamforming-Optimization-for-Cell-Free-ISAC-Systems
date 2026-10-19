function a = benchmark_endpoint_audit(p,z,options)
%BENCHMARK_ENDPOINT_AUDIT Independent NNLS and dense analytic curvature audit.
% This routine does not call either solver's SPSC or HVP implementation.
% Fixed working-set thresholds define the common old/new comparison.
% Optional lambda/active inputs additionally assess a solver's own endpoint.
if nargin<3; options=struct(); end
o=struct('sinr_active_tol',1e-4,'power_active_tol',1e-5, ...
    'epsilon_KKT',1e-5,'epsilon_H',1e-5,'rank_tol',1e-10, ...
    'dense_max_dim',2048,'lambda',[],'active',[]);
names=fieldnames(options);
for j=1:numel(names)
    assert(isfield(o,names{j}),'Unknown audit option: %s',names{j});
    o.(names{j})=options.(names{j});
end
started=tic; cfg=p.cfg; N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U;
z=real(z(:)); d=2*N*U;
assert(numel(z)==d && all(isfinite(z)),'Invalid endpoint.');
W=real_vector_to_complex_matrix(z,N,U);
[f,gc,out]=compute_crlb_objective_grad(cfg,p.smat,W);
g=complex_matrix_to_real_vector(gc);
[c,G]=costar_sosp_constraints(cfg,p.channels,W);
if isempty(o.active)
    active=[c(1:U)>=-o.sinr_active_tol;c(U+1:end)>=-o.power_active_tol];
else
    active=logical(o.active(:)); assert(numel(active)==numel(c));
end
A=G(:,active); scales=max(vecnorm(A),realmin); As=A./scales;
lambda=zeros(size(c)); rankA=0; sigma_raw=inf; Q=zeros(d,0);
if ~isempty(A)
    [Qa,Sa,~]=svd(As,'econ'); sv=diag(Sa);
    rankA=nnz(sv>o.rank_tol); Q=Qa(:,1:rankA);
    sigma_raw=min(svd(A,'econ'));
    if isempty(o.lambda)
        assert(exist('lsqnonneg','file')==2,'Audit requires lsqnonneg.');
        [mu,~,~,exitflag]=lsqnonneg(As,-g);
        assert(exitflag>0,'Independent NNLS did not converge.');
        lambda(active)=mu./scales.';
    end
end
if ~isempty(o.lambda)
    lambda=real(o.lambda(:));
    assert(numel(lambda)==numel(c) && min(lambda)>=-1e-12, ...
        'Invalid solver multipliers.');
end
r=g+G*lambda; primal=max([0;c]);
stationarity=norm(r)/max(1,norm(g));
complementarity=max(abs(lambda.*c))/max(1,norm(lambda,inf));
kkt=max([primal,stationarity,complementarity]);
licq=rankA==nnz(active);
positive=all(lambda(active)>0);
a=struct('dimension',d,'objective',f,'crlb',out.crlb, ...
    'primal',primal,'stationarity',stationarity, ...
    'complementarity',complementarity,'kkt',kkt, ...
    'active',active,'lambda',lambda,'constraints',c, ...
    'active_count',nnz(active),'active_rank',rankA, ...
    'tangent_dimension',d-rankA,'licq',licq, ...
    'all_multipliers_positive',positive, ...
    'min_active_multiplier',min([inf;lambda(active)]), ...
    'lambda_min_tangent',NaN,'lambda_min_ambient',NaN, ...
    'eigenpair_residual',NaN,'hessian_symmetry_error',NaN, ...
    'hessian_norm_bound',NaN,'first_order_pass',kkt<=o.epsilon_KKT, ...
    'second_order_pass',false,'joint_pass',false, ...
    'critical_cone_epsilon2',NaN,'spectral_exact',false, ...
    'spectral_method','not_computed','epsilon_KKT',o.epsilon_KKT, ...
    'epsilon_H',o.epsilon_H,'elapsed_s',NaN, ...
    'note','Dense double-precision eigenvalues are numerical references, not interval bounds.');
if d<=o.dense_max_dim
    [H,symmetry]=analytic_hessian(p,W,lambda);
    a.hessian_symmetry_error=symmetry;
    a.hessian_norm_bound=norm(H,1);
    a.lambda_min_ambient=min(real(eig(H,'vector')));
    if isempty(Q); Z=eye(d);
    else; [Qfull,~]=qr(Q); Z=Qfull(:,rankA+1:end); end
    Ht=Z'*H*Z; Ht=(Ht+Ht')/2;
    if isempty(Ht)
        a.lambda_min_tangent=inf; a.eigenpair_residual=0;
    else
        [V,E]=eig(Ht,'vector'); [minimum,idx]=min(real(E));
        a.lambda_min_tangent=minimum;
        a.eigenpair_residual=norm(Ht*V(:,idx)-minimum*V(:,idx))/ ...
            max(1,norm(Ht,1));
    end
    a.spectral_exact=true; a.spectral_method='dense_independent_analytic';
    a.second_order_pass=a.lambda_min_tangent>=-o.epsilon_H;
    a.joint_pass=a.first_order_pass && a.second_order_pass;
    if licq && positive
        if isempty(A); chi=0;
        else; chi=norm(r)/(min(lambda(active))*sigma_raw); end
        a.critical_cone_epsilon2=max(0,-a.lambda_min_tangent)+ ...
            a.hessian_norm_bound*(2*chi+chi^2);
    end
end
a.elapsed_s=toc(started);
end

function [H,symmetry]=analytic_hessian(p,W,lambda)
cfg=p.cfg; N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U;
d=2*N*U; beta=cfg.radio.effective_bandwidth_sq(:);
matrices={p.smat.Aa,p.smat.Ab,p.smat.Aab};
q=zeros(3,1); J=zeros(d,3);
for j=1:3
    AW=matrices{j}*W;
    q(j)=sum(beta.'.*real(sum(conj(W).*AW,1)));
    J(:,j)=complex_matrix_to_real_vector(2*AW.*beta.');
end
a=q(1); b=q(2); c=q(3); S=a+b; D=a*b-c^2;
assert(S>0 && D>0 && isfinite(S) && isfinite(D), ...
    'Endpoint is outside the smooth positive FIM domain.');
alpha=[1/S-b/D;1/S-a/D;2*c/D];
Bphi=[-1/S^2+b^2/D^2,-1/S^2+c^2/D^2,-2*b*c/D^2; ...
    -1/S^2+c^2/D^2,-1/S^2+a^2/D^2,-2*a*c/D^2; ...
    -2*b*c/D^2,-2*a*c/D^2,2/D+4*c^2/D^2];
objective_type='crlb_ratio';
if isfield(cfg,'alm') && isfield(cfg.alm,'objective_type')
    objective_type=lower(cfg.alm.objective_type);
end
if strcmp(objective_type,'crlb_ratio')
    Bphi=(S/D)*(Bphi+alpha*alpha'); alpha=(S/D)*alpha;
else
    assert(strcmp(objective_type,'log_crlb_ratio'),'Unsupported objective.');
end
Ab=zeros(N);
for j=1:3; Ab=Ab+alpha(j)*matrices{j}; end
H=zeros(d); Hc=p.channels.h_stack;
noise=max(cfg.comm.noise_power,realmin);
Pmax=max(cfg.power.Pmax,realmin);
power=2*kron(lambda(U+1:end)/Pmax,ones(cfg.sys.Nt,1));
for u=1:U
    weights=cfg.comm.Gamma*lambda(1:U); weights(u)=-lambda(u);
    Du=2*beta(u)*Ab+(2/noise)*(Hc.*weights.')*Hc'+diag(power);
    Du=(Du+Du')/2;
    ii=(u-1)*N+(1:N); jj=N*U+ii;
    H(ii,ii)=real(Du); H(ii,jj)=-imag(Du);
    H(jj,ii)=imag(Du); H(jj,jj)=real(Du);
end
H=H+J*((Bphi+Bphi')/2)*J';
symmetry=norm(H-H','fro')/max(1,norm(H,'fro'));
H=real((H+H')/2);
end
