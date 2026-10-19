function [a,op] = simtwc_theorem3_audit(p,z,o,spectral)
%SIMTWC_THEOREM3_AUDIT Independent KKT and dense tangent-spectrum audit.
% Dense matrices are formed ONLY for offline validation, not by COSTAR.
% Eigenvalues below are numerical references, not interval certificates.
if nargin<4; spectral=true; end
cfg=p.cfg; N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U; d=2*N*U;
W=real_vector_to_complex_matrix(z,N,U);
[f,gc,out]=compute_crlb_objective_grad(cfg,p.smat,W);
g=complex_matrix_to_real_vector(gc);
[c,G]=costar_sosp_constraints(cfg,p.channels,W);
active=[c(1:U)>=-cfg.costar.sosp.sinr_active_tol; ...
    c(U+1:end)>=-cfg.costar.sosp.power_active_tol];
A=G(:,active); scale=max(vecnorm(A),realmin);
As=A./scale;
rank_tol=1e-10;
if isfield(cfg.costar.sosp,'rank_tol'); rank_tol=cfg.costar.sosp.rank_tol; end
lambda=zeros(size(c));
if isempty(A)
    Q=zeros(d,0); rankA=0; sigma=inf; scaled_sigma=inf;
    ls_defect=0; positive=true; nnls_exitflag=1;
else
    [Qa,Sa,~]=svd(As,'econ'); sv=diag(Sa);
    rankA=nnz(sv>rank_tol); Q=Qa(:,1:rankA);
    [mu,~,~,nnls_exitflag]=lsqnonneg(As,-g); lambda(active)=mu./scale.';
    mu_ls=-(As\g);
    ls_defect=norm(As*(mu-mu_ls))/max(1,norm(g));
    sigma=min(svd(A,'econ')); scaled_sigma=min(sv);
    positive=all(lambda(active)>0);
end
r=g+G*lambda; gt=g-Q*(Q'*g);
rc=max([0;c]); rg=norm(r)/max(1,norm(g));
rl=max(abs(lambda.*c))/max(1,norm(lambda,inf));
face=max([0;abs(c(active))]);
licq=rankA==nnz(active) && scaled_sigma>=rank_tol;
if any(active); lambda0=min(lambda(active)); else; lambda0=inf; end
a=struct('objective',f,'crlb',out.crlb,'kkt',max([rc,rg,rl]), ...
    'primal',rc,'stationarity',rg,'complementarity',rl, ...
    'tangent_gradient',norm(gt),'active',active,'lambda',lambda,'constraints',c, ...
    'licq',licq,'positive_multipliers',positive,'lambda_min_active',lambda0, ...
    'sigma_min_raw',sigma,'sigma_min_scaled',scaled_sigma, ...
    'nnls_ls_defect',ls_defect,'nnls_exitflag',nnls_exitflag,'face_residual',face, ...
    'inactive_slack',min([inf;-c(~active)]), ...
    'tangent_dimension',d-rankA,'dense_lambda_min',NaN, ...
    'dense_ambient_minimum',NaN,'ambient_joint_pass',false, ...
    'hvp_relative_error',NaN,'beta_L_bound',NaN, ...
    'epsilon2_requested',NaN,'epsilon2_residual',NaN, ...
    'cone_sample_minimum',NaN,'cone_sample_violation',NaN, ...
    'cone_sample_count',0,'cone_bound_observed',false,'joint_ratio',NaN);
op=struct();
if ~spectral; return; end
assert(d<=o.dense_max_dim, ...
    'Dense validation dimension %d exceeds dense_max_dim=%d.',d,o.dense_max_dim);

% Realification uses [real(vec(W)); imag(vec(W))], not interleaved blocks.
[~,hvinfo]=exact_lagrangian_hess_vec(cfg,p.channels,p.smat,W,zeros(N,U),lambda);
alpha=hvinfo.objective.alpha; Bphi=hvinfo.objective.outer_hessian;
matrices={p.smat.Aa,p.smat.Ab,p.smat.Aab}; beta=cfg.radio.effective_bandwidth_sq(:);
Ab=zeros(N); J=zeros(d,3);
for j=1:3
    Ab=Ab+alpha(j)*matrices{j};
    J(:,j)=complex_matrix_to_real_vector(2*(matrices{j}*W).*beta.');
end
H=zeros(d); Hc=p.channels.h_stack;
power=2*kron(lambda(U+1:end)/cfg.power.Pmax,ones(cfg.sys.Nt,1));
for i=1:U
    weights=cfg.comm.Gamma*lambda(1:U); weights(i)=-lambda(i);
    Di=2*beta(i)*Ab+(2/cfg.comm.noise_power)*(Hc.*weights.')*Hc'+diag(power);
    Di=(Di+Di')/2; ii=(i-1)*N+(1:N); jj=N*U+ii;
    H(ii,ii)=real(Di); H(ii,jj)=-imag(Di);
    H(jj,ii)=imag(Di); H(jj,jj)=real(Di);
end
H=H+J*Bphi*J'; H=(H+H')/2;
a.dense_ambient_minimum=min(real(eig(H,'vector')));
a.ambient_joint_pass=a.kkt<=o.epsilon_KKT && ...
    a.dense_ambient_minimum>=-o.epsilon_H;
if isempty(Q)
    Z=eye(d);
else
    [fullQ,~]=qr(Q); Z=fullQ(:,rankA+1:end);
end
Ht=Z'*H*Z; Ht=(Ht+Ht')/2;
if isempty(Ht)
    a.dense_lambda_min=inf; vmin=zeros(d,1);
else
    [V,E]=eig(Ht,'vector'); [a.dense_lambda_min,imin]=min(real(E));
    vmin=Z*V(:,imin);
end
a.beta_L_bound=norm(H,1);
saved_rng=rng; cleanup=onCleanup(@() rng(saved_rng)); %#ok<NASGU>
rng(p.seed+913,'twister'); probe=randn(d,3); errors=zeros(1,3);
for j=1:3
    v=probe(:,j)/norm(probe(:,j));
    hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
        cfg,p.channels,p.smat,W,real_vector_to_complex_matrix(v,N,U),lambda));
    errors(j)=norm(H*v-hv)/max([1,norm(hv),norm(H*v)]);
end
a.hvp_relative_error=max(errors);
assert(a.hvp_relative_error<1e-9,'Dense Hessian and exact HVP disagree.');
a.joint_ratio=max(a.kkt/o.epsilon_KKT,max(0,-a.dense_lambda_min)/o.epsilon_H);
if licq && positive
    if isempty(A)
        chi=0; chi_requested=0;
    else
        chi=norm(r)/(lambda0*sigma);
        chi_requested=max(1,norm(g))*o.epsilon_KKT/(lambda0*sigma);
    end
    a.epsilon2_requested=o.epsilon_H+a.beta_L_bound*(2*chi_requested+chi_requested^2);
    a.epsilon2_residual=o.epsilon_H+a.beta_L_bound*(2*chi+chi^2);
    % Sample feasible linearized critical directions; not a global cone search.
    [a.cone_sample_minimum,a.cone_sample_count]=cone_samples( ...
        H,g,As,Z,vmin,o.cone_samples);
    a.cone_sample_violation=max(0,-a.cone_sample_minimum);
    a.cone_bound_observed=a.cone_sample_violation<= ...
        a.epsilon2_residual+1e-9*max(1,a.beta_L_bound);
end
op=struct('H',H,'Ht',Ht,'Z',Z,'Q',Q,'g',g,'r',r);
end

function [minimum,count]=cone_samples(H,g,As,Z,vmin,n)
minimum=inf; count=0;
if isempty(Z); return; end
T=Z*randn(size(Z,2),n); T(:,1)=vmin;
if ~isempty(As)
    [Ua,Sa,Va]=svd(As,'econ');
    normal_map=-Ua*(Va'./diag(Sa));
end
for j=1:n
    t=T(:,j); if norm(t)==0; continue; end
    t=t/norm(t); if g'*t>0; t=-t; end
    v=t;
    if ~isempty(As)
        aa=rand(size(As,2),1); aa=aa/max(norm(aa),realmin);
        normal=normal_map*aa;
        gain=real(g'*normal);
        if gain>0
            amount=rand*max(0,-real(g'*t))/gain;
        else
            amount=rand;
        end
        v=t+amount*normal;
    end
    v=v/max(norm(v),realmin);
    scale=max(1,norm(g));
    if max([0;As'*v])<=1e-9 && g'*v<=1e-9*scale
        minimum=min(minimum,real(v'*H*v)); count=count+1;
    end
end
if count==0; minimum=NaN; end
end
