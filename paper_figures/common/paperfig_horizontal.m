function a=paperfig_horizontal(p,W,o,geo)
%PAPERFIG_HORIZONTAL Offline finite-Lanczos horizontal Ritz estimate.
% No random curvature-search routine is called by the optimizer here.
N=size(W,1); U=size(W,2); z=complex_matrix_to_real_vector(W); d=numel(z);
[~,gc]=compute_crlb_objective_grad(p.cfg,p.smat,W);
g=complex_matrix_to_real_vector(gc); [c,G]=costar_sosp_constraints(p.cfg,p.channels,W);
if nargin<4 || isempty(geo)
    active=[c(1:U)>=-1e-4;c(U+1:end)>=-1e-5];
    A=G(:,active); lambda=zeros(size(c));
    if isempty(A); Q=zeros(d,0);
    else
        scales=max(vecnorm(A),realmin); As=A./scales;
        Q=orth(As,1e-11); mu=lsqnonneg(As,-g); lambda(active)=mu./scales.';
    end
    zeta=g+G*lambda;
    kkt=max([0;c;norm(zeta)/max(1,norm(g));norm(lambda.*c,inf)/max(1,norm(lambda,inf))]);
else
    active=geo.active; Q=geo.Q; lambda=geo.lambda; kkt=geo.kkt;
end
phase=zeros(d,U);
for u=1:U
    V=zeros(N,U); V(:,u)=1i*W(:,u);
    phase(:,u)=complex_matrix_to_real_vector(V);
end
Q=orth([Q phase],1e-11); nh=d-size(Q,2); project=@(v)v-Q*(Q'*v);
if nh<=0
    a=struct('lambda_min_horizontal',inf,'kkt',kkt,'lanczos_iter',0,'horizontal_dim',0); return;
end
stream=RandStream('mt19937ar','Seed',p.seed+7006);
q=project(randn(stream,d,1)); q=q/norm(q);
m=min(o.lanczos_max_iter,nh); V=zeros(d,m); alpha=zeros(m,1); beta=zeros(m,1);
prev=zeros(d,1); b=0;
for j=1:m
    V(:,j)=q;
    hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
        p.cfg,p.channels,p.smat,W,real_vector_to_complex_matrix(q,N,U),lambda));
    r=project(hv)-b*prev; alpha(j)=q'*r; r=r-alpha(j)*q;
    r=project(r-V(:,1:j)*(V(:,1:j)'*r)); beta(j)=norm(r);
    if beta(j)<=o.lanczos_tol; break; end
    prev=q; q=r/beta(j); b=beta(j);
end
T=diag(alpha(1:j));
if j>1; T=T+diag(beta(1:j-1),1)+diag(beta(1:j-1),-1); end
[B,E]=eig((T+T')/2,'vector'); [value,index]=min(real(E));
v=project(V(:,1:j)*B(:,index)); v=v/norm(v);
hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
    p.cfg,p.channels,p.smat,W,real_vector_to_complex_matrix(v,N,U),lambda));
a=struct('lambda_min_horizontal',real(v'*hv),'ritz_matrix_min',value, ...
    'kkt',kkt,'lanczos_iter',j,'horizontal_dim',nh, ...
    'horizontal_residual',norm(Q'*v),'ritz_residual',norm(project(hv)-value*v), ...
    'lambda',lambda,'active',active);
end
