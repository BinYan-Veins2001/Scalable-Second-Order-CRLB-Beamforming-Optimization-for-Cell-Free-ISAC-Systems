function d = simtwc_unified_kkt(p,Wmat)
%SIMTWC_UNIFIED_KKT Recompute the common paper KKT residual for any method.

s=simtwc_config();
[~,g_complex]=compute_crlb_objective_grad(p.cfg,p.smat,Wmat);
[c,Gcomplex]=compute_constraints_grad(p.cfg,p.channels,p.smat,Wmat);
g=complex_matrix_to_real_vector(g_complex);
G=zeros(numel(g),numel(c));
for j=1:numel(c)
    G(:,j)=complex_matrix_to_real_vector(Gcomplex(:,:,j));
end
active=false(size(c));
active(1:p.cfg.sys.U)=c(1:p.cfg.sys.U)>=-s.kkt_sinr_active_tol;
active(p.cfg.sys.U+1:end)=c(p.cfg.sys.U+1:end)>=-s.kkt_power_active_tol;
lambda=zeros(numel(c),1);
if any(active)
    if exist('lsqnonneg','file')==2
        scales=max(vecnorm(G(:,active)),realmin);
        mu=lsqnonneg(G(:,active)./scales,-g);
        lambda(active)=mu./scales.';
    else
        lambda(active)=projected_nnls(G(:,active),-g);
    end
end
lag_grad=g+G*lambda;
d.r_c=max([0;c]);
d.r_g=norm(lag_grad)/max(1,norm(g));
d.r_lambda=max(abs(lambda.*c))/max(1,norm(lambda,inf));
d.r_kkt=max([d.r_c,d.r_g,d.r_lambda]);
d.active=active;
d.active_sinr_count=nnz(active(1:p.cfg.sys.U));
d.active_power_count=nnz(active(p.cfg.sys.U+1:end));
d.lambda=lambda;
d.constraints=c;
end

function x=projected_nnls(A,b)
x=zeros(size(A,2),1);
alpha=1/max(norm(A,2)^2,realmin);
for k=1:2000
    xn=max(0,x-alpha*A'*(A*x-b));
    if norm(xn-x)<=1e-10*max(1,norm(x)); x=xn; return; end
    x=xn;
end
end
