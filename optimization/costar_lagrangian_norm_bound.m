function L = costar_lagrangian_norm_bound(cfg,channels,smat,z,lambda)
%COSTAR_LAGRANGIAN_NORM_BOUND Structured upper bound on ||Hess L||_2.
% Uses Frobenius bounds on sensing actions and rank-one channel norms.
N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U;
W=real_vector_to_complex_matrix(z,N,U);
[~,info]=exact_lagrangian_hess_vec(cfg,channels,smat,W,zeros(N,U),lambda);
alpha=info.objective.alpha; B=info.objective.outer_hessian;
beta=cfg.radio.effective_bandwidth_sq(:).';
A={smat.Aa,smat.Ab,smat.Aab};
jbound=0; direct=0;
for r=1:3
    Jcomplex=2*(A{r}*W).*beta;
    jbound=jbound+norm(Jcomplex,'fro')^2;
    direct=direct+2*max(beta)*abs(alpha(r))*norm(A{r},'fro');
end
L=direct+norm(B,2)*jbound;
hnorm=sum(abs(channels.h_stack).^2,1).';
L=L+2/max(cfg.comm.noise_power,realmin)*max(1,cfg.comm.Gamma)* ...
    sum(abs(lambda(1:U)).*hnorm);
L=L+2/max(cfg.power.Pmax,realmin)*max(abs(lambda(U+1:end)));
L=real(L)*(1+1e-12);
end
