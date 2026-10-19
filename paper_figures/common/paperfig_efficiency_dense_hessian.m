function [H,diagnostics]=paperfig_efficiency_dense_hessian(p,W,lambda)
%PAPERFIG_EFFICIENCY_DENSE_HESSIAN Independent analytic dense reconstruction.
% Uses the objective chain rule and quadratic constraint Hessians directly.
% No SPSC factors, eigenvectors, inverses, or cached matrices are reused.
% Assembly costs O(U*N^2 + N*U^2 + d^2), not d independent HVP calls.
cfg=p.cfg; N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U; d=2*N*U;
beta=cfg.radio.effective_bandwidth_sq(:); lambda=lambda(:);
matrices={p.smat.Aa,p.smat.Ab,p.smat.Aab};
q=zeros(3,1); J=zeros(d,3);
for j=1:3
    AW=matrices{j}*W;
    q(j)=sum(beta.'.*real(sum(conj(W).*AW,1)));
    J(:,j)=complex_matrix_to_real_vector(2*AW.*beta.');
end
a=q(1); b=q(2); c=q(3); total=a+b; determinant=a*b-c^2;
assert(total>0 && determinant>0 && isfinite(determinant));
alpha=[1/total-b/determinant;1/total-a/determinant;2*c/determinant];
B=[-1/total^2+b^2/determinant^2, ...
   -1/total^2+c^2/determinant^2,-2*b*c/determinant^2; ...
   -1/total^2+c^2/determinant^2, ...
   -1/total^2+a^2/determinant^2,-2*a*c/determinant^2; ...
   -2*b*c/determinant^2,-2*a*c/determinant^2, ...
   2/determinant+4*c^2/determinant^2];
kind=lower(cfg.alm.objective_type);
if strcmp(kind,'crlb_ratio')
    B=(total/determinant)*(B+alpha*alpha');
    alpha=(total/determinant)*alpha;
else
    assert(strcmp(kind,'log_crlb_ratio'),'Unknown objective.');
end
A=zeros(N);
for j=1:3; A=A+alpha(j)*matrices{j}; end
H=zeros(d); channel=p.channels.h_stack;
power=2*kron(lambda(U+1:end)/cfg.power.Pmax,ones(cfg.sys.Nt,1));
common=(channel.*lambda(1:U).')*channel';
for u=1:U
    hu=channel(:,u);
    block=2*beta(u)*A+(2/cfg.comm.noise_power)* ...
        (cfg.comm.Gamma*common-(cfg.comm.Gamma+1)*lambda(u)*(hu*hu'))+diag(power);
    block=(block+block')/2;
    ir=(u-1)*N+(1:N); ii=N*U+ir;
    H(ir,ir)=real(block); H(ir,ii)=-imag(block);
    H(ii,ir)=imag(block); H(ii,ii)=real(block);
end
H=H+J*((B+B')/2)*J'; H=real((H+H')/2);
if nargout>1
    diagnostics=struct('dimension',d,'outer_eigenvalues',eig((B+B')/2), ...
        'assembly','independent analytic block-plus-chain-rule reconstruction');
end
end
