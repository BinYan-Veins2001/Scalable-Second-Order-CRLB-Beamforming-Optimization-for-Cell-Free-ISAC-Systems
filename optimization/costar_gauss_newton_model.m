function [J,Bpsd]=costar_gauss_newton_model(cfg,smat,W,out,eigen_floor)
%COSTAR_GAUSS_NEWTON_MODEL PSD generalized-GN ablation, not a curvature oracle.
N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U;
J=zeros(2*N*U,3); beta2=cfg.radio.effective_bandwidth_sq(:);
matrices={smat.Aa,smat.Ab,smat.Aab};
for j=1:3
    J(:,j)=complex_matrix_to_real_vector(2*(matrices{j}*W).*beta2.');
end
a=max(real(out.ga),realmin); b=max(real(out.gb),realmin); c=real(out.gab);
s=max(a+b,realmin); d=max(a*b-c^2,realmin);
g=[1/s-b/d;1/s-a/d;2*c/d];
B=[-1/s^2+b^2/d^2,-1/s^2+c^2/d^2,-2*b*c/d^2; ...
    -1/s^2+c^2/d^2,-1/s^2+a^2/d^2,-2*a*c/d^2; ...
    -2*b*c/d^2,-2*a*c/d^2,2/d+4*c^2/d^2];
if strcmpi(out.objective_type,'crlb_ratio'); B=out.ratio*(B+g*g'); end
[V,e]=eig((B+B')/2,'vector');
Bpsd=real(V*diag(max(real(e),eigen_floor*max([1;abs(e)])))*V');
Bpsd=(Bpsd+Bpsd')/2;
end
