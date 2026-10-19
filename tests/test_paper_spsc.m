function report=test_paper_spsc()
% Small exact references exercise ambient, tangent, recovery and shift paths.
cfg=struct(); cfg.sys=struct('M',1,'Nt',2,'U',1);
cfg.radio.effective_bandwidth_sq=1;
cfg.comm=struct('Gamma',1,'noise_power',1); cfg.power.Pmax=1;
cfg.alm=struct('objective_type','log_crlb_ratio','det_floor',1e-30);
channels.h_stack=[1;1i];
smat=struct('Aa',[1,1;1,2],'Ab',[1,-1;-1,2],'Aab',zeros(2));
W=[1;0]; lambda=[0;1]; epsilon_H=.1;
H=zeros(4);
for j=1:4
    e=zeros(4,1); e(j)=1;
    H(:,j)=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
        cfg,channels,smat,W,real_vector_to_complex_matrix(e,2,1),lambda));
end
assert(norm(H-diag([4,6,0,-2]),'fro')<1e-12);
[ambient,~]=costar_sosp_structured_psd(cfg,channels,smat,W,lambda,epsilon_H);
assert(~ambient.conditional_psd && strcmp(ambient.status,'structured_ambient_indefinite'));
Q=[eye(4,1),[0;0;0;1]];
[tangent,parts]=costar_sosp_structured_psd(cfg,channels,smat,W,lambda,epsilon_H,Q);
assert(tangent.conditional_psd && strcmp(tangent.mode,'tangent'));
assert(tangent.ambient_negative_inertia==1 && tangent.tangent_negative_inertia==0);
assert(tangent.shift>epsilon_H/2 && tangent.shift<=3*epsilon_H/4);
assert(parts.beta>=norm(H,2)-1e-12);
direction=costar_sosp_structured_psd(cfg,channels,smat,W,lambda,epsilon_H,eye(4,1));
assert(strcmp(direction.status,'negative_curvature'));
assert(abs(norm(direction.negative_direction)-1)<1e-12);
assert(direction.rayleigh<-epsilon_H/2 && direction.tangent_residual<1e-10);
assert(abs(direction.rayleigh-direction.negative_direction'*H*direction.negative_direction)<1e-12);
full=costar_sosp_structured_psd(cfg,channels,smat,W,[0;3],epsilon_H,zeros(4,0));
assert(full.conditional_psd && strcmp(full.mode,'ambient'));
empty=costar_sosp_structured_psd(cfg,channels,smat,W,lambda,epsilon_H,eye(4));
assert(empty.conditional_psd && isinf(empty.eigenvalue_lower_bound));
% The preferred sigma makes an unshifted block eigenvalue exactly zero.
alternative=costar_sosp_structured_psd(cfg,channels,smat,W,[0;.9625],epsilon_H,eye(4,1));
assert(alternative.shift_attempts>1 && alternative.shift<3*epsilon_H/4);
assert(strcmp(alternative.status,'negative_curvature'));
previous=rng; cleanup=onCleanup(@()rng(previous)); rng(8061,'twister'); %#ok<NASGU>
comparisons=0;
for rankQ=0:3
    for trial=1:4
        [Q,~]=qr(randn(4,rankQ),0);
        info=costar_sosp_structured_psd(cfg,channels,smat,W,lambda,epsilon_H,Q);
        Z=null(Q'); eigenvalues=eig(Z'*(H+info.shift*eye(4))*Z);
        expected=nnz(real(eigenvalues)<0);
        assert(info.tangent_negative_inertia==expected);
        assert(info.conditional_psd==(expected==0));
        if expected>0
            assert(strcmp(info.status,'negative_curvature'));
            assert(info.rayleigh<-epsilon_H/2 && norm(Q'*info.negative_direction)<1e-9);
        end
        comparisons=comparisons+1;
    end
end
report=struct('passed',true,'dense_inertia_comparisons',comparisons, ...
    'ambient_first',true,'tangent_certificate',true, ...
    'deterministic_negative_recovery',true,'singular_shift_reselected',true);
disp(report);
end
