function report=test_paperfig_efficiency()
%TEST_PAPERFIG_EFFICIENCY Audit dense reconstruction and all SPSC outcomes.
cfg=struct(); cfg.sys=struct('M',1,'Nt',2,'U',1);
cfg.radio.effective_bandwidth_sq=1; cfg.comm=struct('Gamma',1,'noise_power',1);
cfg.power.Pmax=1; cfg.alm=struct('objective_type','log_crlb_ratio','det_floor',1e-30);
p=struct('cfg',cfg,'channels',struct('h_stack',[1;1i]), ...
    'smat',struct('Aa',[1 1;1 2],'Ab',[1 -1;-1 2],'Aab',zeros(2)));
W=[1;0]; epsH=.1;
normals={zeros(4,0),[eye(4,1),[0;0;0;1]],eye(4,1),eye(4)};
multiplier={[0;3],[0;1],[0;1],[0;1]};
expected={'ambient','tangent','negative_curvature','empty_tangent'};
for trial=1:numel(normals)
    lambda=multiplier{trial}; Q=normals{trial};
    H=paperfig_efficiency_dense_hessian(p,W,lambda);
    for j=1:4
        v=zeros(4,1); v(j)=1;
        exact=complex_matrix_to_real_vector(exact_lagrangian_hess_vec(cfg,p.channels,p.smat,W, ...
            real_vector_to_complex_matrix(v,2,1),lambda));
        assert(norm(H*v-exact)<1e-12);
    end
    s=costar_sosp_structured_psd(cfg,p.channels,p.smat,W,lambda,epsH,Q);
    assert(strcmp(s.mode,expected{trial}));
    sigma=s.shift; if isnan(sigma); sigma=3*epsH/4; end
    e=paperfig_efficiency_dense_evd(p,W,lambda,Q,sigma);
    assert(s.tangent_negative_inertia==e.negative_inertia && s.conditional_psd==e.certified);
    if ~isempty(s.negative_direction)
        assert(norm(Q'*s.negative_direction)<1e-10);
        assert(s.negative_direction'*H*s.negative_direction<-epsH/2);
    end
end
report=struct('passed',true,'dense_hvp_columns',16,'branches',{expected});
disp(report);
end
