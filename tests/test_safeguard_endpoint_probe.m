function report=test_safeguard_endpoint_probe()
% Full-space certificate with one slightly inactive SINR constraint.
% Fine geometry excludes that constraint. Its small nonnegative multiplier
% nevertheless gives a valid complete approximate KKT / ambient PSD pair.
cfg=default_config(); cfg.sys.M=1; cfg.sys.Nt=2; cfg.sys.Nr=1; cfg.sys.U=1;
cfg.power.Pmax=1; cfg.comm.noise_power=1; cfg.comm.Gamma=0.1;
cfg.radio.effective_bandwidth_sq=1; cfg.radio.c=1; cfg.radio.frame_duration=1;
cfg.sensing.N0=8*pi^2; cfg.alm.objective_type='log_crlb_ratio'; cfg.alm.det_floor=1e-30;
cfg.costar.paper=struct(); channels=struct('h_stack',[1;0],'h',{{[1;0]}});
A=diag([1,1+4e-5]); smat=struct('Aa',A,'Ab',A,'Aab',zeros(2));
W=[sqrt(0.1+6e-6);sqrt(0.9-6e-6)];
o=struct('max_iter',0,'epsilon_KKT',1e-5,'epsilon_H',1e-5,'nu0',1e-5, ...
    'ambient_endpoint_probe',false,'feas_tol',1e-10);
[base,~]=solve_costar(cfg,channels,smat,{W},o);
assert(~base.success && base.kkt_residual>o.epsilon_KKT);
o.ambient_endpoint_probe=true;
[sol,H]=solve_costar(cfg,channels,smat,{W},o);
assert(sol.success && sol.kkt_residual<=o.epsilon_KKT);
assert(sol.certificate.oracle.ambient && sol.certificate.oracle.conditional_psd);
assert(sol.certificate.oracle_is_at_final_point);
assert(sol.lambda(1)>0 && all(sol.active) && sol.complementarity<1e-8);
assert(norm(sol.Wmat-W)==0 && sol.endpoint_geometry_nu>sol.nu);
assert(strcmp(H.rtr.step_type{end},'ambient_endpoint_certificate'));
assert(abs(H.rtr.kkt(end)-sol.kkt_residual)<1e-12);
% Independently assembled physical Hessian, using the same multiplier.
gA=real(W'*A*W); M=-2*A/gA+4*(A*W)*(A*W)'/gA^2 ...
    -2*sol.lambda(1)*diag([1,0])+2*sol.lambda(2)*eye(2);
Mi=-2*A/gA-2*sol.lambda(1)*diag([1,0])+2*sol.lambda(2)*eye(2);
assert(min(eig(blkdiag(M,Mi)))>=-o.epsilon_H);
report=struct('passed',true,'fine_kkt',base.kkt_residual, ...
    'probe_kkt',sol.kkt_residual,'primal',sol.primal_residual, ...
    'ambient_min_eigenvalue',min(eig(blkdiag(M,Mi))));
disp(report);
end
