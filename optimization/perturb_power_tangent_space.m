function [z_pert,info] = perturb_power_tangent_space( ...
    cfg,z,N,U,stream)
%PERTURB_POWER_TANGENT_SPACE Random perturbation for projected power balls.
%
% Active AP blocks use the real tangent space
%   Re <W_m,D_m> = 0.
% Interior blocks use their full ambient space. The resulting point is
% projected back onto the product of per-AP power balls for exact feasibility.

W = real_vector_to_complex_matrix(z,N,U);
D = (randn(stream,N,U)+1j*randn(stream,N,U))/sqrt(2);
Nt = cfg.sys.Nt;
active_blocks = 0;
max_tangent_residual = 0;

for m = 1:cfg.sys.M
    rows = (m-1)*Nt+(1:Nt);
    Wm = W(rows,:);
    Dm = D(rows,:);
    power_m = norm(Wm,'fro')^2;
    is_active = power_m >= cfg.power.Pmax*(1-cfg.ppfo.power_active_tol);
    if is_active && power_m > realmin
        radial_coeff = real(sum(conj(Wm(:)).*Dm(:)))/power_m;
        Dm = Dm-radial_coeff*Wm;
        tangent_inner = abs(real(sum(conj(Wm(:)).*Dm(:))));
        max_tangent_residual = max(max_tangent_residual, ...
            tangent_inner/max(1,norm(Wm,'fro')*norm(Dm,'fro')));
        D(rows,:) = Dm;
        active_blocks = active_blocks+1;
    end
end

d_eff = max(1,2*numel(W)-active_blocks);
direction_norm = norm(D,'fro');
if direction_norm <= realmin
    error('PPFO-ALM failed to sample a nonzero tangent perturbation.');
end
D = D/direction_norm;

requested_radius = max(cfg.ppfo.min_perturb_radius, ...
    cfg.ppfo.perturb_radius_relative*max(1,norm(z)));
sampled_radius = requested_radius*rand(stream)^(1/d_eff);
W_raw = W+sampled_radius*D;
W_pert = project_per_ap_power_matrix(cfg,W_raw);
z_pert = complex_matrix_to_real_vector(W_pert);

info.requested_radius = requested_radius;
info.sampled_radius = sampled_radius;
info.actual_step_norm = norm(z_pert-z);
info.active_blocks = active_blocks;
info.effective_dimension = d_eff;
info.max_tangent_residual = max_tangent_residual;
end
