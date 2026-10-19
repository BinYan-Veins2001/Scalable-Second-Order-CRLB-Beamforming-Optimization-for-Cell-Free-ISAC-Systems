function [pg_norm, zproj, step_vec] = projected_grad_mapping_norm(cfg, z, grad_z, N, U, alpha)
%PROJECTED_GRAD_MAPPING_NORM Norm of the projected-gradient mapping.
%
%   G_alpha(z) = (z - Pi_P(z - alpha*grad_z)) / alpha.
%
% This is the stationarity diagnostic for algorithms with hard per-AP
% power projection.

if nargin < 6 || isempty(alpha)
    alpha = 1e-3;
end

ztrial = z - alpha * grad_z;
Wtrial = real_vector_to_complex_matrix(ztrial, N, U);
Wtrial = project_per_ap_power_matrix(cfg, Wtrial);
zproj = complex_matrix_to_real_vector(Wtrial);

step_vec = z - zproj;
pg_norm = norm(step_vec) / max(alpha, realmin);

end
