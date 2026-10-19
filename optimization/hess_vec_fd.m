function [Hv, info] = hess_vec_fd(grad_fun, z, v, eps_base)
%HESS_VEC_FD Central-difference Hessian-vector product.
%
%   Hv ~= [g(z+eps*v) - g(z-eps*v)] / (2*eps).
%
% grad_fun must be a function handle returning a real gradient vector:
%   g = grad_fun(z_query)
% or
%   [~, g] = value_grad_fun(z_query), if wrapped by the caller.

if nargin < 4 || isempty(eps_base)
    eps_base = 1e-5;
end

vnorm = norm(v);
if vnorm <= 0
    Hv = zeros(size(v));
    info.eps = 0;
    info.method = 'zero_direction';
    return;
end

eps_fd = eps_base * (1 + norm(z)) / max(vnorm, 1e-12);

try
    gp = grad_fun(z + eps_fd * v);
    gm = grad_fun(z - eps_fd * v);
    Hv = (gp - gm) / (2 * eps_fd);
    info.method = 'central';
catch
    % Fallback to forward difference if the negative perturbation fails.
    g0 = grad_fun(z);
    gp = grad_fun(z + eps_fd * v);
    Hv = (gp - g0) / eps_fd;
    info.method = 'forward_fallback';
end

if any(~isfinite(Hv))
    g0 = grad_fun(z);
    gp = grad_fun(z + eps_fd * v);
    Hv = (gp - g0) / eps_fd;
    info.method = 'forward_nonfinite_fallback';
end

info.eps = eps_fd;

end
