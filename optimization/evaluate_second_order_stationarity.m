function diag = evaluate_second_order_stationarity(cfg, channels, smat, z, opts)
%EVALUATE_SECOND_ORDER_STATIONARITY KKT and weak second-order diagnostics.
%
% For c(z) <= 0, multipliers are estimated by nonnegative least squares on
% the numerically active constraints. Curvature is measured on the nullspace
% of all active constraint gradients, which is the standard weak critical
% subspace used in practical weak second-order necessary-condition checks.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(opts, numel(z));

N = cfg.sys.M * cfg.sys.Nt;
U = cfg.sys.U;
z = real(z(:));
[grad_f, c, Gz] = objective_constraint_derivatives(cfg, channels, smat, z, N, U);

active = c >= -opts.active_tol;
lambda = zeros(numel(c), 1);
if any(active)
    lambda(active) = solve_nonnegative_least_squares(Gz(:,active), -grad_f, opts);
end

lag_grad = grad_f + Gz * lambda;
primal_residual = max([0; c]);
stationarity_residual = norm(lag_grad) / max(1, norm(grad_f));
complementarity_residual = max(abs(lambda .* c)) / max(1, norm(lambda, inf));
kkt_residual = max([primal_residual, stationarity_residual, ...
    complementarity_residual]);

Q = zeros(numel(z), 0);
if any(active)
    Q = orth(Gz(:,active), opts.rank_tol);
end
tangent_dim = numel(z) - size(Q, 2);

if tangent_dim <= 0
    lambda_min = inf;
    negative_curvature = 0;
    normalized_negative_curvature = 0;
    lanczos_info = struct('iter', 0, 'ritz_values', [], 'T', []);
else
    project_tangent = @(v) v - Q * (Q' * v);
    Hv = @(v) project_tangent(lagrangian_hess_vec( ...
        cfg, channels, smat, z, project_tangent(v), lambda, N, U, opts));
    lanczos_opts.max_iter = min(opts.lanczos_max_iter, numel(z));
    lanczos_opts.tol = opts.lanczos_tol;
    lanczos_opts.seed = opts.seed;
    [lambda_min, lanczos_info] = estimate_min_curvature_lanczos( ...
        Hv, numel(z), lanczos_opts);
    negative_curvature = max(0, -lambda_min);
    spectral_scale = max([1; abs(lanczos_info.ritz_values(:))]);
    normalized_negative_curvature = negative_curvature / spectral_scale;
end

diag.primal_residual = primal_residual;
diag.stationarity_residual = stationarity_residual;
diag.complementarity_residual = complementarity_residual;
diag.kkt_residual = kkt_residual;
diag.lambda_min_weak_critical = lambda_min;
diag.negative_curvature = negative_curvature;
diag.normalized_negative_curvature = normalized_negative_curvature;
diag.sosp_residual = max(kkt_residual, normalized_negative_curvature);
diag.active_count = nnz(active);
diag.tangent_dim = tangent_dim;
diag.lambda = lambda;
diag.constraints = c;
diag.lanczos = lanczos_info;
diag.hessian_model = opts.hessian_model;
end

function Hv = lagrangian_hess_vec(cfg, channels, smat, z, v, lambda, N, U, opts)
if any(strcmpi(opts.hessian_model, ...
        {'exact_lagrangian','lagrangian_exact','exact_matrix_free'}))
    W = real_vector_to_complex_matrix(z,N,U);
    V = real_vector_to_complex_matrix(real(v(:)),N,U);
    Hcomplex = exact_lagrangian_hess_vec( ...
        cfg,channels,smat,W,V,lambda);
    Hv = real(complex_matrix_to_real_vector(Hcomplex));
    return;
end
vnorm = norm(v);
if vnorm <= realmin
    Hv = zeros(size(v));
    return;
end
step = opts.fd_rel_step * max(1, norm(z)) / vnorm;
gp = lagrangian_gradient(cfg, channels, smat, z + step*v, lambda, N, U);
gm = lagrangian_gradient(cfg, channels, smat, z - step*v, lambda, N, U);
Hv = real((gp - gm) / (2*step));
end

function grad_L = lagrangian_gradient(cfg, channels, smat, z, lambda, N, U)
[grad_f, ~, Gz] = objective_constraint_derivatives( ...
    cfg, channels, smat, z, N, U);
grad_L = grad_f + Gz * lambda;
end

function [grad_f, c, Gz] = objective_constraint_derivatives( ...
    cfg, channels, smat, z, N, U)
Wmat = real_vector_to_complex_matrix(z, N, U);
[~, grad_f_complex] = compute_crlb_objective_grad(cfg, smat, Wmat);
[c, Gcomplex] = compute_constraints_grad(cfg, channels, smat, Wmat);
grad_f = complex_matrix_to_real_vector(grad_f_complex);
Gz = zeros(numel(z), numel(c));
for j = 1:numel(c)
    Gz(:,j) = complex_matrix_to_real_vector(Gcomplex(:,:,j));
end
end

function x = solve_nonnegative_least_squares(A, b, opts)
if isempty(A)
    x = zeros(0,1);
    return;
end
if exist('lsqnonneg', 'file') == 2
    x = lsqnonneg(A, b);
    return;
end

x = zeros(size(A,2),1);
L = norm(A,2)^2;
alpha = 1 / max(L, realmin);
for k = 1:opts.nnls_max_iter
    xnew = max(0, x - alpha * (A' * (A*x - b)));
    if norm(xnew-x) <= opts.nnls_tol * max(1,norm(x))
        x = xnew;
        break;
    end
    x = xnew;
end
end

function opts = fill_defaults(opts, dim)
if ~isfield(opts,'active_tol'); opts.active_tol = 1e-3; end
if ~isfield(opts,'rank_tol'); opts.rank_tol = 1e-10; end
if ~isfield(opts,'hessian_model'); opts.hessian_model = 'exact_lagrangian'; end
if ~isfield(opts,'fd_rel_step'); opts.fd_rel_step = 1e-5; end
if ~isfield(opts,'lanczos_max_iter'); opts.lanczos_max_iter = min(20,dim); end
if ~isfield(opts,'lanczos_tol'); opts.lanczos_tol = 1e-8; end
if ~isfield(opts,'seed'); opts.seed = 2026; end
if ~isfield(opts,'nnls_max_iter'); opts.nnls_max_iter = 500; end
if ~isfield(opts,'nnls_tol'); opts.nnls_tol = 1e-8; end
end
