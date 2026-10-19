function [lambda_min_est, info] = estimate_min_curvature_lanczos(Hvfun, dim, opts)
%ESTIMATE_MIN_CURVATURE_LANCZOS Randomized Lanczos estimate of min eigenvalue.
%
% This is a diagnostic routine. It is optional and not needed for the core
% trust-region step because tCG already detects negative curvature.

if nargin < 3
    opts = struct();
end
if ~isfield(opts,'max_iter'); opts.max_iter = min(30, dim); end
if ~isfield(opts,'tol'); opts.tol = 1e-8; end
if ~isfield(opts,'seed'); opts.seed = []; end

if ~isempty(opts.seed)
    rng(opts.seed);
end

kmax = min(opts.max_iter, dim);
Q = zeros(dim, kmax);
alpha = zeros(kmax,1);
beta = zeros(kmax,1);

q = randn(dim,1);
q = q / max(norm(q), realmin);
q_prev = zeros(dim,1);
beta_prev = 0;
actual_k = 0;

for k = 1:kmax
    Q(:,k) = q;
    z = Hvfun(q);
    if k > 1
        z = z - beta_prev * q_prev;
    end
    alpha(k) = q' * z;
    z = z - alpha(k) * q;

    % Full reorthogonalization for robustness in small/medium dimensions.
    for j = 1:k
        z = z - Q(:,j) * (Q(:,j)' * z);
    end

    beta(k) = norm(z);
    actual_k = k;

    if beta(k) <= opts.tol
        break;
    end

    q_prev = q;
    q = z / beta(k);
    beta_prev = beta(k);
end

T = diag(alpha(1:actual_k));
if actual_k >= 2
    off = beta(1:actual_k-1);
    T = T + diag(off,1) + diag(off,-1);
end

[V_T,D_T] = eig((T+T')/2);
eigs_T = diag(D_T);
[lambda_min_est,min_idx] = min(eigs_T);

info.iter = actual_k;
info.ritz_values = eigs_T;
info.T = T;
info.ritz_vector = Q(:,1:actual_k) * V_T(:,min_idx);
info.ritz_vector = info.ritz_vector / max(norm(info.ritz_vector),realmin);

end
