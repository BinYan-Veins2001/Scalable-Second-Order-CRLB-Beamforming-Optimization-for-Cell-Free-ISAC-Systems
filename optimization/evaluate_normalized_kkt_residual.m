function diag = evaluate_normalized_kkt_residual(cfg,channels,smat,z,opts)
%EVALUATE_NORMALIZED_KKT_RESIDUAL Common active-set KKT diagnostic.
% Uses the same first-order residual definition as
% evaluate_second_order_stationarity, without the curvature calculation.

if nargin < 5
    opts = struct();
end
if ~isfield(opts,'active_tol'); opts.active_tol = 1e-3; end

N = cfg.sys.M*cfg.sys.Nt;
U = cfg.sys.U;
z = real(z(:));
Wmat = real_vector_to_complex_matrix(z,N,U);
[~,grad_complex] = compute_crlb_objective_grad(cfg,smat,Wmat);
[c,Gcomplex] = compute_constraints_grad(cfg,channels,smat,Wmat);
grad_f = complex_matrix_to_real_vector(grad_complex);
Gz = zeros(numel(z),numel(c));
for j = 1:numel(c)
    Gz(:,j) = complex_matrix_to_real_vector(Gcomplex(:,:,j));
end

active = c >= -opts.active_tol;
lambda = zeros(numel(c),1);
if any(active)
    lambda(active) = solve_nnls(Gz(:,active),-grad_f);
end

lag_grad = grad_f+Gz*lambda;
primal_residual = max([0;c]);
normalized_stationarity = norm(lag_grad)/max(1,norm(grad_f));
normalized_complementarity = ...
    max(abs(lambda.*c))/max(1,norm(lambda,inf));
normalized_kkt = max([primal_residual,normalized_stationarity, ...
    normalized_complementarity]);

diag.normalized_kkt_residual = normalized_kkt;
diag.kkt_residual = normalized_kkt;
diag.primal_residual = primal_residual;
diag.normalized_stationarity = normalized_stationarity;
diag.normalized_complementarity = normalized_complementarity;
diag.active_count = nnz(active);
diag.lambda = lambda;
diag.constraints = c;
end

function x = solve_nnls(A,b)
if isempty(A)
    x = zeros(0,1);
elseif exist('lsqnonneg','file') == 2
    x = lsqnonneg(A,b);
else
    x = zeros(size(A,2),1);
    alpha = 1/max(norm(A,2)^2,realmin);
    for k = 1:500
        xnew = max(0,x-alpha*(A'*(A*x-b)));
        if norm(xnew-x) <= 1e-8*max(1,norm(x))
            x = xnew;
            break;
        end
        x = xnew;
    end
end
end
