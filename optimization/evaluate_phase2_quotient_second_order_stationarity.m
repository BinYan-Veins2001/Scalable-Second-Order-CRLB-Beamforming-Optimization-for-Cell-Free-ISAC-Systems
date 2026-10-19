function diag = evaluate_phase2_quotient_second_order_stationarity( ...
    cfg,channels,smat,z,opts)
%EVALUATE_PHASE2_QUOTIENT_SECOND_ORDER_STATIONARITY Phase-orbit SOKKT audit.
%
% The objective and constraints are invariant to an independent phase
% rotation of every user beam.  Consequently, an isolated strict SOKKT point
% cannot exist in the unreduced beamformer coordinates: the Lagrangian
% Hessian has U exact zero directions i*w_u.  This diagnostic removes both
% active-constraint normals and those phase-symmetry directions, and measures
% curvature on the resulting horizontal tangent space.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(opts,numel(z));

N = cfg.sys.M*cfg.sys.Nt;
U = cfg.sys.U;
z = real(z(:));
W = real_vector_to_complex_matrix(z,N,U);
[~,grad_complex] = compute_crlb_objective_grad(cfg,smat,W);
[c,Gcomplex] = compute_constraints_grad(cfg,channels,smat,W);
grad_f = complex_matrix_to_real_vector(grad_complex);
G = zeros(numel(z),numel(c));
for idx = 1:numel(c)
    G(:,idx) = complex_matrix_to_real_vector(Gcomplex(:,:,idx));
end

active = false(size(c));
active(1:U) = c(1:U) >= -opts.sinr_active_tol;
active(U+1:end) = c(U+1:end) >= -opts.power_active_tol;
if any(active)
    scales=max(vecnorm(G(:,active)),realmin);
    As=G(:,active)./scales;
    Qactive = orth(As,opts.rank_tol);
else
    Qactive = zeros(numel(z),0);
end

lambda = zeros(numel(c),1);
if any(active)
    mu=lsqnonneg(As,-grad_f);
    lambda(active)=mu./scales.';
end
lag_grad = grad_f+G*lambda;
primal_residual = max([0;c]);
stationarity_residual = norm(lag_grad)/max(1,norm(grad_f));
complementarity_residual = ...
    max(abs(lambda.*c))/max(1,norm(lambda,inf));
kkt_residual = max([primal_residual,stationarity_residual, ...
    complementarity_residual]);

phase_directions = zeros(numel(z),U);
for user_idx = 1:U
    phase_matrix = zeros(N,U);
    phase_matrix(:,user_idx) = 1i*W(:,user_idx);
    phase_directions(:,user_idx) = ...
        complex_matrix_to_real_vector(phase_matrix);
end
phase_directions = phase_directions- ...
    Qactive*(Qactive'*phase_directions);
Qphase = orth(phase_directions,opts.rank_tol);
Qexcluded = orth([Qactive Qphase],opts.rank_tol);
Zhorizontal = null(Qexcluded',opts.rank_tol);
horizontal_dim = size(Zhorizontal,2);

if horizontal_dim == 0
    lambda_min_horizontal = inf;
    lanczos_info = struct('iter',0,'ritz_values',[],'T',[], ...
        'ritz_vector',[]);
else
    Hv_reduced = @(v) Zhorizontal'*lagrangian_hess_vec( ...
        cfg,channels,smat,z,Zhorizontal*v,lambda,N,U,opts);
    lanczos_opts.max_iter = min(opts.lanczos_max_iter,horizontal_dim);
    lanczos_opts.tol = opts.lanczos_tol;
    lanczos_opts.seed = opts.seed;
    [lambda_min_horizontal,lanczos_info] = ...
        estimate_min_curvature_lanczos( ...
        Hv_reduced,horizontal_dim,lanczos_opts);
end
negative_curvature = max(0,-lambda_min_horizontal);
spectral_scale = max([1;abs(lanczos_info.ritz_values(:))]);
normalized_negative_curvature = negative_curvature/spectral_scale;

diag.primal_residual = primal_residual;
diag.stationarity_residual = stationarity_residual;
diag.complementarity_residual = complementarity_residual;
diag.kkt_residual = kkt_residual;
diag.lambda_min_horizontal = lambda_min_horizontal;
diag.negative_curvature = negative_curvature;
diag.normalized_negative_curvature = normalized_negative_curvature;
diag.sosp_residual = max(kkt_residual,normalized_negative_curvature);
diag.constraints = c;
diag.active = active;
diag.active_count = nnz(active);
diag.active_rank = size(Qactive,2);
diag.phase_rank = size(Qphase,2);
diag.horizontal_dim = horizontal_dim;
diag.lambda = lambda;
diag.lanczos = lanczos_info;
diag.hessian_model = opts.hessian_model;
end

function Hv = lagrangian_hess_vec( ...
        cfg,channels,smat,z,v,lambda,N,U,opts)
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
h = opts.fd_rel_step*max(1,norm(z))/vnorm;
gp = lagrangian_gradient( ...
    cfg,channels,smat,z+h*v,lambda,N,U);
gm = lagrangian_gradient( ...
    cfg,channels,smat,z-h*v,lambda,N,U);
Hv = real((gp-gm)/(2*h));
end

function grad_L = lagrangian_gradient( ...
        cfg,channels,smat,z,lambda,N,U)
W = real_vector_to_complex_matrix(z,N,U);
[~,grad_complex] = compute_crlb_objective_grad(cfg,smat,W);
[~,Gcomplex] = compute_constraints_grad(cfg,channels,smat,W);
grad_f = complex_matrix_to_real_vector(grad_complex);
G = zeros(numel(z),numel(lambda));
for idx = 1:numel(lambda)
    G(:,idx) = complex_matrix_to_real_vector(Gcomplex(:,:,idx));
end
grad_L = grad_f+G*lambda;
end

function opts = fill_defaults(opts,dim)
if ~isfield(opts,'sinr_active_tol'); opts.sinr_active_tol=1e-3; end
if ~isfield(opts,'power_active_tol'); opts.power_active_tol=1e-5; end
if ~isfield(opts,'rank_tol'); opts.rank_tol=1e-10; end
if ~isfield(opts,'hessian_model'); opts.hessian_model='exact_lagrangian'; end
if ~isfield(opts,'fd_rel_step'); opts.fd_rel_step=1e-5; end
if ~isfield(opts,'lanczos_max_iter')
    opts.lanczos_max_iter=min(40,dim);
end
if ~isfield(opts,'lanczos_tol'); opts.lanczos_tol=1e-9; end
if ~isfield(opts,'seed'); opts.seed=2026; end
end
