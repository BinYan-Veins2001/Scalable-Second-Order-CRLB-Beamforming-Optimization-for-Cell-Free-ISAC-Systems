function [Hv,info] = exact_lagrangian_hess_vec( ...
    cfg,channels,smat,Wmat,Vmat,lambda)
%EXACT_LAGRANGIAN_HESS_VEC Exact matrix-free Lagrangian Hessian action.
%
% The returned complex matrix represents the real Euclidean Hessian-vector
% product for
%   L(W,lambda) = f(W) + sum_j lambda(j)*c_j(W),
% where f is the CRLB-ratio or log-CRLB-ratio objective and c contains the
% normalized SINR and per-AP power inequalities. No ambient Hessian matrix
% is assembled and no finite differences are used.

N = cfg.sys.M*cfg.sys.Nt;
U = cfg.sys.U;
M = cfg.sys.M;
if ~isequal(size(Wmat),[N,U]) || ~isequal(size(Vmat),[N,U])
    error('exact_lagrangian_hess_vec:DimensionMismatch', ...
        'Wmat and Vmat must both be %d-by-%d.',N,U);
end
lambda = real(lambda(:));
if numel(lambda) ~= U+M
    error('exact_lagrangian_hess_vec:MultiplierDimensionMismatch', ...
        'Expected %d Lagrange multipliers, got %d.',U+M,numel(lambda));
end

[HfV,obj_info] = objective_hess_vec(cfg,smat,Wmat,Vmat);
HcV = constraint_lagrangian_hess_vec(cfg,channels,Vmat,lambda);
Hv = HfV+HcV;

info.objective = obj_info;
info.sinr_multiplier_count = nnz(lambda(1:U));
info.power_multiplier_count = nnz(lambda(U+1:end));
info.used_finite_difference = false;
end

function [Hv,info] = objective_hess_vec(cfg,smat,Wmat,Vmat)
beta2 = cfg.radio.effective_bandwidth_sq(:).';
matrices = {smat.Aa,smat.Ab,smat.Aab};

AW = cell(3,1);
AV = cell(3,1);
g = zeros(3,1);
dg = zeros(3,1);
for r = 1:3
    AW{r} = matrices{r}*Wmat;
    AV{r} = matrices{r}*Vmat;
    quadratic = real(sum(conj(Wmat).*AW{r},1));
    directional = 2*real(sum(conj(Vmat).*AW{r},1));
    g(r) = sum(beta2.*quadratic);
    dg(r) = sum(beta2.*directional);
end

[alpha,Bphi,outer_info] = objective_outer_derivatives(cfg,g);
dalpha = Bphi*dg;
Hv = zeros(size(Wmat));
for r = 1:3
    Hv = Hv+alpha(r)*AV{r}+dalpha(r)*AW{r};
end
Hv = 2*Hv.*beta2;

info.g = g;
info.dg = dg;
info.alpha = alpha;
info.outer_hessian = Bphi;
info.outer = outer_info;
end

function [alpha,Bphi,info] = objective_outer_derivatives(cfg,g)
a = real(g(1));
b = real(g(2));
c = real(g(3));
s = a+b;
d = a*b-c^2;
if isfield(cfg,'alm') && isfield(cfg.alm,'det_floor')
    det_floor = cfg.alm.det_floor;
else
    det_floor = 1e-30;
end

% The paper objective is smooth only inside the positive-definite FIM
% domain. The numerical max safeguards in the value routine are not part of
% that smooth objective and must be inactive for an exact classical HVP.
if ~(isfinite(s) && isfinite(d) && s > realmin && d > det_floor)
    error('exact_lagrangian_hess_vec:NonSmoothObjectiveDomain', ...
        ['Exact HVP requires g_a+g_b>realmin and ' ...
         'g_a*g_b-g_ab^2>det_floor. Got S=%.4e and D=%.4e.'],s,d);
end

alpha_log = [1/s-b/d; 1/s-a/d; 2*c/d];
B_log = [ ...
    -1/s^2+b^2/d^2, -1/s^2+c^2/d^2, -2*b*c/d^2; ...
    -1/s^2+c^2/d^2, -1/s^2+a^2/d^2, -2*a*c/d^2; ...
    -2*b*c/d^2,      -2*a*c/d^2,      2/d+4*c^2/d^2];
B_log = real((B_log+B_log')/2);

if isfield(cfg,'alm') && isfield(cfg.alm,'objective_type')
    objective_type = lower(cfg.alm.objective_type);
else
    objective_type = 'crlb_ratio';
end
switch objective_type
    case 'log_crlb_ratio'
        alpha = alpha_log;
        Bphi = B_log;
    case 'crlb_ratio'
        ratio = s/d;
        alpha = ratio*alpha_log;
        Bphi = ratio*(B_log+alpha_log*alpha_log');
    otherwise
        error('Unknown cfg.alm.objective_type: %s',objective_type);
end
Bphi = real((Bphi+Bphi')/2);

info.objective_type = objective_type;
info.sum_ab = s;
info.det = d;
info.det_floor = det_floor;
info.smooth_domain = true;
end

function Hv = constraint_lagrangian_hess_vec(cfg,channels,Vmat,lambda)
M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
sigma2 = max(cfg.comm.noise_power,realmin);
Pmax = max(cfg.power.Pmax,realmin);

% Every H_u=h_u*h_u' is rank one. The matrix expression below evaluates
% sum_u lambda_u H_u V with stream-dependent SINR signs in O(N U^2).
H = channels.h_stack;
stream_weights = cfg.comm.Gamma*ones(U,U);
stream_weights(1:U+1:end) = -1;
weighted_outputs = (lambda(1:U).*ones(1,U)).*stream_weights.*(H'*Vmat);
Hv = (2/sigma2)*(H*weighted_outputs);

for m = 1:M
    multiplier = lambda(U+m);
    if multiplier == 0
        continue;
    end
    rows = (m-1)*Nt+(1:Nt);
    Hv(rows,:) = Hv(rows,:)+(2*multiplier/Pmax)*Vmat(rows,:);
end
end
