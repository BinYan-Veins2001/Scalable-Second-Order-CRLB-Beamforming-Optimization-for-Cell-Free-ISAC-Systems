function [obj, grad, gout] = compute_crlb_objective_grad(cfg, smat, Wmat)
%COMPUTE_CRLB_OBJECTIVE_GRAD Original-vector CRLB-ratio/log-ratio objective and gradient.
%
% Supported objectives for first-order vector-space algorithms:
%   cfg.alm.objective_type = 'crlb_ratio'
%       obj = (g_a+g_b)/(g_a*g_b-g_ab^2)
%   cfg.alm.objective_type = 'log_crlb_ratio'
%       obj = log(g_a+g_b) - log(g_a*g_b-g_ab^2)
%
% The returned gout.crlb is always the physical CRLB with the constant
% prefactor in the TWC paper, irrespective of the internal objective.
%
% Gradient convention:
%   grad has the same size as Wmat. It is the real Euclidean gradient in
%   complex representation. A gradient step Wmat <- Wmat - alpha*grad is
%   valid for real-valued objectives.

U = cfg.sys.U;
beta2 = cfg.radio.effective_bandwidth_sq(:);

Aa = smat.Aa;
Ab = smat.Ab;
Aab = smat.Aab;

% Compute g-values directly from Wmat.
ga = 0; gb = 0; gab = 0;
for u = 1:U
    wu = Wmat(:,u);
    ga  = ga  + beta2(u) * real(wu' * Aa  * wu);
    gb  = gb  + beta2(u) * real(wu' * Ab  * wu);
    gab = gab + beta2(u) * real(wu' * Aab * wu);
end

ga = max(real(ga), 0);
gb = max(real(gb), 0);
gab = real(gab);
S = ga + gb;
D = ga*gb - gab^2;

if isfield(cfg,'alm') && isfield(cfg.alm,'det_floor')
    det_floor = cfg.alm.det_floor;
else
    det_floor = 1e-30;
end
D_eff = max(D, det_floor);
S_eff = max(S, realmin);

ratio = S_eff / D_eff;

if isfield(cfg,'alm') && isfield(cfg.alm,'objective_type')
    objective_type = lower(cfg.alm.objective_type);
else
    objective_type = 'crlb_ratio';
end

switch objective_type
    case 'crlb_ratio'
        obj = ratio;

        % Derivatives of F = S/D with respect to ga, gb, gab.
        dF_dga  = (D_eff - S_eff*gb) / (D_eff^2);
        dF_dgb  = (D_eff - S_eff*ga) / (D_eff^2);
        dF_dgab = 2*S_eff*gab / (D_eff^2);

    case 'log_crlb_ratio'
        % log F = log(S) - log(D). It has the same minimizers as F over D>0.
        obj = log(S_eff) - log(D_eff);

        dF_dga  = 1/S_eff - gb/D_eff;
        dF_dgb  = 1/S_eff - ga/D_eff;
        dF_dgab = 2*gab/D_eff;

    otherwise
        error('Unknown cfg.alm.objective_type: %s', objective_type);
end

Q = dF_dga * Aa + dF_dgb * Ab + dF_dgab * Aab;
Q = hermitianize(Q);

grad = zeros(size(Wmat));
for u = 1:U
    grad(:,u) = 2 * beta2(u) * Q * Wmat(:,u);
end

crlb_const = cfg.sensing.N0 * cfg.radio.c^2 / (8*pi^2*cfg.sys.Nr*cfg.radio.frame_duration);

gout.ga = ga;
gout.gb = gb;
gout.gab = gab;
gout.det = D;
gout.det_eff = D_eff;
gout.sum_ab = S;
gout.ratio = ratio;
gout.crlb = crlb_const * ratio;
gout.obj_for_alm = obj;
gout.objective_type = objective_type;
end
