function [s, stats] = trust_region_tcg(g, Hvfun, Delta, opts)
%TRUST_REGION_TCG Safeguarded Steihaug-Toint truncated CG solver.
%
% Approximately solves
%
%   min_s  g' s + 0.5 s' H s   subject to ||s|| <= Delta
%
% using Hessian-vector products only.
%
% Compared with the first PTR-ALM package version, this implementation uses
% a stricter negative-curvature test. A direction is marked as true negative
% curvature only when
%
%   p' H p / ||p||^2 <= -opts.nc_tol.
%
% Very small positive/near-zero curvature is treated as a boundary/flat
% curvature event, but is not counted as negative curvature. This avoids
% false negative-curvature detections caused by finite-difference HVP noise.

if nargin < 4
    opts = struct();
end

if ~isfield(opts,'max_iter'); opts.max_iter = min(100, numel(g)); end
if ~isfield(opts,'tol'); opts.tol = 1e-3; end
if ~isfield(opts,'nc_tol'); opts.nc_tol = 1e-8; end
if ~isfield(opts,'flat_tol'); opts.flat_tol = 1e-12; end
if ~isfield(opts,'detect_negative_curvature')
    opts.detect_negative_curvature = true;
end
if ~isfield(opts,'forcing_power'); opts.forcing_power = []; end
if ~isfield(opts,'forcing_constant'); opts.forcing_constant = 1; end
if ~isfield(opts,'absolute_floor'); opts.absolute_floor = 0; end
if ~isfield(opts,'verbose'); opts.verbose = false; end

Delta = max(Delta, 0);
g = real(g(:));
dim = numel(g);

s = zeros(dim,1);
r = g;
p = -r;
rTr = r' * r;
r0 = sqrt(max(rTr,0));

% The released solver retains its original fixed relative/absolute test.
% A positive forcing_power activates an inexact-Newton forcing sequence,
% ||r_CG|| <= forcing_constant*||g||^forcing_power.  In particular,
% forcing_power=2 is the inner-accuracy condition used by the controlled
% local quadratic-rate diagnostic.
if isempty(opts.forcing_power) || ~isfinite(opts.forcing_power)
    residual_target = opts.tol * max(1,r0);
else
    residual_target = opts.forcing_constant * r0^opts.forcing_power;
    residual_target = max(residual_target,opts.absolute_floor);
end

stats.iter = 0;
stats.reason = 'max_iter';
stats.negative_curvature = false;
stats.near_zero_curvature = false;
stats.hit_boundary = false;
stats.min_curvature_est = inf;
stats.final_res_norm = r0;
stats.residual_target = residual_target;
stats.model_value = 0;
stats.model_decrease = 0;
stats.step_norm = 0;

if Delta <= 0 || r0 == 0
    stats.reason = 'zero_delta_or_gradient';
    return;
end

for k = 1:opts.max_iter
    Hp = real(Hvfun(p));
    Hp = Hp(:);

    pHp = p' * Hp;
    pp = p' * p;
    curv = pHp / max(pp, realmin);
    stats.min_curvature_est = min(stats.min_curvature_est, curv);

    % True negative curvature. The sign is important: do not treat small
    % positive curvature as negative curvature.
    if opts.detect_negative_curvature && curv <= -opts.nc_tol
        tau = tau_to_boundary(s, p, Delta);
        s = s + tau * p;
        stats.iter = k;
        stats.reason = 'negative_curvature';
        stats.negative_curvature = true;
        stats.hit_boundary = true;
        break;
    end

    % NC-off must not silently take the same boundary step through the
    % flat-curvature branch. Keep the accumulated positive-curvature step;
    % if no step exists yet, use a conservative gradient Cauchy step.
    if ~opts.detect_negative_curvature && curv <= -opts.nc_tol
        if norm(s)<=realmin
            length=min(0.5*Delta,r0/max(1,norm(Hp)/sqrt(max(pp,realmin))));
            s=-length*g/r0;
        end
        stats.iter=k;
        stats.reason='negative_curvature_disabled';
        break;
    end

    % Near-zero curvature: go to the boundary along the current descent
    % direction, but do not mark it as negative curvature.
    if curv <= opts.flat_tol
        tau = tau_to_boundary(s, p, Delta);
        s = s + tau * p;
        stats.iter = k;
        stats.reason = 'near_zero_curvature_boundary';
        stats.near_zero_curvature = true;
        stats.hit_boundary = true;
        break;
    end

    alpha = rTr / max(pHp, realmin);
    s_next = s + alpha * p;

    if norm(s_next) >= Delta
        tau = tau_to_boundary(s, p, Delta);
        s = s + tau * p;
        stats.iter = k;
        stats.reason = 'boundary';
        stats.hit_boundary = true;
        break;
    end

    s = s_next;
    r_next = r + alpha * Hp;
    r_next_Tr_next = r_next' * r_next;
    res_norm = sqrt(max(r_next_Tr_next,0));

    if res_norm <= residual_target
        stats.iter = k;
        stats.reason = 'residual_tol';
        stats.final_res_norm = res_norm;
        break;
    end

    beta = r_next_Tr_next / max(rTr, realmin);
    p = -r_next + beta * p;
    r = r_next;
    rTr = r_next_Tr_next;
    stats.iter = k;
    stats.final_res_norm = res_norm;
end

Hs = real(Hvfun(s));
stats.model_value = g' * s + 0.5 * s' * Hs;
stats.model_decrease = -stats.model_value;
stats.step_norm = norm(s);

end

function tau = tau_to_boundary(s, p, Delta)
% Solve ||s + tau p|| = Delta for tau >= 0.
a = p' * p;
b = 2 * (s' * p);
c = s' * s - Delta^2;
disc = max(b^2 - 4*a*c, 0);
tau = (-b + sqrt(disc)) / max(2*a, realmin);
tau = max(tau, 0);
end
