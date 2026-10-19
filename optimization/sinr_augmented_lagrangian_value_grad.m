function [Lval, gradL, info] = sinr_augmented_lagrangian_value_grad(cfg, channels, smat, Wmat, mu, rho, obj_scale)
%SINR_AUGMENTED_LAGRANGIAN_VALUE_GRAD ALM value/gradient using only SINR constraints.
%
% This routine is intended for projected-power methods:
%   - SINR constraints are handled by shifted inequality ALM;
%   - per-AP power constraints are handled separately by projection.
%
% It reuses the project's compute_crlb_objective_grad and
% compute_constraints_grad. The first U entries of c are assumed to be SINR
% constraints; entries U+1:U+M are power constraints and are used only for
% diagnostics.

if nargin < 7 || isempty(obj_scale)
    obj_scale = 1;
end
obj_scale = max(obj_scale, realmin);

U = cfg.sys.U;
M = cfg.sys.M;

[F, gradF, gout] = compute_crlb_objective_grad(cfg, smat, Wmat);
[c_all, G_all, raw] = compute_constraints_grad(cfg, channels, smat, Wmat);

c_sinr = c_all(1:U);
mu = mu(:);
if numel(mu) ~= U
    error('sinr_augmented_lagrangian_value_grad expects numel(mu)=U.');
end

shifted = c_sinr + mu / rho;
active = shifted > 0;

Lval = F / obj_scale;
gradL = gradF / obj_scale;
penalty_val = 0;

for u = 1:U
    if active(u)
        coeff = rho * shifted(u);
        penalty_val = penalty_val + 0.5 * rho * shifted(u)^2;
        gradL = gradL + coeff * G_all(:,:,u);
    end
end

Lval = Lval + penalty_val;

if numel(c_all) >= U + M
    c_power = c_all(U+1:U+M);
    max_power_violation = max(max(c_power,0));
else
    max_power_violation = max(0, raw.max_ap_power / max(cfg.power.Pmax, realmin) - 1);
end

info.F = F;
info.optim_obj = F;
info.ratio = gout.ratio;
info.crlb = gout.crlb;
info.objective_type = gout.objective_type;
info.gout = gout;
info.c_sinr = c_sinr;
info.c_all = c_all;
info.raw = raw;
info.penalty_val = penalty_val;
info.active = active;
info.max_sinr_violation = max(max(c_sinr,0));
info.max_power_violation = max_power_violation;
info.max_violation = max(info.max_sinr_violation, info.max_power_violation);
info.grad_norm = norm(gradL(:));
info.obj_scale = obj_scale;

end
