function [Lval, gradL, info] = augmented_lagrangian_value_grad(cfg, channels, smat, Wmat, mu, rho, obj_scale)
%AUGMENTED_LAGRANGIAN_VALUE_GRAD Inequality AL value and gradient.
%
% Uses the shifted Powell-Hestenes-Rockafellar inequality AL term:
%   F(w)/obj_scale + sum_j rho/2 * max(0, c_j(w)+mu_j/rho)^2.
% The omitted constant -mu_j^2/(2rho) does not affect inner minimization.

if nargin < 7 || isempty(obj_scale)
    obj_scale = 1;
end
obj_scale = max(obj_scale, realmin);

[F, gradF, gout] = compute_crlb_objective_grad(cfg, smat, Wmat);
[c, G, raw] = compute_constraints_grad(cfg, channels, smat, Wmat);

shifted = c + mu(:)/rho;
active = shifted > 0;

Lval = F / obj_scale;
gradL = gradF / obj_scale;
penalty_val = 0;

for j = 1:length(c)
    if active(j)
        coeff = rho * shifted(j);
        penalty_val = penalty_val + 0.5 * rho * shifted(j)^2;
        gradL = gradL + coeff * G(:,:,j);
    end
end

Lval = Lval + penalty_val;

info.F = F;                 % internal objective used by ALM: ratio or log-ratio
info.optim_obj = F;
info.ratio = gout.ratio;     % true CRLB ratio, always non-log
info.crlb = gout.crlb;       % physical CRLB with paper constant prefactor
info.objective_type = gout.objective_type;
info.gout = gout;
info.c = c;
info.raw = raw;
info.penalty_val = penalty_val;
info.active = active;
info.max_violation = max(max(c,0));
info.grad_norm = norm(gradL(:));
info.obj_scale = obj_scale;
end
