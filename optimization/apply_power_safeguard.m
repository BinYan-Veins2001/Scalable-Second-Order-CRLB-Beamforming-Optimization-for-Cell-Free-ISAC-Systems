function Wmat = apply_power_safeguard(cfg, smat, Wmat)
%APPLY_POWER_SAFEGUARD Optional safeguard against severe power explosion.
% This is not the main constraint-handling mechanism. It is disabled unless
% cfg.alm.use_power_safeguard is true.

if ~isfield(cfg,'alm') || ~isfield(cfg.alm,'use_power_safeguard') || ~cfg.alm.use_power_safeguard
    return;
end

M = cfg.sys.M;
U = cfg.sys.U;
Pmax = cfg.power.Pmax;
if isfield(cfg.alm,'power_safeguard_factor')
    factor = cfg.alm.power_safeguard_factor;
else
    factor = 5;
end

for m = 1:M
    Bm = smat.Bm{m};
    p_m = 0;
    for u = 1:U
        p_m = p_m + real(Wmat(:,u)' * Bm * Wmat(:,u));
    end
    if p_m > factor * Pmax
        idx = block_indices(m,cfg.sys.Nt);
        scale = sqrt(factor * Pmax / max(p_m, realmin));
        for u = 1:U
            Wmat(idx,u) = scale * Wmat(idx,u);
        end
    end
end
end
