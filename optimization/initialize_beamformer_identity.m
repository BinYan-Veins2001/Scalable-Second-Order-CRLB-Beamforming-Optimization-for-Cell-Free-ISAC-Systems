function wcell = initialize_beamformer_identity(cfg)
%INITIALIZE_BEAMFORMER_IDENTITY Isotropic vector-space counterpart of W=alpha*I.
%
% SDR-SCA initializes lifted variables as:
%   W_u^(0) = alpha I_(M Nt).
%
% A vector method cannot set w_u w_u^H = alpha I exactly because
% w_u w_u^H is rank-one. We therefore draw
%   w_u ~ CN(0, alpha I)
% so that E[w_u w_u^H] = alpha I.
%
% The final per-AP normalization keeps the random initialization power-bounded.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M * Nt;
Pmax = cfg.power.Pmax;

if isfield(cfg, 'alm') && isfield(cfg.alm, 'identity_power_scale')
    scale = cfg.alm.identity_power_scale;
else
    scale = cfg.opt.init_power_scale;
end

alpha = scale * M * Pmax / U / N;

wcell = cell(M,U);
for u = 1:U
    wu = sqrt(alpha) * (randn(N,1) + 1j*randn(N,1)) / sqrt(2);
    for m = 1:M
        idx = block_indices(m,Nt);
        wcell{m,u} = wu(idx);
    end
end

% Match the SCA identity initialization's per-AP power budget behavior:
% E[sum_u ||w_{m,u}||^2] = scale * Pmax. Random draws may exceed Pmax,
% so scale each AP block if necessary.
for m = 1:M
    p_m = 0;
    for u = 1:U
        p_m = p_m + norm(wcell{m,u})^2;
    end
    if p_m > Pmax
        fac = sqrt(Pmax / max(p_m, realmin));
        for u = 1:U
            wcell{m,u} = fac * wcell{m,u};
        end
    end
end
end
