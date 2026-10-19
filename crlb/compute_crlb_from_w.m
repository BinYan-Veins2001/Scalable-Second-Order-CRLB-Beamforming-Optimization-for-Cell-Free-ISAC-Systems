function [crlb, gout] = compute_crlb_from_w(cfg, smat, wcell)
%COMPUTE_CRLB_FROM_W Compute sensing CRLB from beamforming cell array.
% wcell{m,u}: Nt x 1 beamformer from AP m to user u.

U = cfg.sys.U;
M = cfg.sys.M;
Nt = cfg.sys.Nt;
Wcell = cell(U,1);
for u = 1:U
    wu = zeros(M*Nt,1);
    for m = 1:M
        wu(block_indices(m,Nt)) = wcell{m,u};
    end
    Wcell{u} = wu * wu';
end
gout = compute_g_values_from_W(cfg, smat, Wcell);
crlb = gout.crlb;
end
