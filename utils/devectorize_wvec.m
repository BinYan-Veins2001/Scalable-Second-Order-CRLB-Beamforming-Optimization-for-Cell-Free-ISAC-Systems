function wcell = devectorize_wvec(cfg, wvec)
%DEVECTORIZE_WVEC Convert stacked vector into wcell{m,u} beamformers.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;

wcell = cell(M,U);
for u = 1:U
    idxu = block_indices(u,N);
    wu = wvec(idxu);
    for m = 1:M
        idxm = block_indices(m,Nt);
        wcell{m,u} = wu(idxm);
    end
end
end
