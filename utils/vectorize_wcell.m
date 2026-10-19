function wvec = vectorize_wcell(cfg, wcell)
%VECTORIZE_WCELL Convert wcell{m,u} beamformers into one stacked vector.
% Ordering: [w_1; w_2; ...; w_U], where w_u=[w_{1,u};...;w_{M,u}].

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;

wvec = zeros(N*U,1);
for u = 1:U
    idxu = block_indices(u,N);
    wu = zeros(N,1);
    for m = 1:M
        idxm = block_indices(m,Nt);
        wu(idxm) = wcell{m,u};
    end
    wvec(idxu) = wu;
end
end
