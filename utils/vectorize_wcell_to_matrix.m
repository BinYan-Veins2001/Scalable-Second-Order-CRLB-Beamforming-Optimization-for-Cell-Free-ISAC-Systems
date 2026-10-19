function Wmat = vectorize_wcell_to_matrix(cfg, wcell)
%VECTORIZE_WCELL_TO_MATRIX Convert wcell into N-by-U matrix of stacked users.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;

Wmat = zeros(N,U);
for u = 1:U
    for m = 1:M
        Wmat(block_indices(m,Nt),u) = wcell{m,u};
    end
end
end
