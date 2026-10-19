function wcell = matrix_to_wcell(cfg, Wmat)
%MATRIX_TO_WCELL Convert N-by-U stacked beamforming matrix into cell form.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
wcell = cell(M,U);
for u = 1:U
    for m = 1:M
        wcell{m,u} = Wmat(block_indices(m,Nt),u);
    end
end
end
