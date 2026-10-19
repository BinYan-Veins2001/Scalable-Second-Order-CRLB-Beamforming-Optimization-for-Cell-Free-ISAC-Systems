function wcell = recover_beamformers_from_W(cfg, Wcell)
%RECOVER_BEAMFORMERS_FROM_W Recover w_{m,u} from lifted W_u by principal eigenvector.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;
wcell = cell(M,U);

for u = 1:U
    Wu = hermitianize(Wcell{u});
    [V,D] = eig(Wu);
    [lambda_max, idx] = max(real(diag(D)));
    wu = sqrt(max(lambda_max,0)) * V(:,idx);
    for m = 1:M
        wcell{m,u} = wu(block_indices(m,Nt));
    end
end
end
