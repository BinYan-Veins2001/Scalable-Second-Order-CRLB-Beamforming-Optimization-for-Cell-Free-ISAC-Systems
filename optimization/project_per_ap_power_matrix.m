function Wmat = project_per_ap_power_matrix(cfg, Wmat)
%PROJECT_PER_AP_POWER_MATRIX Exact projection onto per-AP transmit-power balls.
%
% Projection set:
%   P = {Wmat : sum_u ||w_{m,u}||^2 <= Pmax, m=1,...,M}.
%
% Wmat is N x U, N = M*Nt. Rows (m-1)*Nt+1 : m*Nt belong to AP m.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
Pmax = cfg.power.Pmax;

for m = 1:M
    idx = (m-1)*Nt + (1:Nt);
    pm = 0;
    for u = 1:U
        pm = pm + norm(Wmat(idx,u))^2;
    end

    if pm > Pmax
        Wmat(idx,:) = sqrt(Pmax / max(pm, realmin)) * Wmat(idx,:);
    end
end

end
