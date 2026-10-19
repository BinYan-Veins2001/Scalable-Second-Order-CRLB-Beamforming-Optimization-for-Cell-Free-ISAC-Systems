function metrics = evaluate_lifted_solution(cfg, channels, smat, Wcell)
%EVALUATE_LIFTED_SOLUTION Evaluate lifted W_u solution.

M = cfg.sys.M;
U = cfg.sys.U;

Hcell = cell(U,1);
for u = 1:U
    hu = channels.h_stack(:,u);
    Hcell{u} = hermitianize(hu*hu');
end

sinr = zeros(U,1);
for u = 1:U
    sig = real(trace(Hcell{u} * Wcell{u}));
    interf = 0;
    for i = 1:U
        if i ~= u
            interf = interf + real(trace(Hcell{u} * Wcell{i}));
        end
    end
    sinr(u) = sig / (interf + cfg.comm.noise_power);
end

ap_power = zeros(M,1);
for m = 1:M
    for u = 1:U
        ap_power(m) = ap_power(m) + real(trace(smat.Bm{m} * Wcell{u}));
    end
end

gout = compute_g_values_from_W(cfg, smat, Wcell);

rank_gap = zeros(U,1);
ranks_num = zeros(U,1);
for u = 1:U
    eigvals = sort(real(eig(hermitianize(Wcell{u}))), 'descend');
    eigvals(eigvals < 0 & eigvals > -1e-8) = 0;
    if isempty(eigvals)
        eigvals = 0;
    end
    rank_gap(u) = max(real(trace(Wcell{u})),0) - max(eigvals(1),0);
    ranks_num(u) = sum(eigvals > max(1e-8, 1e-5*max(eigvals(1),0)));
end

metrics.sinr = sinr;
metrics.ap_power = ap_power;
metrics.crlb = gout.crlb;
metrics.crlb_ratio = gout.ratio;
metrics.ga = gout.ga;
metrics.gb = gout.gb;
metrics.gab = gout.gab;
metrics.det = gout.det;
metrics.rank_gap = rank_gap;
metrics.rank_gap_sum = sum(max(rank_gap,0));
metrics.rank_est = ranks_num;
metrics.sinr_margin = sinr - cfg.comm.Gamma;
metrics.power_margin = cfg.power.Pmax - ap_power;
end
