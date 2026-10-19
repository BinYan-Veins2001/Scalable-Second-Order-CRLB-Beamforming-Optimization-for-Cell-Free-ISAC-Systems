function metrics = evaluate_beamforming_solution(cfg, channels, smat, wcell)
%EVALUATE_BEAMFORMING_SOLUTION Evaluate recovered/original beamformers.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;

sinr = zeros(U,1);
for u = 1:U
    tmp = 0;
    for m = 1:M
        tmp = tmp + channels.h{m,u}' * wcell{m,u};
    end
    sig = abs(tmp)^2;

    interf = 0;
    for i = 1:U
        if i == u
            continue;
        end
        tmpi = 0;
        for m = 1:M
            tmpi = tmpi + channels.h{m,u}' * wcell{m,i};
        end
        interf = interf + abs(tmpi)^2;
    end
    sinr(u) = sig / (interf + cfg.comm.noise_power);
end

ap_power = zeros(M,1);
for m = 1:M
    for u = 1:U
        ap_power(m) = ap_power(m) + norm(wcell{m,u})^2;
    end
end

[crlb, gout] = compute_crlb_from_w(cfg, smat, wcell);

metrics.sinr = sinr;
metrics.ap_power = ap_power;
metrics.crlb = crlb;
metrics.crlb_ratio = gout.ratio;
metrics.ga = gout.ga;
metrics.gb = gout.gb;
metrics.gab = gout.gab;
metrics.det = gout.det;
metrics.sinr_margin = sinr - cfg.comm.Gamma;
metrics.power_margin = cfg.power.Pmax - ap_power;

w_stack = zeros(M*Nt, U);
for u = 1:U
    for m = 1:M
        w_stack(block_indices(m,Nt),u) = wcell{m,u};
    end
end
metrics.w_stack = w_stack;
end
