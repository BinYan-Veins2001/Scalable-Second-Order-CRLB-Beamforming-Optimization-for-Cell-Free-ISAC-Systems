function channels = generate_channels(cfg, scenario)
%GENERATE_CHANNELS Generate cell-free downlink communication channels.
% h{m,u} is Nt x 1, from AP m to user u.
% h_stack(:,u) = [h{1,u}; ...; h{M,u}].

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
K = cfg.comm.rician_K;

h = cell(M,U);
h_stack = zeros(M*Nt, U);
large_scale = zeros(M,U);
dist = zeros(M,U);
theta_mu = zeros(M,U);

for m = 1:M
    ap = scenario.ap_pos(m,:);
    for u = 1:U
        usr = scenario.user_pos(u,:);
        dx = usr(1) - ap(1);
        dy = usr(2) - ap(2);
        d = max(sqrt(dx^2 + dy^2), cfg.comm.pathloss_ref_dist);
        theta = atan2(dy, dx);

        beta0 = db2pow(cfg.comm.pathloss_ref_dB);
        beta = beta0 * (cfg.comm.pathloss_ref_dist / d)^cfg.comm.pathloss_exp;

        a_los = steering_ula(Nt, theta, cfg) / sqrt(Nt);
        h_nlos = (randn(Nt,1) + 1j*randn(Nt,1))/sqrt(2);
        h_mu = sqrt(beta) * (sqrt(K/(K+1))*a_los + sqrt(1/(K+1))*h_nlos);

        h{m,u} = h_mu;
        h_stack(block_indices(m,Nt), u) = h_mu;
        large_scale(m,u) = beta;
        dist(m,u) = d;
        theta_mu(m,u) = theta;
    end
end

channels.h = h;
channels.h_stack = h_stack;
channels.large_scale = large_scale;
channels.dist = dist;
channels.theta = theta_mu;
end
