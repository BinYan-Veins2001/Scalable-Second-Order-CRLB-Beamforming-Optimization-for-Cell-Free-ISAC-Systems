function [c, G, raw] = compute_constraints_grad(cfg, channels, smat, Wmat)
%COMPUTE_CONSTRAINTS_GRAD Normalized inequality constraints and gradients.
%
% Constraints are c_j(Wmat) <= 0.
% First U constraints are SINR constraints:
%   Gamma*(interference + sigma2) - signal <= 0.
% They are normalized by sigma2.
% Last M constraints are per-AP power constraints:
%   p_m - Pmax <= 0.
% They are normalized by Pmax.
%
% G(:,:,j) is the Euclidean gradient of c_j with respect to Wmat.

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
N = M*Nt;
Gamma = cfg.comm.Gamma;
sigma2 = cfg.comm.noise_power;
Pmax = cfg.power.Pmax;

J = U + M;
c = zeros(J,1);
G = zeros(N,U,J);

raw.sinr_numer = zeros(U,1);
raw.sinr_denom = zeros(U,1);
raw.sinr = zeros(U,1);
raw.sinr_c = zeros(U,1);

% SINR constraints.
for u = 1:U
    hu = channels.h_stack(:,u);
    channel_outputs = hu' * Wmat;
    received_powers = abs(channel_outputs).^2;
    sig = received_powers(u);
    interf = sum(received_powers)-sig;
    raw_c = Gamma * (interf + sigma2) - sig;
    c(u) = raw_c / max(sigma2, realmin);

    for i = 1:U
        if i == u
            G(:,i,u) = (-2 * hu * channel_outputs(i)) / ...
                max(sigma2, realmin);
        else
            G(:,i,u) = (2 * Gamma * hu * channel_outputs(i)) / ...
                max(sigma2, realmin);
        end
    end

    raw.sinr_numer(u) = max(sig,0);
    raw.sinr_denom(u) = interf + sigma2;
    raw.sinr(u) = max(sig,0) / max(interf + sigma2, realmin);
    raw.sinr_c(u) = raw_c;
end

% Per-AP power constraints.
raw.ap_power = zeros(M,1);
raw.power_c = zeros(M,1);
for m = 1:M
    rows = (m-1)*Nt+(1:Nt);
    Wm = Wmat(rows,:);
    p_m = norm(Wm,'fro')^2;
    raw_c = p_m - Pmax;
    j = U + m;
    c(j) = raw_c / max(Pmax, realmin);
    G(rows,:,j) = (2 * Wm) / max(Pmax, realmin);
    raw.ap_power(m) = max(p_m,0);
    raw.power_c(m) = raw_c;
end

raw.max_violation = max(max(c,0));
raw.min_sinr = min(raw.sinr);
raw.max_ap_power = max(raw.ap_power);
end
