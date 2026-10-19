function smat = build_sensing_matrices(cfg, scenario)
%BUILD_SENSING_MATRICES Build Aa, Ab, Aab, Bm matrices and geometry coefficients.
% The implementation follows the CRLB quadratic functions:
% ga = sum_u beta_eff_sq(u) * tr(Aa W_u), and similarly for gb/gab.

M = cfg.sys.M;
Nt = cfg.sys.Nt;
N = M*Nt;

ap = scenario.ap_pos;
tgt = scenario.target_pos(:).';

% Target-related AP distances and angles.
d = zeros(M,1);
theta = zeros(M,1);
at = cell(M,1);
for m = 1:M
    dx = tgt(1) - ap(m,1);
    dy = tgt(2) - ap(m,2);
    d(m) = max(sqrt(dx^2 + dy^2), cfg.sensing.d_ref);
    theta(m) = atan2(dy, dx);
    at{m} = steering_ula(Nt, theta(m), cfg);
end

% Geometry and sensing path loss for path transmit l -> target -> receive m.
a_ml = zeros(M,M);
b_ml = zeros(M,M);
beta_ml = zeros(M,M);
xi_ml = zeros(M,M);

for m = 1:M
    for l = 1:M
        dxm = tgt(1) - ap(m,1);
        dym = tgt(2) - ap(m,2);
        dxl = tgt(1) - ap(l,1);
        dyl = tgt(2) - ap(l,2);
        dm2 = max(d(m)^2, eps);
        dl2 = max(d(l)^2, eps);

        a_ml(m,l) = dxm/sqrt(dm2) + dxl/sqrt(dl2);
        b_ml(m,l) = dym/sqrt(dm2) + dyl/sqrt(dl2);

        beta_ml(m,l) = (cfg.sensing.kappa_ref_amp^2) * (cfg.sensing.d_ref^4) / (dm2 * dl2);
        if cfg.sensing.random_rcs_phase
            xi_ml(m,l) = cfg.sensing.rcs_amp * exp(1j*2*pi*rand);
        else
            xi_ml(m,l) = cfg.sensing.rcs_amp;
        end
    end
end

% Coefficients by transmitter l after summing over receiving m.
ca = zeros(M,1);
cb = zeros(M,1);
cab = zeros(M,1);
for l = 1:M
    for m = 1:M
        weight = beta_ml(m,l) * abs(xi_ml(m,l))^2;
        ca(l) = ca(l) + a_ml(m,l)^2 * weight;
        cb(l) = cb(l) + b_ml(m,l)^2 * weight;
        cab(l) = cab(l) + a_ml(m,l)*b_ml(m,l) * weight;
    end
end

Aa = zeros(N,N);
Ab = zeros(N,N);
Aab = zeros(N,N);
for l = 1:M
    idx = block_indices(l,Nt);
    A_l = conj(at{l}) * transpose(at{l});  % gives |a_t^T w_l|^2 = w_l^H A_l w_l
    Aa(idx,idx) = ca(l) * A_l;
    Ab(idx,idx) = cb(l) * A_l;
    Aab(idx,idx) = cab(l) * A_l;
end
Aa = hermitianize(Aa);
Ab = hermitianize(Ab);
Aab = hermitianize(Aab);

Bm = cell(M,1);
for m = 1:M
    B = zeros(N,N);
    idx = block_indices(m,Nt);
    B(idx,idx) = eye(Nt);
    Bm{m} = B;
end

smat.Aa = Aa;
smat.Ab = Ab;
smat.Aab = Aab;
smat.Bm = Bm;
smat.a_ml = a_ml;
smat.b_ml = b_ml;
smat.beta_ml = beta_ml;
smat.xi_ml = xi_ml;
smat.theta_target = theta;
smat.dist_target = d;
smat.at = at;
smat.ca = ca;
smat.cb = cb;
smat.cab = cab;
end
