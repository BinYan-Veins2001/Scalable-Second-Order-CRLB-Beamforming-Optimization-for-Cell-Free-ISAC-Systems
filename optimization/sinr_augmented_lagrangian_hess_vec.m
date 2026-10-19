function Hv = sinr_augmented_lagrangian_hess_vec(cfg, channels, smat, Wmat, Vmat, mu, rho, obj_scale)
%SINR_AUGMENTED_LAGRANGIAN_HESS_VEC Structured HVP for SINR-only ALM.
%
% This is the Hessian-vector product of sinr_augmented_lagrangian_value_grad
% with respect to the complex matrix variable Wmat, using the same complex
% Euclidean gradient convention as the rest of the package.

if nargin < 8 || isempty(obj_scale)
    obj_scale = 1;
end
obj_scale = max(obj_scale, realmin);

U = cfg.sys.U;
N = cfg.sys.M * cfg.sys.Nt;
sigma2 = cfg.comm.noise_power;
Gamma = cfg.comm.Gamma;

Hv = crlb_objective_hess_vec(cfg, smat, Wmat, Vmat) / obj_scale;

[c_all, G_all, ~] = compute_constraints_grad(cfg, channels, smat, Wmat);
c_sinr = c_all(1:U);
shifted = c_sinr + mu(:) / rho;

Hcell = cell(U,1);
for u = 1:U
    hu = channels.h_stack(:,u);
    Hcell{u} = hermitianize(hu * hu');
end

for u = 1:U
    if shifted(u) <= 0
        continue;
    end

    Hu = Hcell{u};
    dc = sinr_constraint_directional_derivative(Hu, Wmat, Vmat, u, U, Gamma, sigma2);
    HcV = zeros(N,U);
    for i = 1:U
        if i == u
            HcV(:,i) = (-2 * Hu * Vmat(:,i)) / max(sigma2, realmin);
        else
            HcV(:,i) = (2 * Gamma * Hu * Vmat(:,i)) / max(sigma2, realmin);
        end
    end

    coeff = rho * shifted(u);
    Hv = Hv + rho * dc * G_all(:,:,u) + coeff * HcV;
end

end

function Hv = crlb_objective_hess_vec(cfg, smat, Wmat, Vmat)
U = cfg.sys.U;
beta2 = cfg.radio.effective_bandwidth_sq(:);

Aa = smat.Aa;
Ab = smat.Ab;
Aab = smat.Aab;

[g, coeff] = crlb_g_and_coeff(cfg, smat, Wmat);

dg = zeros(3,1);
for u = 1:U
    wu = Wmat(:,u);
    vu = Vmat(:,u);
    dg(1) = dg(1) + beta2(u) * 2 * real(vu' * Aa * wu);
    dg(2) = dg(2) + beta2(u) * 2 * real(vu' * Ab * wu);
    dg(3) = dg(3) + beta2(u) * 2 * real(vu' * Aab * wu);
end

dc = coeff_directional_fd(cfg, g, dg);

Q = hermitianize(coeff(1) * Aa + coeff(2) * Ab + coeff(3) * Aab);
dQ = hermitianize(dc(1) * Aa + dc(2) * Ab + dc(3) * Aab);

Hv = zeros(size(Wmat));
for u = 1:U
    Hv(:,u) = 2 * beta2(u) * (Q * Vmat(:,u) + dQ * Wmat(:,u));
end
end

function dc = coeff_directional_fd(cfg, g, dg)
dgnorm = norm(dg);
if dgnorm <= 0
    dc = zeros(3,1);
    return;
end

eps_g = 1e-6 * (1 + norm(g)) / max(dgnorm, 1e-12);
cp = coeff_from_g(cfg, g + eps_g * dg);
cm = coeff_from_g(cfg, g - eps_g * dg);
dc = (cp - cm) / (2 * eps_g);
end

function [g, coeff] = crlb_g_and_coeff(cfg, smat, Wmat)
U = cfg.sys.U;
beta2 = cfg.radio.effective_bandwidth_sq(:);

Aa = smat.Aa;
Ab = smat.Ab;
Aab = smat.Aab;

ga = 0;
gb = 0;
gab = 0;
for u = 1:U
    wu = Wmat(:,u);
    ga  = ga  + beta2(u) * real(wu' * Aa  * wu);
    gb  = gb  + beta2(u) * real(wu' * Ab  * wu);
    gab = gab + beta2(u) * real(wu' * Aab * wu);
end

g = [max(real(ga), 0); max(real(gb), 0); real(gab)];
coeff = coeff_from_g(cfg, g);
end

function coeff = coeff_from_g(cfg, g)
ga = max(real(g(1)), 0);
gb = max(real(g(2)), 0);
gab = real(g(3));

S = max(ga + gb, realmin);
D = ga * gb - gab^2;

if isfield(cfg,'alm') && isfield(cfg.alm,'det_floor')
    det_floor = cfg.alm.det_floor;
else
    det_floor = 1e-30;
end
D = max(D, det_floor);

if isfield(cfg,'alm') && isfield(cfg.alm,'objective_type')
    objective_type = lower(cfg.alm.objective_type);
else
    objective_type = 'crlb_ratio';
end

switch objective_type
    case 'crlb_ratio'
        coeff = [
            (D - S * gb) / (D^2);
            (D - S * ga) / (D^2);
            2 * S * gab / (D^2)
        ];
    case 'log_crlb_ratio'
        coeff = [
            1 / S - gb / D;
            1 / S - ga / D;
            2 * gab / D
        ];
    otherwise
        error('Unknown cfg.alm.objective_type: %s', objective_type);
end
end

function dc = sinr_constraint_directional_derivative(Hu, Wmat, Vmat, u, U, Gamma, sigma2)
dsig = 2 * real(Vmat(:,u)' * Hu * Wmat(:,u));
dinterf = 0;
for i = 1:U
    if i ~= u
        dinterf = dinterf + 2 * real(Vmat(:,i)' * Hu * Wmat(:,i));
    end
end
dc = (Gamma * dinterf - dsig) / max(sigma2, realmin);
end
