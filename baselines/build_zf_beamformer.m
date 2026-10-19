function [wcell,info] = build_zf_beamformer(cfg,channels,smat,opts)
%BUILD_ZF_BEAMFORMER Paper-style local ZF with optimized power allocation.
%
% For Nt >= U, the direction at every AP follows equation (35):
%
% d_{m,u} = [H_m^H(H_mH_m^H)^(-1)]_{:,u} /
%           sqrt([(H_mH_m^H)^(-1)]_{u,u}),
%
% and w_{m,u}=sqrt(P_{m,u})d_{m,u}.  The square-root powers are optimized
% by OPTIMIZE_FIXED_DIRECTION_POWER_ALLOCATION rather than assigned equally.
%
% When Nt < U, the inverse required by (35) does not exist.  A local
% regularized-ZF direction is used as a dimension-safe extension; INFO
% records this explicitly so it is not confused with exact local ZF.

if nargin < 3
    smat = [];
end
if nargin < 4
    opts = struct();
end

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
directions = cell(M,U);
exact_local_zf = Nt >= U;

for m = 1:M
    Hm = zeros(U,Nt);
    for u = 1:U
        Hm(u,:) = channels.h{m,u}';
    end
    gram = Hm*Hm';
    gram_scale = max(real(trace(gram))/max(U,1),realmin);
    if exact_local_zf
        gram_inverse = pinv(hermitianize(gram),1e-12*gram_scale);
        raw_directions = Hm'*gram_inverse;
    else
        regularization = 1e-3*gram_scale;
        raw_directions = Hm'/(gram+regularization*eye(U));
    end

    for u = 1:U
        d = raw_directions(:,u);
        if norm(d) <= realmin
            d = channels.h{m,u};
        end
        directions{m,u} = d/norm(d);
    end
end

opts.method_name = 'ZF-BF';
[wcell,info] = optimize_fixed_direction_power_allocation( ...
    cfg,channels,smat,directions,opts);
info.direction_type = ternary(exact_local_zf, ...
    'equation_35_local_zf','local_regularized_zf_for_Nt_less_than_U');
info.exact_equation_35 = exact_local_zf;
end

function value = ternary(condition,true_value,false_value)
if condition
    value = true_value;
else
    value = false_value;
end
end
