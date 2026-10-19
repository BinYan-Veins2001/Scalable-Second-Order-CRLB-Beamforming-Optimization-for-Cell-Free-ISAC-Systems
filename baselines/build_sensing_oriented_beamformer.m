function [wcell,info] = build_sensing_oriented_beamformer( ...
    cfg,scenario,channels,smat,opts)
%BUILD_SENSING_ORIENTED_BEAMFORMER Eq. (36) with optimized power allocation.
%
% The target-pointing direction at AP m is a_t,m(theta_m)/sqrt(Nt), under
% the steering convention used by this package.  Hence
%
%       w_{m,u}=sqrt(P_{m,u}/Nt) a_t,m(theta_m),
%
% exactly as in equation (36).  Unlike the former equal-power sanity
% baseline, P_{m,u} is optimized subject to SINR and per-AP power limits.

if nargin < 3
    channels = [];
end
if nargin < 4
    smat = [];
end
if nargin < 5
    opts = struct();
end

M = cfg.sys.M;
U = cfg.sys.U;
Nt = cfg.sys.Nt;
directions = cell(M,U);

for m = 1:M
    dx = scenario.target_pos(1)-scenario.ap_pos(m,1);
    dy = scenario.target_pos(2)-scenario.ap_pos(m,2);
    theta = atan2(dy,dx);
    % The sensing model uses a_t^T w, so CONJ(a_t)/sqrt(Nt) is the
    % target-pointing transmit direction for STEERING_ULA's convention.
    d = conj(steering_ula(Nt,theta,cfg))/sqrt(Nt);
    for u = 1:U
        directions{m,u} = d;
    end
end

if isempty(channels)
    % Backward-compatible construction for geometry-only visualization.
    wcell = cell(M,U);
    for m = 1:M
        for u = 1:U
            wcell{m,u} = sqrt(cfg.power.Pmax/U)*directions{m,u};
        end
    end
    info.method = 'SO-BF';
    info.direction_type = 'equation_36_target_steering';
    info.power_allocation = 'equal_power_fallback_without_channels';
    return;
end

opts.method_name = 'SO-BF';
[wcell,info] = optimize_fixed_direction_power_allocation( ...
    cfg,channels,smat,directions,opts);
info.direction_type = 'equation_36_target_steering';
info.exact_equation_36 = true;
end
