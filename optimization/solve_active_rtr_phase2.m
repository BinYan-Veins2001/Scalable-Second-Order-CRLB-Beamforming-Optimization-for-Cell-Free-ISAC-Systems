function [sol,hist] = solve_active_rtr_phase2( ...
    cfg,channels,smat,init_wcell,rtr_opts)
%SOLVE_ACTIVE_RTR_PHASE2 Dispatch Full or Core active-set RTR.
%
% Direct AS-RTR baselines keep the established Full implementation unless
% phase2_variant is explicitly supplied. solve_costar sets the package
% default to the KKT-triggered Core implementation.

if nargin < 5
    rtr_opts = struct();
end
variant = "full";
if isfield(rtr_opts,'phase2_variant')
    variant = lower(strtrim(string(rtr_opts.phase2_variant)));
    rtr_opts = rmfield(rtr_opts,'phase2_variant');
end

switch variant
    case {"full","full_as_rtr","periodic"}
        [sol,hist] = polish_beamformer_active_rtr( ...
            cfg,channels,smat,init_wcell,rtr_opts);
        method = 'AS RTR';
        version = 'full_as_rtr_periodic_v1';
    case {"core","core_kkt_triggered","terminal"}
        if isfield(cfg,'costar') && isfield(cfg.costar,'core_rtr')
            rtr_opts = fill_missing(rtr_opts,cfg.costar.core_rtr);
        end
        [sol,hist] = polish_beamformer_active_rtr_core( ...
            cfg,channels,smat,init_wcell,rtr_opts);
        method = 'Core AS-RTR';
        if isfield(rtr_opts,'hessian_model') && any(strcmpi( ...
                rtr_opts.hessian_model,{'lagrangian_fd','full','full_lagrangian'}))
            version = 'core_as_rtr_fd_hvp_v1';
        else
            version = 'core_as_rtr_exact_hvp_v2';
        end
    otherwise
        error('Unknown active-RTR phase2_variant: %s',variant);
end

sol.phase2_method = method;
sol.phase2_variant = char(variant);
sol.phase2_version = version;
hist.phase2_method = method;
hist.phase2_variant = char(variant);
hist.phase2_version = version;
sol.hist = hist;
end

function target = fill_missing(target,defaults)
names = fieldnames(defaults);
for idx = 1:numel(names)
    name = names{idx};
    if ~isfield(target,name)
        target.(name) = defaults.(name);
    end
end
end
