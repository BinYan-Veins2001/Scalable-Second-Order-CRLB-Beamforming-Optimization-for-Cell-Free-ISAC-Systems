function [wcell,info] = build_shared_theorem_aware_zf_initializer( ...
    cfg,channels,smat)
%BUILD_SHARED_THEOREM_AWARE_ZF_INITIALIZER Common paper-baseline entry.
%
% PFO-ALM, PPFO-ALM, and PTR-ALM must receive the same beamformer in each
% realization. The default paper policy is centralized aggregate exact ZF
% whenever MNt >= U. Local ZF/RZF modes are explicit legacy ablations and
% are never selected silently by the default policy.

if isfield(cfg,'costar') && isfield(cfg.costar,'initializer_mode')
    configured_mode = cfg.costar.initializer_mode;
else
    configured_mode = 'auto';
end
meta = describe_costar_initializer(configured_mode,cfg);

timer = tic;
switch string(meta.mode)
    case "global_aggregate_zf"
        if isfield(cfg.costar,'global_zf_init')
            opts = cfg.costar.global_zf_init;
        else
            opts = struct();
        end
        [wcell,raw_info] = build_global_aggregate_zf_initializer( ...
            cfg,channels,smat,opts);
    case "low_complexity_theorem_aware"
        if isfield(cfg.costar,'low_complexity_init')
            opts = cfg.costar.low_complexity_init;
        else
            opts = struct();
        end
        [wcell,raw_info] = build_low_complexity_theorem_aware_initializer( ...
            cfg,channels,smat,opts);
    case "zf_bf"
        [wcell,raw_info] = build_zf_beamformer(cfg,channels,smat);
    otherwise
        error(['Paper-facing ALM baselines require a theorem-aware ' ...
            'ZF-family entry, not initializer mode %s.'],meta.mode);
end

Wmat = vectorize_wcell_to_matrix(cfg,wcell);
metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
entry_feasible = all(metrics.sinr >= cfg.comm.Gamma*(1-1e-4)) && ...
    all(metrics.ap_power <= cfg.power.Pmax*(1+1e-4));
if ~entry_feasible
    error(['The shared %s entry failed its original-constraint check: ' ...
        'min SINR %.3f dB, max AP power %.6g W.'],meta.mode, ...
        pow2db(max(min(metrics.sinr),realmin)),max(metrics.ap_power));
end

info = raw_info;
info.method = meta.name;
info.configured_mode = meta.configured_mode;
info.resolved_mode = meta.mode;
info.policy_token = ['shared_ALM_' meta.policy_token];
info.initializer_version = meta.initializer_version;
info.entry_feasible = entry_feasible;
info.metrics = metrics;
info.total_runtime = toc(timer);
info.initial_Wmat = Wmat;
info.fingerprint = [norm(Wmat,'fro'),sum(real(Wmat(:))), ...
    sum(imag(Wmat(:))),sum(abs(Wmat(:)).^2)];
end
