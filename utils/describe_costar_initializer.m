function meta = describe_costar_initializer(initializer_mode,cfg)
%DESCRIBE_COSTAR_INITIALIZER Canonical name and traceable COSTAR metadata.

configured_mode = lower(strtrim(string(initializer_mode)));
mode = configured_mode;
if mode == "auto"
    if nargin < 2 || ~isfield(cfg,'sys') || ...
            ~isfield(cfg.sys,'M') || ~isfield(cfg.sys,'Nt') || ...
            ~isfield(cfg.sys,'U')
        error(['Resolving COSTAR initializer_mode=''auto'' requires ' ...
            'cfg.sys.M, cfg.sys.Nt, and cfg.sys.U.']);
    end
    if cfg.sys.M*cfg.sys.Nt < cfg.sys.U
        error(['The default global aggregate-ZF entry requires MNt >= U; ' ...
            'received M=%d, Nt=%d, U=%d. Select a legacy initializer ' ...
            'explicitly only for an ablation study.'], ...
            cfg.sys.M,cfg.sys.Nt,cfg.sys.U);
    end
    mode = "global_aggregate_zf";
end

switch mode
    case {"global_aggregate_zf","global_zf","network_zf", ...
            "aggregate_zf"}
        meta.mode = 'global_aggregate_zf';
        meta.name = 'Global aggregate theorem-aware ZF';
        refine_iter = 0;
        if nargin >= 2
            refine_iter = get_global_refine_iter(cfg);
        end
        if refine_iter > 0
            meta.description = [ ...
                'global aggregate exact ZF: pseudoinverse strict-SINR ' ...
                'floor + complex Dykstra AP-power projection + projected ' ...
                'CRLB redistribution (no CVX)'];
            meta.version = sprintf( ...
                'global_aggregate_zf_refine%d_as_rtr_v2',refine_iter);
            meta.policy_token = sprintf( ...
                'GLOBAL_aggregate_ZF_refine%d_no_CVX',refine_iter);
        else
            meta.description = [ ...
                'global aggregate exact ZF: pseudoinverse strict-SINR ' ...
                'floor + complex Dykstra AP-power feasibility projection ' ...
                '(no CRLB pre-refinement, no CVX)'];
            meta.version = 'global_aggregate_zf_feasible_entry_as_rtr_v2';
            meta.policy_token = 'GLOBAL_aggregate_ZF_feasible_entry_no_CVX';
        end

    case {"low_complexity_theorem_aware","theorem_aware", ...
            "low_complexity","lc_theorem_aware"}
        meta.mode = 'low_complexity_theorem_aware';
        meta.name = 'LC theorem-aware ZF';
        meta.description = [ ...
            'low-complexity theorem-aware ZF: closed-form strict-SINR ' ...
            'power floor + Dykstra projection + projected CRLB ' ...
            'redistribution (no CVX)'];
        meta.version = 'lc_theorem_aware_zf_as_rtr_v1';
        meta.policy_token = 'LC_theorem_aware_ZF_no_CVX';

    case {"zf_bf","zf","eq35_zf"}
        meta.mode = 'zf_bf';
        meta.name = 'ZF-BF';
        meta.description = ...
            'equation-(35) ZF-BF with CVX/SCA power allocation';
        meta.version = 'eq35_power_as_rtr_v2';
        meta.policy_token = 'EQ35_ZF_power_alloc';

    case {"random_rank_one","random","rank_one"}
        meta.mode = 'random_rank_one';
        meta.name = 'Random rank-one';
        meta.description = ...
            'random rank-one diagnostic point (no feasibility restoration)';
        meta.version = 'random_rank_one_as_rtr_v1';
        meta.policy_token = 'random_rank_one_direct_AS_RTR';

    otherwise
        error('Unknown cfg.costar.initializer_mode: %s',mode);
end
meta.initializer_version = meta.version;
if nargin >= 2 && isfield(cfg,'costar') && ...
        isfield(cfg.costar,'phase2_version')
    meta.version = sprintf('%s__%s',meta.initializer_version, ...
        cfg.costar.phase2_version);
end
meta.configured_mode = char(configured_mode);
meta.was_auto_selected = configured_mode == "auto";
end

function refine_iter = get_global_refine_iter(cfg)
refine_iter = 0;
if nargin >= 1 && isfield(cfg,'costar') && ...
        isfield(cfg.costar,'global_zf_init') && ...
        isfield(cfg.costar.global_zf_init,'max_refine_iter')
    refine_iter = cfg.costar.global_zf_init.max_refine_iter;
end
end
