function cfg_stage = configure_costar_tr_stage(cfg)
%CONFIGURE_COSTAR_TR_STAGE Synchronize the fair paper iteration budget.
%
% cfg.alm is the budget ceiling shared by PFO-ALM and PTR-ALM. The
% safeguarded wrapper may impose a smaller second-stage refinement budget.
% COSTAR does not call this helper because it has no ALM Phase I.

required = {'phase1_tcg_start_outer','phase1_verbose'};
if ~isfield(cfg,'alm') || ...
        ~all(isfield(cfg.alm,{'max_outer','max_inner'}))
    error('The common iterative budget cfg.alm.max_outer/max_inner is incomplete.');
end
if ~isfield(cfg,'costar') || ~all(isfield(cfg.costar,required))
    error('PTR-ALM trust-region configuration is incomplete.');
end

cfg_stage = cfg;
cfg_stage.tr.max_outer = cfg.alm.max_outer;
cfg_stage.tr.max_inner = cfg.alm.max_inner;
cfg_stage.tr.tcg_start_outer = cfg.costar.phase1_tcg_start_outer;
cfg_stage.tr.verbose = cfg.costar.phase1_verbose;
end
