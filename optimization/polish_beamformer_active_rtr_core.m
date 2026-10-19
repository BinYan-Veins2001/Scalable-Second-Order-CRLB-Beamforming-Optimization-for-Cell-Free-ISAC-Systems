function [sol,hist] = polish_beamformer_active_rtr_core( ...
    cfg,channels,smat,init_wcell,opts)
%POLISH_BEAMFORMER_ACTIVE_RTR_CORE Simplified, SOSP-aware AS-RTR phase II.
%
% The core variant keeps active-set geometry, NNLS multipliers, tCG,
% feasibility retraction and exact trust-region acceptance. It removes
% duplicate HVP/trial evaluations and replaces periodic Lanczos sweeps by a
% second-order certificate triggered at a genuine first-order candidate.

if nargin < 5
    opts = struct();
end
opts = set_default(opts,'curvature_policy','terminal');
opts = set_default(opts,'reuse_tcg_model',true);
opts = set_default(opts,'fuse_trial_evaluation',true);
opts = set_default(opts,'enable_sinr_correction',true);
opts = set_default(opts,'multiplier_mode','nnls');
opts = set_default(opts,'cg_max_iter',30);
opts = set_default(opts,'sinr_retraction_max_iter',4);
opts = set_default(opts,'lanczos_max_iter',20);
opts = set_default(opts,'hessian_model','exact_lagrangian');
opts = set_default(opts,'enable_tcg_negative_curvature',true);
% The formal audit found no accepted explicit Ritz trial. Keep Lanczos as
% a terminal certificate and leave the dormant escape branch opt-in.
opts = set_default(opts,'enable_ritz_escape',false);
[sol,hist] = polish_beamformer_active_rtr_simplified( ...
    cfg,channels,smat,init_wcell,opts);
end

function s = set_default(s,name,value)
if ~isfield(s,name)
    s.(name) = value;
end
end
