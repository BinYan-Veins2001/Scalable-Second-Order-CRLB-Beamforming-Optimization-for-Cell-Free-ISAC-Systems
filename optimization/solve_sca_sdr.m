function [sol,hist] = solve_sca_sdr(cfg,channels,smat,init_wcell)
%SOLVE_SCA_SDR Rank-continuation SDR-SCA benchmark.
%
% This remains a lifted SDR-SCA method: no COSTAR, RTR, MRFP-SCA, or vector
% post-polishing is used. Independent PSD blocks are exactly equivalent to
% the former block-diagonal cone. A decreasing affine rank-leakage budget
% avoids the CRLB loss caused by the previous fixed rank-penalty jump.

if isfield(cfg,'sdr')
    opts = cfg.sdr;
else
    opts = struct();
end
opts.cvx_quiet = cfg.opt.cvx_quiet;
opts.cvx_solver = cfg.opt.cvx_solver;
opts.verbose = cfg.opt.verbose;
[sol,hist] = solve_sca_sdr_rank_continuation( ...
    cfg,channels,smat,init_wcell,opts);
if ~logical(sol.success)
    % Never expose the entry or an uncertified covariance as an SDR result.
    sol.Wcell = {};
end
end
