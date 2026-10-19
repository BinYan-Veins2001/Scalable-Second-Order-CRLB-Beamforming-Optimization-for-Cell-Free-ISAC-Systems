function curve = evaluate_second_order_stationarity_history( ...
    cfg, channels, smat, name, snapshot_iter, z_snapshots, opts)
%EVALUATE_SECOND_ORDER_STATIONARITY_HISTORY Diagnose selected iterates.

if nargin < 7
    opts = struct();
end
if ~isfield(opts,'max_points'); opts.max_points = 21; end
if ~isfield(opts,'seed'); opts.seed = 2026; end

snapshot_iter = snapshot_iter(:).';
if isempty(snapshot_iter) || isempty(z_snapshots)
    error('No stored iterates are available for %s.', name);
end
if size(z_snapshots,2) ~= numel(snapshot_iter)
    error('Snapshot iteration and state counts do not match for %s.', name);
end

[snapshot_iter, unique_idx] = unique(snapshot_iter, 'stable');
z_snapshots = z_snapshots(:,unique_idx);
n = numel(snapshot_iter);
keep = unique(round(linspace(1, n, min(opts.max_points,n))));

curve.name = name;
curve.iter = snapshot_iter(keep);
curve.sosp_residual = zeros(size(keep));
curve.kkt_residual = zeros(size(keep));
curve.negative_curvature = zeros(size(keep));
curve.normalized_negative_curvature = zeros(size(keep));
curve.lambda_min = zeros(size(keep));
curve.primal_residual = zeros(size(keep));
curve.stationarity_residual = zeros(size(keep));
curve.complementarity_residual = zeros(size(keep));
curve.active_count = zeros(size(keep));

for k = 1:numel(keep)
    opts_k = opts;
    opts_k.seed = opts.seed + k - 1;
    d = evaluate_second_order_stationarity( ...
        cfg, channels, smat, z_snapshots(:,keep(k)), opts_k);
    curve.sosp_residual(k) = d.sosp_residual;
    curve.kkt_residual(k) = d.kkt_residual;
    curve.negative_curvature(k) = d.negative_curvature;
    curve.normalized_negative_curvature(k) = d.normalized_negative_curvature;
    curve.lambda_min(k) = d.lambda_min_weak_critical;
    curve.primal_residual(k) = d.primal_residual;
    curve.stationarity_residual(k) = d.stationarity_residual;
    curve.complementarity_residual(k) = d.complementarity_residual;
    curve.active_count(k) = d.active_count;
end
end
