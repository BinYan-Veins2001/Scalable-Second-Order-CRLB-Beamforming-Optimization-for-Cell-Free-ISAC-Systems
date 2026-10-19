function scenario = generate_scenario(cfg)
%GENERATE_SCENARIO Generate AP, user, and target positions.

M = cfg.sys.M;
U = cfg.sys.U;

if cfg.area.use_paper_like_ap_layout && M == 4
    scenario.ap_pos = [50, 0; -50, 0; 0, 50*sqrt(3); 0, -50*sqrt(3)];
else
    scenario.ap_pos = sample_points_with_min_distance( ...
        M, cfg.area.ap_x, cfg.area.ap_y, cfg.area.min_ap_distance, []);
end

scenario.user_pos = sample_points_with_min_distance( ...
    U, cfg.area.user_x, cfg.area.user_y, cfg.area.min_user_ap_distance, scenario.ap_pos);

scenario.target_pos = sample_points_with_min_distance( ...
    1, cfg.area.target_x, cfg.area.target_y, cfg.area.min_target_ap_distance, scenario.ap_pos);

scenario.M = M;
scenario.U = U;
end

function pts = sample_points_with_min_distance(K, xrange, yrange, min_dist, ref_pts)
pts = zeros(K,2);
max_trials = 10000;

for k = 1:K
    accepted = false;
    for t = 1:max_trials
        cand = [xrange(1) + diff(xrange)*rand, yrange(1) + diff(yrange)*rand];
        ok_self = true;
        if k > 1
            d_self = sqrt(sum((pts(1:k-1,:) - cand).^2, 2));
            ok_self = all(d_self >= min_dist);
        end
        ok_ref = true;
        if ~isempty(ref_pts)
            d_ref = sqrt(sum((ref_pts - cand).^2, 2));
            ok_ref = all(d_ref >= min_dist);
        end
        if ok_self && ok_ref
            pts(k,:) = cand;
            accepted = true;
            break;
        end
    end
    if ~accepted
        error('Failed to sample positions. Relax minimum distance or enlarge region.');
    end
end
end
