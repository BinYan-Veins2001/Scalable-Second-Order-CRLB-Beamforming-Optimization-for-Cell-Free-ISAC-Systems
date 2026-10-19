function target_info = compute_target_geometry_for_ula(scenario,reference_ap)
%COMPUTE_TARGET_GEOMETRY_FOR_ULA Target geometry seen from one reference AP.

if ~isfield(scenario,'ap_pos') || ~isfield(scenario,'target_pos')
    error('Scenario must contain ap_pos and target_pos.');
end
if ~isscalar(reference_ap) || reference_ap ~= round(reference_ap) || ...
        reference_ap < 1 || reference_ap > size(scenario.ap_pos,1)
    error('reference_ap must index one row of scenario.ap_pos.');
end
if numel(scenario.target_pos) ~= 2 || ...
        any(~isfinite(scenario.target_pos(:)))
    error('scenario.target_pos must be one finite two-dimensional point.');
end

ap_position = reshape(scenario.ap_pos(reference_ap,:),1,2);
target_position = reshape(scenario.target_pos,1,2);
delta = target_position-ap_position;
distance_m = hypot(delta(1),delta(2));
if distance_m <= sqrt(eps)
    error('Target and reference AP positions must be distinct.');
end

physical_bearing_rad = atan2(delta(2),delta(1));
% A ULA with steering phase sin(theta) cannot distinguish front/back
% bearings. Map the physical bearing to its visible theta in [-pi/2,pi/2].
sin_bearing = min(max(sin(physical_bearing_rad),-1),1);
visible_theta_rad = asin(sin_bearing);
plot_direction_deg = rad2deg(visible_theta_rad)+90;
plot_direction_deg = min(max(plot_direction_deg,0),180);

target_info = struct( ...
    'reference_ap',reference_ap, ...
    'reference_ap_position_m',ap_position, ...
    'target_position_m',target_position, ...
    'delta_m',delta, ...
    'distance_m',distance_m, ...
    'physical_bearing_rad',physical_bearing_rad, ...
    'physical_bearing_deg',rad2deg(physical_bearing_rad), ...
    'visible_theta_rad',visible_theta_rad, ...
    'visible_theta_deg',rad2deg(visible_theta_rad), ...
    'plot_direction_deg',plot_direction_deg);
end
