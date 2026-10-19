function t = simtwc_first_hit_time(history,metric,threshold)
%SIMTWC_FIRST_HIT_TIME First feasible CRLB or KKT threshold crossing.

t=NaN;
switch lower(metric)
    case 'crlb'
        hit=find(history.feasible & history.crlb<=threshold,1,'first');
    case 'kkt'
        hit=find(history.kkt_residual<=threshold,1,'first');
    otherwise
        error('Unknown first-hit metric: %s',metric);
end
if ~isempty(hit); t=history.time_s(hit); end
end
