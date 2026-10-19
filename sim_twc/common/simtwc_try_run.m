function record = simtwc_try_run(key,p,capture_history,opts)
%SIMTWC_TRY_RUN Convert solver failures into explicit records.

if nargin<4; opts=struct(); end
try
    record=simtwc_timed_method(key,p,capture_history,opts);
    record.status='ok';
catch ME
    record=struct('key',key,'seed',p.seed,'status','failed', ...
        'message',ME.message,'identifier',ME.identifier,'feasible',false, ...
        'runtime_median_s',NaN,'elapsed_s',NaN, ...
        'metrics',struct('crlb',NaN),'kkt',struct('r_kkt',NaN), ...
        'history',struct());
    if isfield(p,'initialization')
        record.initialization_success=p.initialization.success;
        record.initialization_runtime_s=p.initialization.runtime_s;
        record.initializer_message=p.initialization.message;
    end
    fprintf('    FAILED %s: %s\n',key,ME.message);
end
end
