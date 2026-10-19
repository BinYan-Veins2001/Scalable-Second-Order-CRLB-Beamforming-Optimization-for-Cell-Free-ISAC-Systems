function result = simtwc_timed_method(key,p,capture_history,run_opts)
%SIMTWC_TIMED_METHOD End-to-end runtime with a history-bearing primary run.

if nargin<3; capture_history=true; end
if nargin<4; run_opts=struct(); end
s=simtwc_config();
run_opts.capture_history=capture_history;
result=simtwc_run_method(key,p,run_opts);
runtimes=result.elapsed_s;
optimizer_times=result.optimization_runtime_s;
for rep=2:s.timing_repetitions
    qopts=run_opts; qopts.capture_history=false;
    q=simtwc_run_method(key,p,qopts);
    runtimes(end+1)=q.elapsed_s; %#ok<AGROW>
    optimizer_times(end+1)=q.optimization_runtime_s; %#ok<AGROW>
end
result.runtime_repetitions_s=runtimes;
result.runtime_median_s=median(runtimes,'omitnan');
result.optimization_repetitions_s=optimizer_times;
result.optimization_median_s=median(optimizer_times,'omitnan');
result.runtime_definition=sprintf( ...
    'median of %d (shared stacked-ZF initialization + optimization) run(s)',numel(runtimes));
end
