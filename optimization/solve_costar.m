function [sol,hist]=solve_costar(cfg,channels,smat,initial_wcell,options)
%SOLVE_COSTAR One safeguarded loop with feasible first-order acceleration.
if nargin<4 || isempty(initial_wcell)
    error('COSTAR:ExternalInitializationRequired','Supply the shared feasible stacked-ZF point.');
end
if nargin<5; options=struct(); end
o=struct();
if isfield(cfg.costar,'paper'); o=cfg.costar.paper; end
names=fieldnames(options);
for j=1:numel(names)
    if ~strcmp(names{j},'terminal_options'); o.(names{j})=options.(names{j}); end
end
if isfield(options,'terminal_options')
    names=fieldnames(options.terminal_options);
    for j=1:numel(names); o.(names{j})=options.terminal_options.(names{j}); end
end
o=costar_paper_options(o);
if ~strcmpi(cfg.alm.objective_type,'log_crlb_ratio')
    error('COSTAR:PaperObjective','The paper-aligned solver uses log_crlb_ratio.');
end
[sol,trace]=costar_paper_rtr(cfg,channels,smat,initial_wcell,o);
sol.method='COSTAR'; sol.algorithm_version='COSTAR_Safeguard_Adaptive_20261019';
sol.phase2_method='Single safeguarded exact-HVP RTR with deterministic SPSC';
sol.initialization_runtime_s=0; sol.fast_runtime_s=0;
sol.terminal_runtime_s=sol.runtime; sol.fast_iterations=0;
sol.terminal_iterations=trace.iter(end); sol.terminal_options=o; sol.options=o;
sol.work_fast=struct(); sol.work_terminal=sol.work;
hist=struct('method','COSTAR','algorithm_version',sol.algorithm_version, ...
    'rtr',trace,'terminal',trace,'fast',struct('iter',[]), ...
    'terminal_start_time_s',0,'work',sol.work);
sol.hist=hist;
end
