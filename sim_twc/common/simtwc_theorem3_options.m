function o = simtwc_theorem3_options(options)
%SIMTWC_THEOREM3_OPTIONS User-editable Theorem 3 experiment specification.
if nargin<1; options=struct(); end
s=simtwc_config();
o=struct('mc_trials',10,'base_seed',2026,'seeds',[], ...
    'regimes',s.fig06.regimes,'epsilon_KKT',1e-5,'epsilon_H',1e-5, ...
    'failure_probability',1e-3,'terminal_max_iter',200, ...
    'use_structured_psd',true,'curvature_stride',8, ...
    'cone_samples',128,'random_trials',30,'random_audit_budget',32, ...
    'oracle_trace_max_iter',inf,'dense_max_dim',4096, ...
    'resume',s.resume,'output_tag','');
if s.quick
    o.mc_trials=1; o.terminal_max_iter=12; o.random_trials=3;
    o.random_audit_budget=16; o.oracle_trace_max_iter=24;
    o.curvature_stride=12; o.cone_samples=16;
end
names=fieldnames(options);
for j=1:numel(names)
    assert(isfield(o,names{j}),'Unknown Fig. 20 option: %s',names{j});
    o.(names{j})=options.(names{j});
end
if isempty(o.seeds); o.seeds=o.base_seed+(0:o.mc_trials-1); end
o.seeds=o.seeds(:).'; o.mc_trials=numel(o.seeds);
assert(o.mc_trials>=1 && all(o.seeds>=0) && all(mod(o.seeds,1)==0));
assert(o.epsilon_KKT>0 && o.epsilon_H>0);
assert(o.failure_probability>0 && o.failure_probability<1);
assert(o.curvature_stride>=1 && o.random_trials>=1);
assert(o.terminal_max_iter>=1 && o.random_audit_budget>=1);
assert(isempty(o.output_tag) || ...
    ~isempty(regexp(o.output_tag,'^[A-Za-z0-9_-]+$','once')));
end
