function o=paperfig_options(options)
%PAPERFIG_OPTIONS Small, shared controls for resumable paper experiments.
if nargin<1; options=struct(); end
s=paperfig_setup();
o=struct('seeds',s.seeds,'resume',true,'export',true, ...
    'solver_options',struct(),'dense_max_dim',4096, ...
    'lanczos_max_iter',40,'lanczos_tol',1e-9, ...
    'epsilon_KKT',s.epsilon_KKT,'epsilon_H',s.epsilon_H);
names=fieldnames(options);
for j=1:numel(names); o.(names{j})=options.(names{j}); end
o.seeds=o.seeds(:).';
assert(~isempty(o.seeds) && all(o.seeds>=0 & mod(o.seeds,1)==0));
end
