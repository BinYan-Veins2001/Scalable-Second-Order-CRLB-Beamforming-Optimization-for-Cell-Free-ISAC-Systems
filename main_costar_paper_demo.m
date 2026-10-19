function result=main_costar_paper_demo(seed,options)
%MAIN_COSTAR_PAPER_DEMO Run the paper-aligned solver from shared stacked ZF.
if nargin<1; seed=2026; end
if nargin<2; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
s=simtwc_setup(); p=simtwc_make_problem(s.default,seed);
assert(p.initialization.success,p.initialization.message);
[sol,hist]=solve_costar(p.cfg,p.channels,p.smat,p.initialization.wcell,options);
audit=benchmark_endpoint_audit(p,sol.z,struct('lambda',sol.lambda,'active',sol.active));
result=struct('problem',p,'solution',sol,'history',hist,'audit',audit);
folder=fullfile(root,'results','demo'); if ~exist(folder,'dir'); mkdir(folder); end
save(fullfile(folder,sprintf('seed%d.mat',seed)),'result','-v7.3');
fprintf('COSTAR paper M=%d U=%d Nt=%d: CRLB=%.10g, KKT=%.3g, eigT=%.3g, time=%.3fs\n', ...
    p.cfg.sys.M,p.cfg.sys.U,p.cfg.sys.Nt,sol.metrics.crlb,audit.kkt, ...
    audit.lambda_min_tangent,sol.runtime);
fprintf('Numerical joint certificate: %d. %s\n',sol.success,sol.message);
end
