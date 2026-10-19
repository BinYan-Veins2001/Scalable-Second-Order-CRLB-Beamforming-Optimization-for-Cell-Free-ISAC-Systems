function report=test_latest_alignment()
%TEST_LATEST_ALIGNMENT Mechanism checks for the current single-loop solver.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root),'-begin');
setenv('COSTAR_QUICK',''); setenv('COSTAR_SMOKE','');
maxNumCompThreads(1);
assert(startsWith(which('solve_costar'),root));
o=costar_paper_options();
assert(~o.guard_directional && ~o.guard_scaled_inward && ~o.ambient_endpoint_probe);
assert(o.guard_cached && o.guard_initial_scaling && o.accelerator_norm_cap);
assert(o.epsilon_KKT==1e-5 && o.epsilon_H==1e-5 && o.feas_tol==1e-10);
report=struct('passed',false,'solver_source',which('solve_costar'), ...
    'matlab_version',version,'threads',1,'defaults',o);
report.geometry=test_paper_geometry();
report.spsc=test_paper_spsc();
report.stopping=test_paper_stopping();
report.path_cache=test_guard_path_cache();
report.physical=test_paper_alignment();
assert(report.physical.solver_success,'Physical solver did not certify endpoint.');
assert(report.physical.work.fd_hvp_calls==0);
assert(report.physical.work.meo_calls==0);
assert(report.physical.work.meo_hvp_calls==0);
report.passed=true;
save(fullfile(root,'audit','alignment_mechanisms.mat'),'report','-v7.3');
fid=fopen(fullfile(root,'audit','alignment_mechanisms.json'),'w','n','UTF-8');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('LATEST ALIGNMENT: 5 mechanism suites passed.\n');
end
