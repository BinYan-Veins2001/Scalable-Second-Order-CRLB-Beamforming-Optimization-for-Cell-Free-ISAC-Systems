function report=test_costar_regimes()
%TEST_COSTAR_REGIMES Full settings, independent entries, several network regimes.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
quick=getenv('COSTAR_QUICK'); smoke=getenv('COSTAR_SMOKE');
cleanup=onCleanup(@() restore(quick,smoke)); %#ok<NASGU>
setenv('COSTAR_QUICK',''); setenv('COSTAR_SMOKE','');
s=simtwc_setup(); report=cell(1,4);
seeds=[2027 2035 2026 2026]; M=[16 16 16 32]; Gamma=[10 10 20 10];
for j=1:numel(seeds)
    spec=s.default; spec.M=M(j); spec.gamma_dB=Gamma(j);
    p=simtwc_make_problem(spec,seeds(j));
    q=simtwc_run_method('costar',p); report{j}=q;
    assert(q.feasible && q.solution.work.fd_hvp_calls==0);
    assert(isequal(q.initial_Wmat,p.initialization.Wmat));
    assert(norm(q.history.raw.rtr.z_snapshots(:,end)- ...
        complex_matrix_to_real_vector(q.Wmat))==0);
    fprintf('REGIME seed=%d M=%d SINR=%g: CRLB %.8g KKT %.3g time %.3fs SOSP %d [%s]\n', ...
        seeds(j),M(j),Gamma(j),q.metrics.crlb,q.kkt.r_kkt,q.elapsed_s, ...
        q.conditional_sosp_reached,q.message);
end
save(fullfile(root,'tests','costar_regimes_report.mat'),'report','-v7.3');
end

function restore(quick,smoke)
setenv('COSTAR_QUICK',quick); setenv('COSTAR_SMOKE',smoke);
end
