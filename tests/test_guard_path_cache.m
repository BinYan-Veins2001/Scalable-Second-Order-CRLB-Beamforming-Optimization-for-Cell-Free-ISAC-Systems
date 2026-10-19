function report=test_guard_path_cache()
%TEST_GUARD_PATH_CACHE Compare quartic cache with independent direct values.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
spec=struct('M',4,'U',2,'Nt',2,'Nr',2,'gamma_dB',10,'min_ap_distance',10);
p=simtwc_make_problem(spec,2026);
assert(p.initialization.success,p.initialization.message);
W=p.initialization.Wmat;
rng(781,'twister');
P=randn(size(W))+1i*randn(size(W)); P=.1*norm(W,'fro')*P/norm(P,'fro');
V=randn(size(W))+1i*randn(size(W)); V=.01*norm(W,'fro')*V/norm(V,'fro');
cache=costar_guard_path_cache(p.cfg,p.channels,p.smat,W,P,V);
times=[0,1e-7,1e-3,.05]; errors=zeros(numel(times),4);
for k=1:numel(times)
    t=times(k); DW=t*P+t*t*V; Wt=W+DW;
    [f,c,actual,noise,out]=costar_guard_path_evaluate(cache,t);
    [f_direct,direct]=costar_paper_value(p.cfg,p.smat,Wt);
    c_direct=costar_sosp_constraints(p.cfg,p.channels,Wt);
    [actual_direct,noise_direct]=costar_sosp_stable_decrease(p.cfg,p.smat,W,DW);
    errors(k,:)=[abs(f-f_direct),max(abs(c-c_direct)), ...
        abs(actual-actual_direct),abs(out.crlb-direct.crlb)];
    assert(out.identifiable==direct.identifiable,'Path domain mismatch.');
    assert(errors(k,1)<=5e-12*max(1,abs(f_direct)),'Path objective mismatch.');
    assert(errors(k,2)<=5e-12*max(1,max(abs(c_direct))),'Path constraint mismatch.');
    assert(errors(k,3)<=max(1e-14,5*(noise+noise_direct)), ...
        'Stable decrease mismatch.');
    assert(errors(k,4)<=5e-12*max(1,abs(direct.crlb)),'CRLB mismatch.');
end
% Check the coefficient of t matches the real-vector Euclidean derivative.
[~,G]=costar_sosp_constraints(p.cfg,p.channels,W);
pv=complex_matrix_to_real_vector(P);
assert(norm(cache.ccoef(:,2)-G'*pv,inf)<=5e-11*max(1,norm(G'*pv,inf)), ...
    'The real-coordinate path derivative has a factor/conjugation error.');
% A zero path must return exactly zero increments, not subtract close logs.
zero=costar_guard_path_cache(p.cfg,p.channels,p.smat,W,zeros(size(W)),zeros(size(W)));
[~,~,a0,n0]=costar_guard_path_evaluate(zero,1);
assert(a0==0 && n0==0,'Zero path produces nonzero decrease/noise.');
report=struct('passed',true,'times',times,'absolute_errors',errors, ...
    'description','Quartic sensing/constraints, stable decrease, and real convention');
save(fullfile(root,'tests','guard_path_cache_report.mat'),'report');
fprintf('GUARD PATH CACHE passed: max objective %.3g, constraint %.3g, decrease %.3g.\n', ...
    max(errors(:,1)),max(errors(:,2)),max(errors(:,3)));
end
