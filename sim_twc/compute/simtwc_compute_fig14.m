function archive = simtwc_compute_fig14()
%SIMTWC_COMPUTE_FIG14 Separate user-count and antenna-count scaling tests.

s=simtwc_setup(); spec=s.fig14; methods=cellstr(string(spec.methods));
data=simtwc_initialize_data(14,spec,methods); path=simtwc_archive_path(14);
nA=numel(methods); nT=numel(spec.seeds);

nU=numel(spec.U_values);
if ~isfield(data,'user_runs') || ~same_size3(data.user_runs,nA,nU,nT)
    data.user_runs=cell(nA,nU,nT);
end
for u=1:nU
    pspec=spec.base; pspec.U=spec.U_values(u);
    for t=1:nT
        p=simtwc_make_problem(pspec,spec.seeds(t));
        for a=1:nA
            if ok(data.user_runs{a,u,t}); continue; end
            fprintf('[Fig14 users] U=%d, seed=%d, %s\n', ...
                pspec.U,spec.seeds(t),methods{a});
            data.user_runs{a,u,t}=simtwc_try_run(methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end

nN=numel(spec.Nt_values);
if ~isfield(data,'antenna_runs') || ~same_size3(data.antenna_runs,nA,nN,nT)
    data.antenna_runs=cell(nA,nN,nT);
end
for n=1:nN
    pspec=spec.base; pspec.Nt=spec.Nt_values(n); pspec.Nr=spec.Nt_values(n);
    for t=1:nT
        p=simtwc_make_problem(pspec,spec.seeds(t));
        for a=1:nA
            if ok(data.antenna_runs{a,n,t}); continue; end
            fprintf('[Fig14 antennas] Nt=Nr=%d, seed=%d, %s\n', ...
                pspec.Nt,spec.seeds(t),methods{a});
            data.antenna_runs{a,n,t}=simtwc_try_run( ...
                methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end

data.methods=methods; data.seeds=spec.seeds;
data.U_values=spec.U_values; data.Nt_values=spec.Nt_values;
data.user_dimension=2*spec.base.M*spec.base.Nt*spec.U_values;
data.antenna_dimension=2*spec.base.M*spec.base.U*spec.Nt_values;
[data.user_crlb,data.user_runtime_s,data.user_feasible, ...
    data.user_kkt,data.user_hvp_calls,data.user_tcg_iterations]= ...
    extract(data.user_runs);
[data.antenna_crlb,data.antenna_runtime_s,data.antenna_feasible, ...
    data.antenna_kkt,data.antenna_hvp_calls,data.antenna_tcg_iterations]= ...
    extract(data.antenna_runs);
data.complete=all(cellfun(@ok,data.user_runs(:))) && ...
    all(cellfun(@ok,data.antenna_runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function [C,T,F,K,H,I]=extract(runs)
sz=size(runs); C=nan(sz); T=nan(sz); F=false(sz); K=nan(sz);
H=nan(sz); I=nan(sz);
for k=1:numel(runs)
    q=runs{k}; if ~ok(q); continue; end
    C(k)=q.metrics.crlb; T(k)=q.runtime_median_s;
    F(k)=q.feasible; K(k)=q.kkt.r_kkt;
    if isfield(q,'solution') && isfield(q.solution,'work')
        w=q.solution.work;
        if isfield(w,'hvp_calls'); H(k)=w.hvp_calls; end
        if isfield(w,'tcg_iterations'); I(k)=w.tcg_iterations; end
    end
end
end

function yes=same_size3(x,n1,n2,n3)
yes=size(x,1)==n1 && size(x,2)==n2 && size(x,3)==n3;
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
