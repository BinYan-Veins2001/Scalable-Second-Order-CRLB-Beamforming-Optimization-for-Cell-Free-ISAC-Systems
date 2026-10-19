function archive = simtwc_compute_fig09()
%SIMTWC_COMPUTE_FIG09 Power-budget sweep and COSTAR active-set composition.

s=simtwc_setup(); spec=s.fig09; methods=spec.methods;
data=simtwc_initialize_data(9,spec,methods); path=simtwc_archive_path(9);
nA=numel(methods); nP=numel(spec.power_dBm); nT=numel(spec.seeds);
if ~isfield(data,'runs') || ~same_size3(data.runs,nA,nP,nT)
    data.runs=cell(nA,nP,nT);
end
active_sinr=nan(nP,nT); active_power=nan(nP,nT);
if isfield(data,'active_sinr'); active_sinr=data.active_sinr; end
if isfield(data,'active_power'); active_power=data.active_power; end
for q=1:nP
    problem_spec=spec; problem_spec.power_dBm=spec.power_dBm(q);
    for t=1:nT
        p=simtwc_make_problem(problem_spec,spec.seeds(t));
        for a=1:nA
            if ok(data.runs{a,q,t}); continue; end
            fprintf('[Fig09] Pmax=%g dBm, seed=%d, %s\n', ...
                spec.power_dBm(q),spec.seeds(t),methods{a});
            data.runs{a,q,t}=simtwc_try_run(methods{a},p,false,struct());
            if strcmp(methods{a},'costar') && ok(data.runs{a,q,t})
                active_sinr(q,t)=data.runs{a,q,t}.kkt.active_sinr_count;
                active_power(q,t)=data.runs{a,q,t}.kkt.active_power_count;
            end
            data.active_sinr=active_sinr; data.active_power=active_power;
            simtwc_save_atomic(path,data);
        end
    end
end
data.power_dBm=spec.power_dBm; data.seeds=spec.seeds;
data.active_sinr=active_sinr; data.active_power=active_power;
[data.final_crlb,data.runtime_s,data.feasible]=extract_results(data.runs);
data.complete=all(cellfun(@ok,data.runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function yes=ok(x)
yes=isstruct(x)&&isfield(x,'status')&&strcmp(x.status,'ok');
end

function yes=same_size3(x,n1,n2,n3)
yes=size(x,1)==n1 && size(x,2)==n2 && size(x,3)==n3;
end

function [C,T,F]=extract_results(runs)
sz=size(runs); C=nan(sz); T=nan(sz); F=false(sz);
for k=1:numel(runs)
    q=runs{k}; if ~ok(q); continue; end
    C(k)=q.metrics.crlb; T(k)=q.runtime_median_s; F(k)=q.feasible;
end
end
