function archive = simtwc_compute_fig07()
%SIMTWC_COMPUTE_FIG07 Seven-method target-SINR sweep.

s=simtwc_setup(); spec=s.fig07; methods=s.methods;
data=simtwc_initialize_data(7,spec,methods); path=simtwc_archive_path(7);
nA=numel(methods); nG=numel(spec.gamma_dB); nT=numel(spec.seeds);
if ~isfield(data,'runs') || ~same_size3(data.runs,nA,nG,nT)
    data.runs=cell(nA,nG,nT);
end
for g=1:nG
    q=spec; q.gamma_dB=spec.gamma_dB(g);
    for t=1:nT
        p=simtwc_make_problem(q,spec.seeds(t));
        for a=1:nA
            if ok(data.runs{a,g,t}); continue; end
            fprintf('[Fig07] SINR=%g dB, seed=%d, %s\n', ...
                q.gamma_dB,spec.seeds(t),methods{a});
            data.runs{a,g,t}=simtwc_try_run(methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end
data.gamma_dB=spec.gamma_dB; data.seeds=spec.seeds;
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
