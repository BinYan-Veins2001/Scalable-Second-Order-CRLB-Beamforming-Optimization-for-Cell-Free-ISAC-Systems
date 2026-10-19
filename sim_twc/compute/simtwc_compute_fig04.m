function archive = simtwc_compute_fig04()
%SIMTWC_COMPUTE_FIG04 Derive equal-target and matched-KKT timing from Fig. 3.

s=simtwc_setup(); spec=s.fig04; methods=s.iterative_methods;
path3=simtwc_archive_path(3);
source=simtwc_load_data(3);
data=simtwc_initialize_data(4,spec,methods); path=simtwc_archive_path(4);
all_methods=source.methods; seeds=source.seeds;
nA=numel(methods); nT=numel(seeds);
time_target=nan(nA,nT); time_kkt=nan(nA,nT);
gap=nan(nA,nT); feasible=false(nA,nT);
target=nan(1,nT);
isdr=find(strcmp(all_methods,'sdr_sca'),1);
for t=1:nT
    finals=nan(1,numel(all_methods)); valid=false(size(finals));
    for j=1:numel(all_methods)
        r=source.runs{j,t};
        if ok(r) && r.feasible
            finals(j)=r.metrics.crlb; valid(j)=true;
        end
    end
    if any(valid)
        target(t)=spec.crlb_factor*min(finals(valid));
    end
    for a=1:nA
        j=find(strcmp(all_methods,methods{a}),1); r=source.runs{j,t};
        if ~ok(r); continue; end
        feasible(a,t)=r.feasible;
        time_target(a,t)=simtwc_first_hit_time(r.history,'crlb',target(t));
        time_kkt(a,t)=simtwc_first_hit_time(r.history,'kkt',spec.kkt_tol);
        rsdr=source.runs{isdr,t};
        if r.feasible && ok(rsdr) && rsdr.feasible
            gap(a,t)=100*(r.metrics.crlb/rsdr.metrics.crlb-1);
        end
    end
end
data.seeds=seeds; data.source_fig03=path3;
data.time_to_target_s=time_target; data.matched_kkt_time_s=time_kkt;
data.paired_crlb_gap_pct=gap; data.feasible=feasible;
data.crlb_target=target;
data.target_reach_rate=mean(isfinite(time_target),2);
data.kkt_reach_rate=mean(isfinite(time_kkt),2);
data.feasibility_rate=mean(feasible,2);
data.complete=true; data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function yes=ok(r)
yes=isstruct(r)&&isfield(r,'status')&&strcmp(r.status,'ok');
end
