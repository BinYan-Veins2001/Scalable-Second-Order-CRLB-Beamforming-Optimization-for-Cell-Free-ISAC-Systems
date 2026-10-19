function archive=simtwc_compute_fig11()
%SIMTWC_COMPUTE_FIG11 Postprocess Fig03/Fig07; no independent solver runs.
s=simtwc_setup(); d=simtwc_load_data(7); b=simtwc_load_data(3);
nA=numel(d.methods); nG=numel(d.gamma_dB); nT=numel(d.seeds);
G=nan(nA,nG*nT); refok=false(1,nG*nT);
ir=find(strcmp(d.methods,'sdr_sca'),1);
for g=1:nG
    for t=1:nT
        j=(g-1)*nT+t; ref=d.runs{ir,g,t};
        if ~valid(ref); continue; end
        refok(j)=true;
        for a=1:nA
            q=d.runs{a,g,t};
            if valid(q); G(a,j)=100*(q.metrics.crlb/ref.metrics.crlb-1); end
        end
    end
end
if ~any(refok); error('Simulation:EmptyFigure','No feasible paired SDR references.'); end
epsgrid=logspace(-4,4,180); success=zeros(nA,numel(epsgrid));
for a=1:nA
    for j=1:numel(epsgrid)
        success(a,j)=100*sum(G(a,refok)<=epsgrid(j))/sum(refok);
    end
end
times=inf(nA,numel(b.seeds)); defaultok=false(1,numel(b.seeds));
for t=1:numel(b.seeds)
    best=Inf;
    for a=1:numel(b.methods)
        q=b.runs{a,t}; if valid(q); best=min(best,q.metrics.crlb); end
    end
    if ~isfinite(best); continue; end
    defaultok(t)=true;
    for a=1:nA
        k=find(strcmp(b.methods,d.methods{a}),1); if isempty(k); continue; end
        q=b.runs{k,t}; if ~ok(q); continue; end
        if isfield(q.history,'time_s')
            h=q.history;
            hit=h.feasible(:)&isfinite(h.crlb(:))&h.crlb(:)>0& ...
                h.crlb(:)<=1.01*best&isfinite(h.time_s(:));
            if any(hit); times(a,t)=min(h.time_s(hit)); end
        end
        if valid(q) && q.metrics.crlb<=1.01*best
            times(a,t)=min(times(a,t),q.runtime_median_s);
        end
    end
end
if ~any(defaultok); error('Simulation:EmptyFigure','No feasible default reference.'); end
tg=logspace(-2,4,220); reached=zeros(nA,numel(tg));
for a=1:nA
    for j=1:numel(tg); reached(a,j)=100*mean(times(a,defaultok)<=tg(j)); end
end
summary=table(string(d.methods(:)),sum(isfinite(G(:,refok)),2), ...
    sum(G(:,refok)<=1,2),mean(times(:,defaultok)<=1,2)*100, ...
    mean(times(:,defaultok)<=10,2)*100, ...
    'VariableNames',{'method','feasible_pairs','within_one_percent', ...
    'default_target_by_1s_pct','default_target_by_10s_pct'});
data=struct('schema_version',s.schema_version,'code_version',s.code_version, ...
    'figure_id',11,'complete',true,'methods',{d.methods}, ...
    'paired_gap_pct',G,'paired_reference_eligible',refok, ...
    'time_to_target_s',times,'default_reference_eligible',defaultok, ...
    'accuracy_grid_pct',epsgrid,'accuracy_success_pct',success, ...
    'time_grid_s',tg,'target_reached_pct',reached,'summary',summary, ...
    'completed_at',char(datetime('now')));
archive=simtwc_archive_path(11); simtwc_save_atomic(archive,data);
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end

function yes=valid(q)
yes=ok(q)&&q.feasible&&isfinite(q.metrics.crlb)&&q.metrics.crlb>0;
end
