function h = simtwc_canonical_history( ...
    key,p,sol,hist,Wfinal,elapsed,entry_w,entry_time)
%SIMTWC_CANONICAL_HISTORY Build time/CRLB/KKT traces from stored iterates.

N=p.cfg.sys.M*p.cfg.sys.Nt;
U=p.cfg.sys.U;
s=simtwc_config();
h=struct('time_s',[],'crlb',[],'kkt_residual',[], ...
    'feasible',[],'phase',strings(0,1),'z_snapshots',[]);
zfinal=complex_matrix_to_real_vector(Wfinal);

switch key
    case {'zf_bf','so_bf'}
        append(elapsed,zfinal,"static");
    case 'pfo_alm'
        append_entry();
        append_snapshot_history(hist,entry_time,"PFO-ALM");
    case 'ptr_alm'
        append_entry();
        if isfield(hist,'stage1_pfo')
            append_snapshot_history(hist.stage1_pfo,entry_time,"PFO globalization");
        end
        if isfield(hist,'stage2_ptr')
            append_snapshot_history(hist.stage2_ptr, ...
                entry_time+field_or(hist,'ptr_start_time_s',0),"PTR refinement");
        end
        h.ptr_start_time_s=entry_time+field_or(hist,'ptr_start_time_s',0);
    case 'mrfp_sca'
        append_entry();
        append_snapshot_history(hist,entry_time,"MRFP-SCA");
    case 'costar'
        append_entry();
        if isfield(hist,'rtr')
            append_snapshot_history(hist.rtr,entry_time,"COSTAR");
        end
        h.rtr_start_time_s=entry_time;
        h.terminal_start_time_s=entry_time+hist.terminal_start_time_s;
        h.raw=hist;
    case 'sdr_sca'
        append_entry();
        append_snapshot_history(hist,entry_time,"SDR-SCA");
end

if isempty(h.z_snapshots) || ...
        norm(h.z_snapshots(:,end)-zfinal)>1e-12*max(1,norm(zfinal))
    append(elapsed,zfinal,"final");
elseif h.time_s(end)<elapsed
    append(elapsed,zfinal,"termination");
end
[h.time_s,order]=sort(real(h.time_s(:)),'ascend');
h.z_snapshots=h.z_snapshots(:,order);
h.phase=h.phase(order);

n=size(h.z_snapshots,2);
h.crlb=nan(n,1); h.kkt_residual=nan(n,1); h.feasible=false(n,1);
h.active_sinr_count=zeros(n,1); h.active_power_count=zeros(n,1);
for k=1:n
    W=real_vector_to_complex_matrix(h.z_snapshots(:,k),N,U);
    metrics=evaluate_beamforming_solution( ...
        p.cfg,p.channels,p.smat,matrix_to_wcell(p.cfg,W));
    d=simtwc_unified_kkt(p,W);
    h.crlb(k)=metrics.crlb;
    h.kkt_residual(k)=d.r_kkt;
    h.feasible(k)=d.r_c<=s.feasibility_tol;
    h.active_sinr_count(k)=d.active_sinr_count;
    h.active_power_count(k)=d.active_power_count;
end

    function append_entry()
        if ~isempty(entry_w)
            append(entry_time,complex_matrix_to_real_vector( ...
                vectorize_wcell_to_matrix(p.cfg,entry_w)),"global-ZF entry");
        end
    end

    function append_snapshot_history(q,offset,label)
        if ~isfield(q,'z_snapshots') || isempty(q.z_snapshots); return; end
        zs=q.z_snapshots;
        if isfield(q,'snapshot_iter') && isfield(q,'iter') && isfield(q,'time_s')
            tq=map_times(q.snapshot_iter,q.iter,q.time_s);
        elseif isfield(q,'cumulative_time_s')
            tq=q.cumulative_time_s(:);
        elseif isfield(q,'time_s')
            tq=q.time_s(:);
        else
            error('Simulation:MissingTimestamps', ...
                'Stored iterates of %s have no wall-clock timestamps.',key);
        end
        count=min(size(zs,2),numel(tq));
        for j=1:count; append(offset+tq(j),zs(:,j),label); end
    end

    function append(t,z,label)
        h.time_s(end+1,1)=max(real(t),0);
        h.z_snapshots(:,end+1)=real(z(:));
        h.phase(end+1,1)=string(label);
    end
end

function times=map_times(snapshot_iter,iter,time_s)
snapshot_iter=snapshot_iter(:); iter=iter(:); time_s=time_s(:);
times=nan(size(snapshot_iter));
for j=1:numel(snapshot_iter)
    idx=find(iter==snapshot_iter(j),1,'last');
    if isempty(idx); idx=min(j,numel(time_s)); end
    times(j)=time_s(idx);
end
end

function value=field_or(s,name,fallback)
if isfield(s,name); value=s.(name); else; value=fallback; end
end
