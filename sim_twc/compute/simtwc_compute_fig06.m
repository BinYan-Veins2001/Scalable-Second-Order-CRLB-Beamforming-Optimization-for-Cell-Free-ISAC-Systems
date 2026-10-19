function archive = simtwc_compute_fig06()
%SIMTWC_COMPUTE_FIG06 Two-by-two Hessian/negative-curvature factor ablation.

s=simtwc_setup(); spec=s.fig06;
variants=struct( ...
    'key',{'exact_nc','exact_no_nc','gn_nc','gn_no_nc'}, ...
    'name',{'Exact Hessian + NC','Exact Hessian, NC off', ...
    'Gauss-Newton + NC','Gauss-Newton, NC off'}, ...
    'hessian',{'exact_lagrangian','exact_lagrangian', ...
    'gauss_newton_psd','gauss_newton_psd'}, ...
    'negative_curvature',{true,false,true,false});
data=simtwc_initialize_data(6,spec,{variants.key}); path=simtwc_archive_path(6);
nV=numel(variants); nR=numel(spec.regimes); nT=numel(spec.seeds);
if ~isfield(data,'runs') || ~same_size3(data.runs,nV,nR,nT)
    data.runs=cell(nV,nR,nT);
end
for r=1:nR
    regime=spec.regimes(r);
    for t=1:nT
        p=simtwc_make_problem(regime,spec.seeds(t));
        for v=1:nV
            if ok(data.runs{v,r,t}); continue; end
            ropts=struct('rtr_opts',struct( ...
                'hessian_model',variants(v).hessian, ...
                'enable_tcg_negative_curvature',variants(v).negative_curvature, ...
                'enable_ritz_escape',false));
            fprintf('[Fig06] %s, seed=%d, %s\n', ...
                regime.name,spec.seeds(t),variants(v).name);
            data.runs{v,r,t}=simtwc_try_run('costar',p,false,ropts);
            simtwc_save_atomic(path,data);
        end
    end
end
data.variants=variants; data.regimes=spec.regimes; data.seeds=spec.seeds;
data.ablation_scope=['Hessian/tCG factors apply to both optimization stages; ' ...
    'the independent terminal curvature oracle always uses the exact Lagrangian Hessian.'];
[data.final_crlb,data.kkt_residual,data.runtime_s, ...
    data.tcg_nc_detected,data.tcg_nc_accepted]=extract_results(data.runs);
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

function [C,K,T,ND,NA]=extract_results(runs)
sz=size(runs); C=nan(sz); K=nan(sz); T=nan(sz); ND=nan(sz); NA=nan(sz);
for k=1:numel(runs)
    q=runs{k}; if ~ok(q); continue; end
    C(k)=q.metrics.crlb; K(k)=q.kkt.r_kkt; T(k)=q.runtime_median_s;
    if isfield(q,'solution') && isfield(q.solution,'work')
        w=q.solution.work;
        if isfield(w,'tcg_negative_curvature_detections')
            ND(k)=w.tcg_negative_curvature_detections;
        end
        if isfield(w,'tcg_negative_curvature_accepted')
            NA(k)=w.tcg_negative_curvature_accepted;
        end
    end
end
end
