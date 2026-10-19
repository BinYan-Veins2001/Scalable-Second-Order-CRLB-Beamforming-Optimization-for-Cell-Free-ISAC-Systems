function [sol,hist]=costar_paper_rtr(cfg,channels,smat,initial_wcell,o)
%COSTAR_PAPER_RTR Feasible safeguarded RTR with a cheap first-order accelerator.
% Implements the coordinated NNLS safeguard, first-order acceptance, and
% deterministic SPSC in the current paper. Curvature trials retain the
% stricter restored-model accuracy test.
% All finite resource or precision exits are explicitly uncertified.
N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U; d=2*N*U;
W=vectorize_wcell_to_matrix(cfg,initial_wcell); z=complex_matrix_to_real_vector(W);
Delta=o.delta0; nu=o.nu0; level=0; kappa_n=o.normal_second_order_cap;
Lc=2/cfg.power.Pmax;
if U>1
    Lc=max(Lc,max(2*cfg.comm.Gamma*sum(abs(channels.h_stack).^2,1)/cfg.comm.noise_power));
end
bc=2*Lc+1;
names={'derivative_calls','value_calls','constraint_calls','hvp_calls','exact_hvp_calls', ...
    'fd_hvp_calls','meo_calls','meo_hvp_calls','nnls_calls','inward_qp_calls', ...
    'inward_qp_sweeps','threshold_refinements','safeguard_calls','safeguard_trials', ...
    'safeguard_accepted','tcg_calls','tcg_iterations','tcg_negative_curvature_detections', ...
    'tcg_negative_curvature_accepted','structured_checks','structured_passes', ...
    'structured_tangent_passes','tangent_checks','negative_curvature_steps', ...
    'normal_newton_calls','normal_newton_accepted','restoration_iterations', ...
    'restoration_failures','accepted_steps','rejected_steps','stage_refinements', ...
    'precision_limited_trials','rtr_iterations','guard_cache_builds', ...
    'guard_cache_trials','guard_direct_fallbacks','feasible_accelerator_accepted'};
work=struct(); for j=1:numel(names); work.(names{j})=0; end
timing=struct('geometry_s',0,'safeguard_s',0,'tcg_s',0,'restoration_s',0,'spsc_s',0);
hist=struct('iter',[],'crlb',[],'objective',[],'time_s',[],'kkt',[], ...
    'stationarity',[],'primal_residual',[],'complementarity',[], ...
    'tangent_grad_norm',[],'active_count',[],'delta',[],'accepted',[], ...
    'ratio',[],'step_type',{{}},'model_error',[],'restoration_defect',[],'normal_ratio',[], ...
    'restoration_iterations',[],'lambda_min',[],'curvature_lower_bound',[], ...
    'critical_cone_tolerance',[],'meo_status',{{}},'step_norm',[],'cg_iter',[], ...
    'tcg_final_res_norm',[],'level',[],'nu',[],'tau',[],'safeguard_decrease',[], ...
    'safeguard_alpha',[],'active_membership',[],'active_change_count',[], ...
    'z_snapshots',[],'snapshot_iter',[]);
success=false; message='iteration budget exhausted; uncertified';
oracle=struct('status','not_called','ambient',false,'eigenvalue_lower_bound',-inf);
cert_z=[]; cert_geometry_nu=NaN; cone_tolerance=inf; lambda_min=NaN; timer=tic;
for iter=0:o.max_iter
    [f,g,c,G,out]=evaluate(z);
    if max([0;c])>o.feas_tol || ~isfinite(f)
        message='nonfeasible or nonidentifiable input/iterate; uncertified'; break;
    end
    t=tic; geo=costar_paper_geometry(g,c,G,U,nu,o); timing.geometry_s=timing.geometry_s+toc(t);
    work=add_counts(work,geo.work);
    % A full-space certificate may use any nonnegative multiplier satisfying
    % the same complete KKT tolerance. No tangent certificate is transferred
    % from the larger numerical working set used by this optional probe.
    if o.ambient_endpoint_probe && geo.success && geo.kkt>o.epsilon_KKT && ...
            geo.kkt<=o.endpoint_probe_factor*o.epsilon_KKT
        probe_nu=max(nu,2*o.tau_bar/o.kappa_tau);
        t=tic; probe=costar_paper_geometry(g,c,G,U,probe_nu,o);
        timing.geometry_s=timing.geometry_s+toc(t); work=add_counts(work,probe.work);
        if probe.success && probe.kkt<=o.epsilon_KKT
            t=tic; pi=costar_sosp_structured_psd(cfg,channels,smat, ...
                real_vector_to_complex_matrix(z,N,U),probe.lambda,o.epsilon_H,zeros(d,0));
            timing.spsc_s=timing.spsc_s+toc(t);
            work.structured_checks=work.structured_checks+1;
            work.hvp_calls=work.hvp_calls+pi.hvp_calls; work.exact_hvp_calls=work.exact_hvp_calls+pi.hvp_calls;
            work.structured_passes=work.structured_passes+(pi.conditional_psd&&pi.ambient);
            if pi.conditional_psd && pi.ambient
                geo=probe; oracle=pi; cert_z=z; cert_geometry_nu=probe_nu;
                lambda_min=pi.rayleigh; cone_tolerance=max(0,-pi.eigenvalue_lower_bound);
                success=true; message='full-space numerical KKT/curvature endpoint certificate passed';
                record(iter,f,c,out,geo); hist.step_type{end}='ambient_endpoint_certificate';
                hist.lambda_min(end)=lambda_min; hist.curvature_lower_bound(end)=pi.eigenvalue_lower_bound;
                hist.critical_cone_tolerance(end)=cone_tolerance; hist.meo_status{end}=pi.status;
                break;
            end
        end
    end
    record(iter,f,c,out,geo);
    if ~geo.success; message=['working-set coordination failed: ' geo.message]; break; end
    if o.verbose && (iter==0 || mod(iter,o.print_every)==0)
        fprintf('COSTAR paper %d f=%.9g r=%.3g nu=%.3g tau=%.3g Delta=%.3g\n', ...
            iter,f,geo.kkt,nu,geo.tau,Delta);
    end
    if iter==o.max_iter; break; end
    project=@(v) v-geo.Q*(geo.Q'*v);
    Hv=@(v) projected_hvp(z,v,geo.lambda,geo.Q);
    if geo.kkt>nu
        t=tic; [zf,DF,alpha,ok]=safeguard(z,geo); timing.safeguard_s=timing.safeguard_s+toc(t);
        hist.safeguard_decrease(end)=DF; hist.safeguard_alpha(end)=alpha;
        if ~ok; message='safeguard precision/backtracking budget reached; uncertified'; break; end
        oldDelta=Delta; t=tic; work.tcg_calls=work.tcg_calls+1;
        cgopts=struct('max_iter',o.cg_max_iter,'forcing_power',1, ...
            'forcing_constant',o.cg_relative_tol,'absolute_floor',0, ...
            'nc_tol',0,'flat_tol',0,'detect_negative_curvature',true);
        [s,cg]=trust_region_tcg(geo.gt,Hv,Delta,cgopts); s=project(s);
        timing.tcg_s=timing.tcg_s+toc(t); work.tcg_iterations=work.tcg_iterations+cg.iter;
        work.tcg_negative_curvature_detections=work.tcg_negative_curvature_detections+cg.negative_curvature;
        hist.cg_iter(end)=cg.iter; hist.tcg_final_res_norm(end)=cg.final_res_norm;
        Hs=Hv(s); pred=-(g'*s+0.5*s'*Hs);
        t=tic;
        if o.fast_retraction
            [zt,~,rest]=accelerated_restore(z,s,geo,c,pred,kappa_n,Delta,DF);
        else
            rest=struct('model_ok',false,'accelerator_ok',false);
        end
        if ~(rest.model_ok || rest.accelerator_ok)
            [zt,~,rest]=ordinary_restore(z,s,geo,c,pred,kappa_n,Delta,DF);
        end
        timing.restoration_s=timing.restoration_s+toc(t);
        [actual,noise]=decrease(z,zt); ratio=actual/pred;
        quality_ok=rest.model_ok || (o.feasible_accelerator && rest.accelerator_ok);
        accepted=pred>0 && isfinite(pred) && quality_ok && actual>noise && ratio>=o.eta_accept;
        if ~accepted && o.first_order_normal_restore && geo.full_column_rank && norm(s)>0 && pred>0
            t=tic;
            [zn,~,rn,trialcap]=normal_restore(z,s,geo,c,pred,kappa_n,norm(s));
            [an,en]=decrease(z,zn); ratn=an/pred;
            if rn.model_ok && norm(zn-z-s)<=kappa_n*norm(s)^2 && an>en && ratn>=o.eta_accept
                zt=zn; rest=rn; actual=an; noise=en; ratio=ratn; accepted=true;
                work.normal_newton_accepted=work.normal_newton_accepted+1;
            end
            timing.restoration_s=timing.restoration_s+toc(t);
        end
        hist.ratio(end)=ratio; hist.model_error(end)=abs(actual-pred);
        hist.normal_ratio(end)=norm(zt-z-s)/max(norm(s)^2,realmin);
        hist.restoration_defect(end)=rest.defect; hist.restoration_iterations(end)=rest.iterations;
        if accepted && actual>=o.omega_F*DF
            z=zt; hist.step_type{end}='tangent_tcg'; hist.step_norm(end)=norm(s);
            if ~rest.model_ok
                hist.step_type{end}='feasible_safeguarded_tcg';
                work.feasible_accelerator_accepted=work.feasible_accelerator_accepted+1;
            end
            Delta=update_radius(oldDelta,ratio,norm(s));
            work.tcg_negative_curvature_accepted=work.tcg_negative_curvature_accepted+cg.negative_curvature;
        else
            z=zf; Delta=oldDelta; hist.step_type{end}='nnls_safeguard';
            if o.shrink_failed_tangent && ~accepted
                Delta=o.gamma_shrink*oldDelta;
            end
            hist.step_norm(end)=norm(z-hist.z_snapshots(:,end));
            work.safeguard_accepted=work.safeguard_accepted+1;
            work.rejected_steps=work.rejected_steps+1;
            work.restoration_failures=work.restoration_failures+~rest.model_ok;
        end
        hist.accepted(end)=true; work.accepted_steps=work.accepted_steps+1;
    else
        t=tic; work.structured_checks=work.structured_checks+1;
        oracle=costar_sosp_structured_psd(cfg,channels,smat, ...
            real_vector_to_complex_matrix(z,N,U),geo.lambda,o.epsilon_H,geo.Q);
        timing.spsc_s=timing.spsc_s+toc(t);
        work.hvp_calls=work.hvp_calls+oracle.hvp_calls;
        work.exact_hvp_calls=work.exact_hvp_calls+oracle.hvp_calls;
        work.structured_passes=work.structured_passes+(oracle.conditional_psd&&oracle.ambient);
        work.structured_tangent_passes=work.structured_tangent_passes+(oracle.conditional_psd&&~oracle.ambient);
        work.tangent_checks=work.tangent_checks+(~oracle.ambient&&~strcmp(oracle.mode,'empty_tangent'));
        cert_z=z; lambda_min=oracle.rayleigh; cone_tolerance=calibrate(oracle,geo);
        hist.lambda_min(end)=lambda_min; hist.curvature_lower_bound(end)=oracle.eigenvalue_lower_bound;
        hist.critical_cone_tolerance(end)=cone_tolerance; hist.meo_status{end}=oracle.status;
        if oracle.conditional_psd && geo.kkt<=o.epsilon_KKT && cone_tolerance<=o.epsilon_H
            success=true; message='deterministic numerical first/second-order stopping tests passed';
            hist.step_type{end}='certified_endpoint'; break;
        end
        if strcmp(oracle.status,'negative_curvature') && geo.full_column_rank && geo.positive_multipliers
            v=oracle.negative_direction; if geo.gt'*v>0; v=-v; end
            Delta=o.delta0; ok=false;
            for bt=1:o.nc_max_backtracks
                s=Delta*v; Hs=Hv(s); pred=-(g'*s+0.5*s'*Hs);
                t=tic; [zt,~,rest]=ordinary_restore(z,s,geo,c,pred,kappa_n,Delta);
                trial_cap=kappa_n; used_normal=false;
                [actual,noise]=decrease(z,zt); ratio=actual/pred;
                accepted=rest.model_ok && pred>0 && actual>noise && ratio>=o.eta_accept;
                if ~accepted
                    [zt,~,rest,trial_cap]=normal_restore(z,s,geo,c,pred,kappa_n,Delta);
                    used_normal=true; [actual,noise]=decrease(z,zt); ratio=actual/pred;
                    accepted=rest.model_ok && pred>0 && actual>noise && ratio>=o.eta_accept;
                end
                timing.restoration_s=timing.restoration_s+toc(t);
                if accepted
                    z=zt; kappa_n=trial_cap; ok=true;
                    work.normal_newton_accepted=work.normal_newton_accepted+used_normal;
                    work.negative_curvature_steps=work.negative_curvature_steps+1;
                    work.accepted_steps=work.accepted_steps+1;
                    hist.accepted(end)=true; hist.step_type{end}='spsc_negative_curvature';
                    hist.step_norm(end)=norm(s); hist.ratio(end)=ratio;
                    hist.restoration_defect(end)=rest.defect; hist.restoration_iterations(end)=rest.iterations;
                    hist.model_error(end)=abs(actual-pred);
                    Delta=update_radius(Delta,ratio,norm(s)); break;
                end
                work.rejected_steps=work.rejected_steps+1; Delta=o.gamma_shrink*Delta;
                if Delta==0 || (o.minimum_radius>0 && Delta<o.minimum_radius); break; end
            end
            if ~ok; message='negative-curvature numerical search exhausted; uncertified'; break; end
        elseif ~oracle.conditional_psd && ~strcmp(oracle.status,'negative_curvature')
            message=['structured factorization unresolved: ' oracle.status]; break;
        else
            hist.step_type{end}='accuracy_refinement';
        end
        level=level+1; work.stage_refinements=work.stage_refinements+1;
        nu=nu*o.refinement_factor;
        if level>=o.max_refinements || nu==0
            message='accuracy refinement budget reached; uncertified'; break;
        end
    end
    work.rtr_iterations=work.rtr_iterations+1;
end

% A final-point recomputation prevents stale pre-step residuals/certificates.
[f,g,c,G,out]=evaluate(z);
final_nu=nu;
if success && isequal(cert_z,z) && isfinite(cert_geometry_nu); final_nu=cert_geometry_nu; end
t=tic; final=costar_paper_geometry(g,c,G,U,final_nu,o); timing.geometry_s=timing.geometry_s+toc(t);
work=add_counts(work,final.work);
if isempty(hist.iter) || ~isequal(z,hist.z_snapshots(:,end)); record(iter+1,f,c,out,final); end
endpoint_bound=-inf;
if ~isequal(cert_z,z) || ~final.success || final.kkt>o.epsilon_KKT || ...
        ~isfield(oracle,'conditional_psd') || ~oracle.conditional_psd
    success=false; cone_tolerance=inf;
    if ~isequal(cert_z,z); lambda_min=NaN; end
elseif success
    cone_tolerance=calibrate(oracle,final); success=cone_tolerance<=o.epsilon_H;
    endpoint_bound=oracle.eigenvalue_lower_bound;
end
W=real_vector_to_complex_matrix(z,N,U); wcell=matrix_to_wcell(cfg,W);
sol=struct('z',z,'Wmat',W,'wcell',{wcell},'runtime',toc(timer), ...
    'metrics',evaluate_beamforming_solution(cfg,channels,smat,wcell), ...
    'success',success,'message',message,'conditional_sosp_reached',success, ...
    'second_order_verified',false,'numerical_joint_pass',success, ...
    'lambda',final.lambda,'active',final.active, ...
    'working_thresholds',[final.sinr_threshold final.power_threshold], ...
    'kkt_residual',final.kkt,'stationarity',final.rg,'complementarity',final.r_lambda, ...
    'primal_residual',max([0;c]),'lambda_min',lambda_min, ...
    'curvature_lower_bound',endpoint_bound,'critical_cone_tolerance',cone_tolerance, ...
    'delta',Delta,'work',work,'timing',timing,'options',o,'level',level,'nu',nu, ...
    'endpoint_geometry_nu',final_nu, ...
    'constraint_curvature_bound',Lc,'safeguard_coefficient',bc,'normal_second_order_cap',kappa_n);
sol.certificate=struct('oracle',oracle,'oracle_is_at_final_point',isequal(cert_z,z), ...
    'deterministic_numerical',success,'floating_point_validated',false, ...
    'conditional_high_probability',false,'epsilon_KKT',o.epsilon_KKT, ...
    'epsilon_H',o.epsilon_H,'feasibility_tolerance',o.feas_tol, ...
    'critical_cone_tolerance',cone_tolerance,'geometry',final, ...
    'scope','ambient or calibrated numerical critical cone; see oracle mode', ...
    'failure_probability',0);
hist.work=work; hist.options=o; hist.phase=repmat("paper_aligned",size(hist.iter));
hist.terminal_start_iter=1;

    function [fv,gv,cv,Gv,ov]=evaluate(zv)
        work.derivative_calls=work.derivative_calls+1;
        Wv=real_vector_to_complex_matrix(zv,N,U);
        [fv,gc,ov]=compute_crlb_objective_grad(cfg,smat,Wv);
        if ~(ov.det>cfg.alm.det_floor && ov.sum_ab>0); fv=inf; end
        gv=complex_matrix_to_real_vector(gc); [cv,Gv]=costar_sosp_constraints(cfg,channels,Wv);
        work.constraint_calls=work.constraint_calls+1;
    end
    function [fv,cv]=value(zv)
        work.value_calls=work.value_calls+1; work.constraint_calls=work.constraint_calls+1;
        Wv=real_vector_to_complex_matrix(zv,N,U);
        fv=costar_paper_value(cfg,smat,Wv); cv=costar_sosp_constraints(cfg,channels,Wv);
    end
    function hv=projected_hvp(zv,v,lambda,Q)
        work.hvp_calls=work.hvp_calls+1; work.exact_hvp_calls=work.exact_hvp_calls+1;
        v=v-Q*(Q'*v);
        hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec(cfg,channels,smat, ...
            real_vector_to_complex_matrix(zv,N,U),real_vector_to_complex_matrix(v,N,U),lambda));
        hv=hv-Q*(Q'*hv);
    end
    function [actual,noise]=decrease(base,trial)
        try
            [actual,noise]=costar_sosp_stable_decrease(cfg,smat, ...
                real_vector_to_complex_matrix(base,N,U), ...
                real_vector_to_complex_matrix(trial-base,N,U));
        catch
            actual=-inf; noise=inf;
        end
    end
    function [zf,DF,alpha,ok]=safeguard(base,geo)
        work.safeguard_calls=work.safeguard_calls+1;
        alpha=o.alpha0; ok=false; DF=-inf; zf=base; zeta=geo.zeta;
        if o.guard_directional
            zn=norm(zeta); p=-zeta/zn;
            Pw=real_vector_to_complex_matrix(p,N,U);
            YY=channels.h_stack'*Pw;
            qdir=zeros(U+cfg.sys.M,1);
            for uu=1:U
                others=[1:uu-1 uu+1:U];
                qdir(uu)=(cfg.comm.Gamma*sum(abs(YY(uu,others)).^2)-abs(YY(uu,uu))^2)/cfg.comm.noise_power;
            end
            for ap=1:cfg.sys.M
                rows=(ap-1)*cfg.sys.Nt+(1:cfg.sys.Nt);
                qdir(U+ap)=sum(abs(Pw(rows,:)).^2,'all')/cfg.power.Pmax;
            end
            inward_slopes=-G(:,geo.active)'*geo.normal_inward;
            beta=max([0;(qdir(geo.active)+o.guard_direction_margin)./inward_slopes]);
            v=beta*geo.normal_inward; alpha=o.guard_step0;
            if o.guard_scaled_inward
                [vs,vi]=costar_guard_scaled_inward(G(:,geo.active),qdir(geo.active),o.guard_direction_margin,o);
                if vi.success; v=vs; end
            end
        else
            p=-zeta; v=bc*(zeta'*zeta)*geo.normal_inward;
            if o.guard_initial_scaling
                correction_scale=bc*norm(zeta)*norm(geo.normal_inward);
                if correction_scale>0
                    alpha=min(alpha,o.guard_correction_fraction/correction_scale);
                end
            end
        end
        initial_alpha=alpha; cache=[];
        % If a cached screen is inconclusive, repeat the same search with
        % direct evaluations. Polynomial rounding cannot remove the backup.
        for attempt=1:2
            alpha=initial_alpha;
            for bt=1:o.guard_max_backtracks
                work.safeguard_trials=work.safeguard_trials+1;
                trial=base+alpha*p+alpha^2*v;
                if o.guard_directional; target=o.eta_F*alpha*zn;
                else; target=o.eta_F*alpha*(zeta'*zeta); end
                cached=o.guard_cached && attempt==1 && bt>=o.guard_cache_min_trials;
                possible=true;
                if cached
                    if isempty(cache)
                        cache=costar_guard_path_cache(cfg,channels,smat, ...
                            real_vector_to_complex_matrix(base,N,U), ...
                            real_vector_to_complex_matrix(p,N,U),real_vector_to_complex_matrix(v,N,U));
                        work.guard_cache_builds=work.guard_cache_builds+1;
                    end
                    [fv,cv,DF,noise,cscreen]=costar_guard_path_evaluate(cache,alpha);
                    work.guard_cache_trials=work.guard_cache_trials+1;
                    possible=max([0;cv-cscreen.constraint_noise])<=o.feas_tol && ...
                        isfinite(fv) && DF+noise>=target;
                end
                if possible
                    Wtrial=real_vector_to_complex_matrix(trial,N,U);
                    cv=costar_sosp_constraints(cfg,channels,Wtrial);
                    work.constraint_calls=work.constraint_calls+1;
                    if max([0;cv])<=o.feas_tol
                        fv=costar_paper_value(cfg,smat,Wtrial); work.value_calls=work.value_calls+1;
                        if isfinite(fv)
                            [DF,noise]=decrease(base,trial);
                            if DF>=target && DF>noise
                                zf=trial; ok=true; return;
                            end
                        end
                    end
                end
                alpha=alpha/2;
            end
            if isempty(cache); break; end
            work.guard_direct_fallbacks=work.guard_direct_fallbacks+1;
        end
    end
    function [trial,cv,result]=ordinary_restore(base,s,geo,cb,pred,cap,trust,DF)
        if nargin<8; DF=inf; end
        trial=base+s; result=struct('model_ok',false,'accelerator_ok',false,'defect',inf,'iterations',0);
        [fv,cv]=value(trial);
        for r=0:o.restoration_max_iter
            normal=trial-base-s;
            result.defect=abs(geo.zeta'*normal)+sum(geo.lambda(geo.active).*abs(cv(geo.active)-cb(geo.active)));
            result.model_ok=pred>0 && isfinite(fv) && max([0;cv])<=o.feas_tol && ...
                norm(normal)<=cap*norm(s)^2 && result.defect<=o.restoration_fraction*pred;
            if o.feasible_accelerator && isfinite(DF) && isfinite(fv) && max([0;cv])<=o.feas_tol
                [da,dn]=decrease(base,trial);
                result.accelerator_ok=da>dn && da>=o.omega_F*DF && pred>0 && da/pred>=o.eta_accept && ...
                    (~o.accelerator_norm_cap || norm(normal)<=cap*norm(s)^2);
            end
            if result.model_ok || result.accelerator_ok || r==o.restoration_max_iter; return; end
            [cv,Gv]=costar_sosp_constraints(cfg,channels,real_vector_to_complex_matrix(trial,N,U));
            work.constraint_calls=work.constraint_calls+1;
            ids=find(cv(1:U)>-geo.sinr_threshold);
            pids=U+find(cv(U+1:end)>-geo.power_threshold);
            if isempty(pids); Qp=zeros(d,0); else; Qp=orth(Gv(:,pids),o.rank_tol); end
            if ~isempty(ids)
                A=Gv(:,ids); A=A-Qp*(Qp'*A);
                K=A'*A; mu=o.restoration_regularization*max(1,norm(K,2));
                correction=-A*((K+mu*eye(numel(ids)))\cv(ids));
            else; correction=zeros(d,1); end
            % Largest xi in [0,1] preserving the total trust-region bound.
            a=correction'*correction; b=2*((trial-base)'*correction);
            cc=norm(trial-base)^2-trust^2; xi=1;
            if o.restoration_total_radius && norm(trial+correction-base)>trust && a>0
                xi=max(0,min(1,(-b+sqrt(max(0,b*b-4*a*cc)))/(2*a)));
            end
            trial=trial+xi*correction; Wv=real_vector_to_complex_matrix(trial,N,U);
            for ap=1:cfg.sys.M
                rows=(ap-1)*cfg.sys.Nt+(1:cfg.sys.Nt); power=sum(abs(Wv(rows,:)).^2,'all');
                if power>cfg.power.Pmax
                    Wv(rows,:)=Wv(rows,:)*sqrt(cfg.power.Pmax/power)*(1-4*eps);
                end
            end
            trial=complex_matrix_to_real_vector(Wv); [fv,cv]=value(trial);
            result.iterations=result.iterations+1; work.restoration_iterations=work.restoration_iterations+1;
        end
    end
    function [trial,cv,result]=accelerated_restore(base,s,geo,cb,pred,cap,trust,DF)
        trial=retract_power(base+s,cb,geo);
        inward_margin=min(o.sinr_inward_margin,o.sinr_margin_fraction*geo.sinr_threshold);
        result=struct('model_ok',false,'accelerator_ok',false,'defect',inf,'iterations',0);
        for jj=0:o.fast_retraction_iterations
            [fv,cv]=value(trial); normal=trial-base-s;
            result.defect=abs(geo.zeta'*normal)+sum(geo.lambda(geo.active).*abs(cv(geo.active)-cb(geo.active)));
            feasible=isfinite(fv) && max([0;cv])<=o.feas_tol;
            norm_ok=norm(normal)<=cap*norm(s)^2;
            result.model_ok=pred>0 && feasible && norm_ok && result.defect<=o.restoration_fraction*pred;
            if feasible
                [da,dn]=decrease(base,trial);
                result.accelerator_ok=o.feasible_accelerator && da>dn && da>=o.omega_F*DF && ...
                    pred>0 && da/pred>=o.eta_accept && (~o.accelerator_norm_cap || norm_ok);
            end
            if result.model_ok || result.accelerator_ok || jj==o.fast_retraction_iterations; return; end
            [cv,Gv]=costar_sosp_constraints(cfg,channels,real_vector_to_complex_matrix(trial,N,U));
            work.constraint_calls=work.constraint_calls+1;
            ids=find(cv(1:U)>-inward_margin);
            if isempty(ids); return; end
            pids=U+find(cb(U+1:end)>-geo.power_threshold);
            if isempty(pids); Qp=zeros(d,0); else; Qp=orth(Gv(:,pids),o.rank_tol); end
            A=Gv(:,ids); A=A-Qp*(Qp'*A);
            K=A'*A;
            correction=-A*((K+o.restoration_regularization*eye(numel(ids)))\(cv(ids)+inward_margin));
            if norm(correction)>trust; correction=correction*(trust/norm(correction)); end
            trial=retract_power(trial+correction,cb,geo);
            result.iterations=result.iterations+1; work.restoration_iterations=work.restoration_iterations+1;
        end
    end
    function trial=retract_power(trial,cb,geo)
        Wv=real_vector_to_complex_matrix(trial,N,U);
        for ap=1:cfg.sys.M
            rows=(ap-1)*cfg.sys.Nt+(1:cfg.sys.Nt); power=sum(abs(Wv(rows,:)).^2,'all');
            if power>0 && (cb(U+ap)>-geo.power_threshold || power>cfg.power.Pmax)
                Wv(rows,:)=Wv(rows,:)*sqrt(cfg.power.Pmax/power)*(1-4*eps);
            end
        end
        trial=complex_matrix_to_real_vector(Wv);
    end
    function [trial,cv,result,trialcap]=normal_restore(base,s,geo,cb,pred,cap,trust)
        work.normal_newton_calls=work.normal_newton_calls+1; trialcap=cap;
        result=struct('model_ok',false,'defect',inf,'iterations',1); trial=base+s;
        [cs,Gs]=costar_sosp_constraints(cfg,channels,real_vector_to_complex_matrix(trial,N,U));
        work.constraint_calls=work.constraint_calls+1;
        if isempty(geo.Q); normal=zeros(d,1);
        else
            J=Gs(:,geo.active)'*geo.Q;
            if rcond(J)<eps; cv=cs; return; end
            rhs=-(cs(geo.active)-cb(geo.active)+o.normal_inward_margin*trust^3);
            normal=geo.Q*(J\rhs); trial=trial+normal;
        end
        [fv,cv]=value(trial);
        residual=cv(geo.active)-cb(geo.active)+o.normal_inward_margin*trust^3;
        nonlinear_ok=max([0;abs(residual)])<=o.normal_inward_margin*trust^3/2;
        trialcap=max(cap,2*norm(normal)/trust^2);
        result.defect=abs(geo.zeta'*normal)+sum(geo.lambda(geo.active).*abs(cv(geo.active)-cb(geo.active)));
        result.model_ok=nonlinear_ok && isfinite(fv) && max([0;cv])<=o.feas_tol && ...
            result.defect<=o.restoration_fraction*pred && norm(normal)<=trialcap*norm(s)^2;
    end
    function val=calibrate(info,geo)
        if ~info.conditional_psd; val=inf;
        elseif info.ambient || ~any(geo.active); val=max(0,-info.eigenvalue_lower_bound);
        elseif geo.full_column_rank && geo.positive_multipliers
            chi=norm(geo.zeta)/(min(geo.lambda(geo.active))*geo.sigma_min);
            val=max(0,-info.eigenvalue_lower_bound)+info.norm_bound*chi*(2+chi);
        else; val=inf; end
    end
    function radius=update_radius(radius,ratio,stepnorm)
        if ratio<o.eta_shrink; radius=o.gamma_shrink*radius;
        elseif ratio>=o.eta_expand && stepnorm>=o.eta_boundary*radius
            radius=min(o.gamma_expand*radius,o.delta_max);
        end
    end
    function record(k,fv,cv,ov,geo)
        hist.iter(end+1)=k; hist.crlb(end+1)=ov.crlb; hist.objective(end+1)=fv;
        hist.time_s(end+1)=toc(timer); hist.kkt(end+1)=geo.kkt;
        hist.stationarity(end+1)=geo.rg; hist.complementarity(end+1)=geo.r_lambda;
        hist.primal_residual(end+1)=max([0;cv]); hist.tangent_grad_norm(end+1)=norm(geo.gt);
        hist.active_count(end+1)=nnz(geo.active); hist.delta(end+1)=Delta;
        hist.accepted(end+1)=false; hist.ratio(end+1)=NaN; hist.step_type{end+1}='diagnostic';
        hist.model_error(end+1)=NaN; hist.restoration_defect(end+1)=NaN;
        hist.normal_ratio(end+1)=NaN;
        hist.restoration_iterations(end+1)=0; hist.lambda_min(end+1)=NaN;
        hist.curvature_lower_bound(end+1)=-inf; hist.critical_cone_tolerance(end+1)=inf;
        hist.meo_status{end+1}='not_called'; hist.step_norm(end+1)=0;
        hist.cg_iter(end+1)=0; hist.tcg_final_res_norm(end+1)=NaN;
        hist.level(end+1)=level; hist.nu(end+1)=nu; hist.tau(end+1)=geo.tau;
        hist.safeguard_decrease(end+1)=NaN; hist.safeguard_alpha(end+1)=NaN;
        hist.active_membership(:,end+1)=geo.active;
        if size(hist.active_membership,2)>1
            hist.active_change_count(end+1)=nnz(xor(geo.active,hist.active_membership(:,end-1)));
        else; hist.active_change_count(end+1)=0; end
        hist.z_snapshots(:,end+1)=z; hist.snapshot_iter(end+1)=k;
    end
end

function a=add_counts(a,b)
names=fieldnames(b);
for k=1:numel(names)
    if ~isfield(a,names{k}); a.(names{k})=0; end
    a.(names{k})=a.(names{k})+b.(names{k});
end
end
