function [sol,hist] = polish_costar_sosp(cfg,channels,smat,initial_wcell,options)
%POLISH_COSTAR_SOSP Accuracy-controlled, matrix-free terminal RTR.
% A conditional active-tangent SOSP flag is distinct from a validated
% floating-point certificate and from unconditional critical-cone optimality.
o=defaults(options); N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U;
W=vectorize_wcell_to_matrix(cfg,initial_wcell);
z=complex_matrix_to_real_vector(W); dim=numel(z); Delta=o.delta0;
work=struct('derivative_calls',0,'hvp_calls',0,'meo_calls',0, ...
    'exact_hvp_calls',0,'fd_hvp_calls',0,'gn_model_builds',0,'tcg_iterations',0, ...
    'meo_hvp_calls',0,'norm_bound_builds',0,'structured_checks',0, ...
    'structured_passes',0,'restoration_iterations',0, ...
    'restoration_failures',0,'accepted_steps',0,'rejected_steps',0, ...
    'tcg_calls',0,'tcg_negative_curvature_detections',0, ...
    'tcg_negative_curvature_accepted',0,'ritz_negative_curvature_detections',0, ...
    'ritz_escape_accepted',0,'precision_limited_trials',0);
hist=struct('iter',[],'crlb',[],'objective',[],'time_s',[], ...
    'kkt',[],'tangent_grad_norm',[],'active_count',[],'delta',[], ...
    'accepted',[],'ratio',[],'step_type',{{}},'model_error',[], ...
    'restoration_defect',[],'restoration_iterations',[], ...
    'lambda_min',[],'curvature_lower_bound',[],'meo_status',{{}}, ...
    'snapshot_iter',[],'z_snapshots',[]);
hist.active_membership=[]; hist.active_change_count=[];
hist.step_norm=[]; hist.cg_iter=[]; hist.tcg_final_res_norm=[];
sol.success=false; sol.message='terminal iteration budget exhausted';
sol.conditional_sosp_reached=false; info=struct('status','not_called');
lambda_min=NaN; meo_z=[]; geometry=struct(); timer=tic;
spent_failure_probability=0;
for iter=0:o.max_iter
    [f,g,c,G,gout]=evaluate(z);
    active=false(size(c));
    active(1:U)=c(1:U)>=-o.sinr_active_tol;
    active(U+1:end)=c(U+1:end)>=-o.power_active_tol;
    [lambda,Q,geometry]=multipliers_and_geometry(G,g,c,active);
    laggrad=g+G*lambda; project=@(v) v-Q*(Q'*v);
    gt=project(g);
    kkt=max([max([0;c]),norm(laggrad)/max(1,norm(g)), ...
        max(abs(lambda.*c))/max(1,norm(lambda,inf))]);
    hist.iter(end+1)=iter; hist.crlb(end+1)=gout.crlb;
    hist.objective(end+1)=f; hist.time_s(end+1)=toc(timer);
    hist.kkt(end+1)=kkt; hist.tangent_grad_norm(end+1)=norm(gt);
    hist.active_count(end+1)=nnz(active); hist.delta(end+1)=Delta;
    hist.accepted(end+1)=false; hist.ratio(end+1)=NaN;
    hist.step_type{end+1}='diagnostic'; hist.model_error(end+1)=NaN;
    hist.restoration_defect(end+1)=NaN; hist.restoration_iterations(end+1)=0;
    hist.lambda_min(end+1)=NaN; hist.meo_status{end+1}='not_called';
    hist.curvature_lower_bound(end+1)=NaN;
    hist.snapshot_iter(end+1)=iter; hist.z_snapshots(:,end+1)=z;
    hist.active_membership(:,end+1)=active;
    if size(hist.active_membership,2)==1
        hist.active_change_count(end+1)=0;
    else
        hist.active_change_count(end+1)=nnz(xor( ...
            hist.active_membership(:,end),hist.active_membership(:,end-1)));
    end
    hist.step_norm(end+1)=0; hist.cg_iter(end+1)=0;
    hist.tcg_final_res_norm(end+1)=NaN;

    % Identify a face numerically before relying on its local tangent model.
    if geometry.face_residual>o.face_tol && iter<o.max_iter
        [zs,~,cs,rest]=restore(z,z,active,zeros(size(c)), ...
            lambda,laggrad,inf,true);
        [ared,noise]=decrease(z,zs);
        if rest.face_ok && max([0;cs])<=o.feas_tol && ared>=-noise
            z=zs; hist.step_type{end}='face_restoration';
            hist.accepted(end)=true; hist.restoration_iterations(end)=rest.iterations;
            work.accepted_steps=work.accepted_steps+1;
            continue;
        end
    end

    Hv=@(v) project(hvp(z,project(v),lambda));
    % Certify a SUPERSPACE of the critical cone, not the smaller tangent
    % space formed by imposing equality on zero-multiplier constraints.
    strong=active & lambda.*vecnorm(G).' > o.dual_margin;
    if any(strong)
        Astrong=G(:,strong)./max(vecnorm(G(:,strong)),realmin);
        Qcert=orth(Astrong,o.rank_tol);
    else
        Qcert=zeros(dim,0);
    end
    if strcmp(o.curvature_space,'ambient'); Qcert=zeros(dim,0); end
    if strcmp(o.curvature_space,'active_tangent'); Qcert=Q; end
    project_cert=@(v) v-Qcert*(Qcert'*v);
    Hv_cert=@(v) project_cert(hvp(z,project_cert(v),lambda));
    candidate=kkt<=o.epsilon_KKT && geometry.face_residual<=o.face_tol;
    nc_direction=[];
    if candidate
        if isequal(z,meo_z)
            oracle=info;
        else
            oracle=struct('conditional_psd',false);
            internal_epsilon=o.epsilon_H;
            if isfinite(o.trusted_hvp_error_bound)
                internal_epsilon=internal_epsilon-o.trusted_hvp_error_bound;
            end
            if o.use_structured_psd
                work.structured_checks=work.structured_checks+1;
                oracle=costar_sosp_structured_psd(cfg,channels,smat, ...
                    real_vector_to_complex_matrix(z,N,U),lambda,internal_epsilon);
                work.hvp_calls=work.hvp_calls+oracle.hvp_calls;
                work.exact_hvp_calls=work.exact_hvp_calls+oracle.hvp_calls;
                work.structured_passes=work.structured_passes+oracle.conditional_psd;
            end
            if ~oracle.conditional_psd
                structured_info=oracle;
                work.meo_calls=work.meo_calls+1;
                call_probability=o.failure_probability/ ...
                    (work.meo_calls*(work.meo_calls+1));
                spent_failure_probability=spent_failure_probability+call_probability;
                work.norm_bound_builds=work.norm_bound_builds+1;
                L=costar_lagrangian_norm_bound(cfg,channels,smat,z,lambda);
                before_hvp=work.hvp_calls;
                oracle=costar_sosp_meo(Hv_cert,dim,project_cert,struct( ...
                    'epsilon_H',internal_epsilon,'norm_bound',L, ...
                    'failure_probability',call_probability, ...
                    'max_iter',o.meo_max_iter,'effective_dim',dim-size(Qcert,2), ...
                    'seed',o.seed+100003*work.meo_calls, ...
                    'record_history',o.meo_record_history, ...
                    'overlap_safety_factor',o.meo_overlap_safety_factor));
                oracle.ambient=strcmp(o.curvature_space,'ambient');
                oracle.structured_attempt=structured_info;
                work.meo_hvp_calls=work.meo_hvp_calls+work.hvp_calls-before_hvp;
            end
            info=oracle; meo_z=z;
        end
        lambda_min=oracle.rayleigh;
        hist.lambda_min(end)=lambda_min; hist.meo_status{end}=oracle.status;
        hist.curvature_lower_bound(end)=oracle.eigenvalue_lower_bound;
        if oracle.conditional_psd
            sol.conditional_sosp_reached=oracle.ambient || geometry.licq;
            if strcmp(o.curvature_space,'active_tangent') && ~oracle.ambient
                sol.conditional_sosp_reached=geometry.licq && ...
                    geometry.positive_multipliers;
            end
            sol.success=sol.conditional_sosp_reached;
            if sol.success
                sol.message='conditional approximate SOSP stopping contract passed';
            else
                sol.message='tangent curvature passed; critical-cone regularity unresolved';
            end
            break;
        elseif strcmp(oracle.status,'negative_curvature')
            % A negative direction in the larger certification space need
            % not obey weak active inequalities. Try its all-active tangent
            % projection; otherwise declare the cone search unresolved.
            nc_direction=project(oracle.negative_direction);
            if norm(nc_direction)>o.rank_tol
                nc_direction=nc_direction/norm(nc_direction);
            end
            if norm(nc_direction)<=o.rank_tol || ...
                    nc_direction'*Hv(nc_direction)>-o.epsilon_H/2
                sol.message='negative curvature in superspace; feasible cone direction unresolved';
                break;
            end
            work.ritz_negative_curvature_detections= ...
                work.ritz_negative_curvature_detections+1;
        else
            sol.message=['curvature unresolved: ' oracle.status]; break;
        end
    end
    if iter==o.max_iter; break; end

    model_hvp=Hv;
    if strcmpi(o.hessian_model,'gauss_newton_psd')
        work.gn_model_builds=work.gn_model_builds+1;
        [J,Bpsd]=costar_gauss_newton_model(cfg,smat, ...
            real_vector_to_complex_matrix(z,N,U),gout,o.gn_eigen_floor);
        model_hvp=@(v) project(gn_hvp(project(v),J,Bpsd));
    end
    if ~isempty(nc_direction)
        % The independent curvature oracle always uses the true Hessian,
        % including for tCG/GN ablations; an oracle escape uses its model.
        model_hvp=Hv;
        if g'*nc_direction>0; nc_direction=-nc_direction; end
        step=Delta*nc_direction; step_type='meo_negative_curvature'; tcg_nc=false;
    elseif norm(gt)>o.gradient_floor
        work.tcg_calls=work.tcg_calls+1;
        cgopts=struct('max_iter',o.cg_max_iter,'tol',o.cg_tol, ...
            'nc_tol',o.cg_nc_tol,'flat_tol',1e-15, ...
            'forcing_power',2,'forcing_constant',0.2, ...
            'absolute_floor',o.cg_absolute_floor, ...
            'detect_negative_curvature',o.enable_tcg_negative_curvature);
        % Regularization is used only in the inner solve, not in the
        % unshifted Lagrangian model or terminal curvature certificate.
        [step,cg]=trust_region_tcg(gt, ...
            @(v) model_hvp(v)+o.newton_regularization*o.epsilon_H*project(v),Delta,cgopts);
        work.tcg_iterations=work.tcg_iterations+cg.iter;
        hist.cg_iter(end)=cg.iter; hist.tcg_final_res_norm(end)=cg.final_res_norm;
        tcg_nc=cg.negative_curvature; step_type='accurate_tcg';
        work.tcg_negative_curvature_detections= ...
            work.tcg_negative_curvature_detections+tcg_nc;
    else
        sol.message='gradient numerically zero but full stopping gate unresolved'; break;
    end
    step=project(step); Hstep=model_hvp(step);
    hist.step_norm(end)=norm(step);
    linear=g'*step; quadratic=0.5*step'*Hstep;
    predicted=-(linear+quadratic);
    pred_noise=64*eps*(abs(linear)+abs(quadratic)+ ...
        norm(g)*norm(step)+0.5*norm(step)*norm(Hstep));
    if predicted<=pred_noise || ~isfinite(predicted)
        sol.message='predicted decrease is numerically unresolved'; break;
    end
    [trial,~,ct,rest]=restore(z,z+step,active,c, ...
        lambda,laggrad,predicted,false);
    [actual,noise]=decrease(z,trial);
    ratio=actual/predicted;
    accepted=max([0;ct])<=o.feas_tol && rest.model_ok && ...
        actual>noise && ratio>=o.eta_accept;
    hist.ratio(end)=ratio; hist.step_type{end}=step_type;
    hist.accepted(end)=accepted; hist.model_error(end)=abs(predicted-actual);
    hist.restoration_defect(end)=rest.model_defect;
    hist.restoration_iterations(end)=rest.iterations;
    if accepted
        z=trial; work.accepted_steps=work.accepted_steps+1;
        if tcg_nc
            work.tcg_negative_curvature_accepted=work.tcg_negative_curvature_accepted+1;
        elseif ~isempty(nc_direction)
            work.ritz_escape_accepted=work.ritz_escape_accepted+1;
        end
    else
        work.rejected_steps=work.rejected_steps+1;
        work.restoration_failures=work.restoration_failures+~rest.model_ok;
        work.precision_limited_trials=work.precision_limited_trials+(actual<=noise);
    end
    if ~rest.model_ok || max([0;ct])>o.feas_tol || ratio<o.eta_shrink
        Delta=o.gamma_shrink*Delta;
    elseif ratio>o.eta_expand && norm(step)>=0.9*Delta
        Delta=min(o.gamma_expand*Delta,o.delta_max);
    end
    if Delta<o.minimum_radius
        sol.message='numerical radius limit; no second-order success'; break;
    end
end

% Always recompute endpoint quantities; no stale pre-step residuals.
[~,g,c,G,gout]=evaluate(z);
active=false(size(c)); active(1:U)=c(1:U)>=-o.sinr_active_tol;
active(U+1:end)=c(U+1:end)>=-o.power_active_tol;
[lambda,~,geometry]=multipliers_and_geometry(G,g,c,active);
kkt=max([max([0;c]),norm(g+G*lambda)/max(1,norm(g)), ...
    max(abs(lambda.*c))/max(1,norm(lambda,inf))]);
same_point=isequal(meo_z,z);
if ~same_point || kkt>o.epsilon_KKT
    sol.success=false; sol.conditional_sosp_reached=false;
end
% A radius/budget exit immediately after a step must retain that endpoint.
if ~isequal(hist.z_snapshots(:,end),z)
    row=numel(hist.iter);
    names=fieldnames(hist);
    for jj=1:numel(names)
        name=names{jj}; value=hist.(name);
        if any(strcmp(name,{'z_snapshots','active_membership'})); continue; end
        if isvector(value) && numel(value)==row
            hist.(name)(end+1)=value(end);
        end
    end
    hist.iter(end)=hist.iter(end-1)+1; hist.snapshot_iter(end)=hist.iter(end);
    hist.time_s(end)=toc(timer); hist.crlb(end)=gout.crlb;
    hist.objective(end)=f_endpoint(z); hist.kkt(end)=kkt;
    [~,Qend]=multipliers_and_geometry(G,g,c,active);
    hist.tangent_grad_norm(end)=norm(g-Qend*(Qend'*g));
    hist.active_count(end)=nnz(active); hist.delta(end)=Delta;
    hist.active_change_count(end)=nnz(xor(active,hist.active_membership(:,end)));
    hist.active_membership(:,end+1)=active; hist.z_snapshots(:,end+1)=z;
    hist.accepted(end)=false; hist.ratio(end)=NaN;
    hist.step_type{end}='endpoint'; hist.step_norm(end)=0; hist.cg_iter(end)=0;
end
W=real_vector_to_complex_matrix(z,N,U); wcell=matrix_to_wcell(cfg,W);
sol.Wmat=W; sol.wcell=wcell; sol.z=z;
sol.metrics=evaluate_beamforming_solution(cfg,channels,smat,wcell);
sol.runtime=toc(timer); sol.lambda=lambda; sol.kkt_residual=kkt;
sol.stationarity=norm(g+G*lambda)/max(1,norm(g));
sol.primal_residual=max([0;c]); sol.lambda_min=lambda_min; sol.delta=Delta;
sol.work=work; sol.options=o; sol.phase2_method='Accuracy-controlled SOSP RTR';
sol.second_order_verified=false;
scope='positive-multiplier superspace; critical-cone interpretation requires correct identification';
if strcmp(o.curvature_space,'active_tangent')
    scope='all-active tangent; critical-cone interpretation requires identification and positive multipliers';
end
if isfield(info,'ambient') && info.ambient
    scope='ambient Lagrangian Hessian; sufficient for every critical cone without strict complementarity';
end
sol.certificate=struct('conditional_high_probability',sol.conditional_sosp_reached, ...
    'floating_point_validated',false,'geometry',geometry,'oracle',info, ...
    'oracle_is_at_final_point',same_point,'epsilon_KKT',o.epsilon_KKT, ...
    'epsilon_H',o.epsilon_H,'feasibility_tolerance',o.feas_tol, ...
    'identified_face_tolerance',o.face_tol,'total_failure_probability_bound', ...
    spent_failure_probability,'requested_failure_probability',o.failure_probability, ...
    'scope',scope,'assumptions',{{'smooth positive-definite FIM domain', ...
    'exact analytic Hessian and valid structured decomposition', ...
    'orthogonal projection and independent Gaussian starts for random fallback', ...
    'correct active identification only for nonambient critical-cone interpretation'}});
hist.work=work; hist.options=o; sol.hist=hist;

    function [fv,gv,cv,Gv,out]=evaluate(zv)
        work.derivative_calls=work.derivative_calls+1;
        Wv=real_vector_to_complex_matrix(zv,N,U);
        [fv,gcomplex,out]=compute_crlb_objective_grad(cfg,smat,Wv);
        gv=complex_matrix_to_real_vector(gcomplex);
        [cv,Gv]=costar_sosp_constraints(cfg,channels,Wv);
    end

    function h=hvp(zv,v,lv)
        work.hvp_calls=work.hvp_calls+1;
        work.exact_hvp_calls=work.exact_hvp_calls+1;
        Wv=real_vector_to_complex_matrix(zv,N,U);
        V=real_vector_to_complex_matrix(v,N,U);
        h=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
            cfg,channels,smat,Wv,V,lv));
    end

    function h=gn_hvp(v,J,Bpsd)
        work.hvp_calls=work.hvp_calls+1;
        h=J*(Bpsd*(J'*v))+o.gn_damping*v;
    end

    function fv=f_endpoint(zv)
        fv=compute_crlb_objective_grad(cfg,smat,real_vector_to_complex_matrix(zv,N,U));
    end

    function [lv,Qv,geo]=multipliers_and_geometry(Gv,gv,cv,av)
        lv=zeros(size(cv)); A=Gv(:,av);
        scales=max(vecnorm(A),realmin); As=A./scales;
        if isempty(A)
            Qv=zeros(dim,0); singular=inf; positive=true; dual_margin=inf;
        else
            [Qall,S,~]=svd(As,'econ'); sv=diag(S);
            rankv=nnz(sv>o.rank_tol); Qv=Qall(:,1:rankv);
            mu=lsqnonneg(As,-gv); lv(av)=mu./scales.';
            singular=min(sv); positive=all(mu>o.dual_margin);
            dual_margin=min(mu);
        end
        geo=struct('licq',size(Qv,2)==nnz(av) && singular>=o.rank_tol, ...
            'scaled_sigma_min',singular,'positive_multipliers',positive, ...
            'scaled_multiplier_min',dual_margin,'face_residual',max([0;abs(cv(av))]), ...
            'inactive_slack',min([inf;-cv(~av)]),'active',av);
    end

    function [zv,fv,cv,result]=restore(base,trial,av,target,lv,lg,pred,entry)
        zv=trial; result=struct('iterations',0,'model_defect',inf, ...
            'model_ok',false,'face_ok',false);
        tangent_step=trial-base;
        for correction=0:o.restoration_max_iter
            [fv,~,cv,Gv]=evaluate(zv);
            delta_c=cv-target;
            normal=zv-trial;
            weighted=sum(abs(lv(av)).*abs(delta_c(av)));
            result.model_defect=weighted+abs(lg'*normal);
            result.face_ok=max([0;abs(cv(av))])<=o.face_tol;
            result.model_ok=max([0;cv])<=o.feas_tol && ...
                result.model_defect<=o.restoration_fraction*pred && ...
                norm(normal)<=o.normal_second_order_cap*norm(tangent_step)^2;
            if (entry && result.face_ok && max([0;cv])<=o.feas_tol) || ...
                    (~entry && result.model_ok)
                return;
            end
            if correction==o.restoration_max_iter; return; end
            violated=~av & cv>0;
            ids=find(av | violated);
            if isempty(ids); return; end
            residual=cv(ids)-target(ids);
            residual(violated(ids))=cv(ids(violated(ids)))+o.feas_tol/10;
            A=Gv(:,ids); scales=max(vecnorm(A),realmin);
            As=A./scales; rs=residual./scales.';
            [Qa,Sa,Va]=svd(As,'econ'); sv=diag(Sa);
            keep=sv>o.rank_tol;
            correction_step=-Qa(:,keep)*((Va(:,keep)'*rs)./sv(keep));
            if entry
                cap=o.delta0/4;
            else
                cap=o.normal_second_order_cap*norm(tangent_step)^2;
            end
            if norm(correction_step)>cap
                correction_step=correction_step*cap/norm(correction_step);
            end
            zv=zv+correction_step;
            result.iterations=correction+1;
            work.restoration_iterations=work.restoration_iterations+1;
        end
    end

    function [actual,noise]=decrease(base,trial)
        Wb=real_vector_to_complex_matrix(base,N,U);
        DW=real_vector_to_complex_matrix(trial-base,N,U);
        [actual,noise]=costar_sosp_stable_decrease(cfg,smat,Wb,DW);
    end
end

function o=defaults(o)
names={'max_iter','epsilon_KKT','epsilon_H','failure_probability','meo_max_iter', ...
    'delta0','delta_max','minimum_radius','eta_accept','eta_shrink','eta_expand', ...
    'gamma_shrink','gamma_expand','feas_tol','face_tol','sinr_active_tol', ...
    'power_active_tol','dual_margin','rank_tol','restoration_max_iter', ...
    'restoration_fraction','normal_second_order_cap','cg_max_iter','cg_tol', ...
    'cg_nc_tol','cg_absolute_floor','gradient_floor','newton_regularization', ...
    'trusted_hvp_error_bound','seed','use_structured_psd','curvature_space', ...
    'hessian_model','enable_tcg_negative_curvature','gn_eigen_floor','gn_damping', ...
    'meo_record_history','meo_overlap_safety_factor'};
values={120,1e-5,1e-5,1e-3,inf, ...
    0.05,0.5,1e-12,0.05,0.25,0.75, ...
    0.5,1.8,1e-7,1e-7,1e-4, ...
    1e-5,1e-10,1e-10,20, ...
    0.05,20,120,1e-12, ...
    1e-10,1e-14,1e-13,2,NaN,2026,true,'ambient', ...
    'exact_lagrangian',true,0,1e-10,false,0.5};
for j=1:numel(names)
    if ~isfield(o,names{j}); o.(names{j})=values{j}; end
end
assert(o.epsilon_H>0 && o.epsilon_KKT>0 && o.failure_probability>0 && ...
    o.failure_probability<1);
assert(~isfinite(o.trusted_hvp_error_bound) || ...
    (o.trusted_hvp_error_bound>=0 && o.trusted_hvp_error_bound<o.epsilon_H));
assert(any(strcmp(o.curvature_space, ...
    {'ambient','positive_multiplier_superspace','active_tangent'})));
assert(any(strcmp(o.hessian_model,{'exact_lagrangian','gauss_newton_psd'})));
end
