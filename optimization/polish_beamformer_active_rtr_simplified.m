function [sol, hist] = polish_beamformer_active_rtr_simplified( ...
    cfg, channels, smat, init_wcell, opts)
%POLISH_BEAMFORMER_ACTIVE_RTR Active-set Riemannian trust-region polish.
%
% Active AP power balls are treated as a product of complex spheres. Nearly
% active SINR inequalities are added to the critical tangent space. Each
% subproblem is solved by Steihaug-Toint tCG using matrix-free Hessian-vector
% products of the original constrained Lagrangian. The default is the exact
% analytic action; the former central-difference model remains an ablation.

if nargin < 5
    opts = struct();
end
opts = fill_defaults(opts);

work.derivative_calls = 0;
work.hvp_calls = 0;
work.exact_hvp_calls = 0;
work.fd_hvp_calls = 0;
work.lanczos_calls = 0;
work.lanczos_hvp_calls = 0;
work.tcg_calls = 0;
work.tcg_iterations = 0;
work.tcg_negative_curvature_detections = 0;
work.tcg_negative_curvature_accepted = 0;
work.ritz_negative_curvature_detections = 0;
work.ritz_escape_candidates = 0;
work.ritz_escape_accepted = 0;
work.gn_model_builds = 0;
work.nnls_calls = 0;
work.multiplier_calls = 0;
work.normal_correction_steps = 0;
work.accepted_steps = 0;
work.rejected_steps = 0;

N = cfg.sys.M * cfg.sys.Nt;
U = cfg.sys.U;
dim = 2*N*U;
Wmat = vectorize_wcell_to_matrix(cfg,init_wcell);
Wmat = project_per_ap_power_matrix(cfg,Wmat);
z = complex_matrix_to_real_vector(Wmat);
Delta = opts.delta0;

hist = init_history();
sol.success = false;
sol.message = 'maximum RTR iterations reached';
lambda_min = NaN;
negative_curvature = inf;

timer = tic;
for iter = 0:opts.max_iter
    [f,grad_f,c,Gz,crlb,objective_info] = derivatives(z);
    active = select_active_constraints(c);
    lambda = estimate_multipliers(Gz,grad_f,active);
    lag_grad = grad_f + Gz*lambda;
    Q = active_normal_basis(Gz,active);
    project_tangent = @(v) v - Q*(Q'*v);
    rgrad = project_tangent(grad_f);
    stationarity = norm(lag_grad) / max(1,norm(grad_f));
    primal = max([0;c]);

    switch lower(string(opts.hessian_model))
        case {"exact_lagrangian","lagrangian_exact","exact_matrix_free"}
            Hv_raw = @(v) exact_lagrangian_hess_vec_real(z,v,lambda);
        case {"lagrangian_fd","full","full_lagrangian"}
            Hv_raw = @(v) lagrangian_hess_vec(z,v,lambda);
        case {"gauss_newton_psd","gn_psd","generalized_gauss_newton"}
            [gn_jacobian,gn_outer_hessian] = ...
                build_generalized_gauss_newton_model(z,objective_info);
            Hv_raw = @(v) generalized_gauss_newton_hess_vec( ...
                v,gn_jacobian,gn_outer_hessian);
        otherwise
            error('Unknown hessian_model: %s',opts.hessian_model);
    end
    Hv_tangent = @(v) project_tangent(Hv_raw(project_tangent(v)));

    switch lower(string(opts.curvature_policy))
        case "periodic"
            do_curvature_check = iter == 0 || ...
                mod(iter,opts.curvature_check_every) == 0 || ...
                stationarity <= opts.curvature_check_stationarity || ...
                iter == opts.max_iter;
        case "lazy"
            do_curvature_check = iter == opts.max_iter || ...
                norm(rgrad) <= opts.grad_floor || ...
                (primal <= opts.feas_tol && stationarity <= ...
                max(opts.stationarity_tol,opts.curvature_check_stationarity));
        case "terminal"
            % Certify curvature only when the current point is an actual
            % first-order stopping candidate. tCG remains responsible for
            % negative-curvature detection during ordinary RTR iterations.
            do_curvature_check = iter == opts.max_iter || ...
                norm(rgrad) <= opts.grad_floor || ...
                (primal <= opts.feas_tol && ...
                stationarity <= opts.stationarity_tol);
        case "none"
            do_curvature_check = false;
            negative_curvature = 0;
        otherwise
            error('Unknown curvature_policy: %s',opts.curvature_policy);
    end
    curvature_info = struct('ritz_vector',[],'iter',0,'ritz_values',[],'T',[]);
    ritz_negative_detected = false;
    if do_curvature_check
        work.lanczos_calls = work.lanczos_calls+1;
        hvp_before_lanczos = work.hvp_calls;
        lanczos_opts.max_iter = min(opts.lanczos_max_iter,dim);
        lanczos_opts.tol = opts.lanczos_tol;
        lanczos_opts.seed = opts.seed + iter;
        [lambda_min,curvature_info] = estimate_min_curvature_lanczos( ...
            Hv_tangent,dim,lanczos_opts);
        work.lanczos_hvp_calls = work.lanczos_hvp_calls + ...
            (work.hvp_calls-hvp_before_lanczos);
        negative_curvature = max(0,-lambda_min);
        ritz_negative_detected = negative_curvature > opts.curvature_tol;
        if ritz_negative_detected
            work.ritz_negative_curvature_detections = ...
                work.ritz_negative_curvature_detections+1;
        end
    end

    hist = store_history(hist,iter,f,crlb,toc(timer),c,stationarity,primal, ...
        lambda_min,negative_curvature,Delta,z,active,false,NaN,0,'diagnostic');
    hist.tangent_grad_norm(end) = norm(rgrad);
    hist.lanczos_checked(end) = do_curvature_check;
    hist.ritz_negative_detected(end) = ritz_negative_detected;

    if primal <= opts.feas_tol && stationarity <= opts.stationarity_tol && ...
            negative_curvature <= opts.curvature_tol
        sol.success = true;
        sol.message = 'KKT and weak second-order stationarity reached';
        break;
    end
    if iter == opts.max_iter
        break;
    end

    tcg_opts.max_iter = min(opts.cg_max_iter,dim);
    tcg_opts.tol = opts.cg_tol;
    tcg_opts.nc_tol = opts.cg_nc_tol;
    tcg_opts.flat_tol = opts.cg_flat_tol;
    tcg_opts.forcing_power = opts.cg_forcing_power;
    tcg_opts.forcing_constant = opts.cg_forcing_constant;
    tcg_opts.absolute_floor = opts.cg_absolute_floor;
    tcg_opts.detect_negative_curvature = ...
        opts.enable_tcg_negative_curvature;
    tcg_opts.verbose = false;

    used_tcg = false;
    used_ritz_escape = false;
    if norm(rgrad) > opts.grad_floor
        work.tcg_calls = work.tcg_calls+1;
        [step,cgstats] = trust_region_tcg(rgrad,Hv_tangent,Delta,tcg_opts);
        work.tcg_iterations = work.tcg_iterations+cgstats.iter;
        used_tcg = true;
        if cgstats.negative_curvature
            work.tcg_negative_curvature_detections = ...
                work.tcg_negative_curvature_detections+1;
        end
        step_type = 'tcg';
    elseif opts.enable_ritz_escape && ...
            negative_curvature > opts.curvature_tol && ...
            ~isempty(curvature_info.ritz_vector)
        direction = project_tangent(curvature_info.ritz_vector);
        direction = direction/max(norm(direction),realmin);
        step = Delta*direction;
        if grad_f'*step > 0
            step = -step;
        end
        cgstats = empty_cgstats('lanczos_negative_curvature');
        cgstats.negative_curvature = true;
        step_type = 'lanczos_nc';
        used_ritz_escape = true;
        work.ritz_escape_candidates = work.ritz_escape_candidates+1;
    else
        sol.success = primal <= opts.feas_tol && stationarity <= opts.stationarity_tol;
        sol.message = 'no first- or second-order RTR direction found';
        break;
    end

    step = project_tangent(step);
    if opts.reuse_tcg_model && strcmp(step_type,'tcg') && ...
            isfinite(cgstats.model_decrease)
        predicted = cgstats.model_decrease;
    else
        Hstep = Hv_tangent(step);
        predicted = -(grad_f'*step + 0.5*step'*Hstep);
    end
    ztrial = retract_active_power(z,step,c);
    if opts.fuse_trial_evaluation
        [ztrial,ftrial,ctrial] = ...
            correct_sinr_retraction_evaluated(ztrial,c,Delta);
    else
        if opts.enable_sinr_correction
            ztrial = correct_sinr_retraction(ztrial,c,Delta);
        end
        [ftrial,~,ctrial,~] = derivatives(ztrial);
    end
    actual = f-ftrial;
    feasible_trial = max([0;ctrial]) <= opts.feas_tol;
    if predicted > opts.pred_tol && isfinite(actual)
        ratio = actual/predicted;
    else
        ratio = -inf;
    end
    accepted = feasible_trial && actual > opts.min_actual_decrease && ...
        ratio >= opts.eta_accept;

    boundary_step = norm(step) >= opts.boundary_tol*Delta;
    if ~feasible_trial || ratio < opts.eta_shrink
        Delta = max(opts.gamma_shrink*Delta,opts.delta_min);
    elseif ratio > opts.eta_expand && boundary_step
        Delta = min(opts.gamma_expand*Delta,opts.delta_max);
    end
    if accepted
        z = ztrial;
        work.accepted_steps = work.accepted_steps+1;
        if used_tcg && cgstats.negative_curvature
            work.tcg_negative_curvature_accepted = ...
                work.tcg_negative_curvature_accepted+1;
        end
        if used_ritz_escape
            work.ritz_escape_accepted = work.ritz_escape_accepted+1;
        end
    else
        work.rejected_steps = work.rejected_steps+1;
    end

    hist.accepted(end) = accepted;
    hist.ratio(end) = ratio;
    hist.step_norm(end) = norm(step);
    hist.step_type{end} = step_type;
    hist.cg_iter(end) = cgstats.iter;
    hist.cg_negative_curvature(end) = cgstats.negative_curvature;
    hist.tcg_called(end) = used_tcg;
    hist.tcg_reason{end} = cgstats.reason;
    hist.tcg_final_res_norm(end) = cgstats.final_res_norm;
    hist.tcg_min_curvature_est(end) = cgstats.min_curvature_est;
    hist.tcg_nc_accepted(end) = ...
        used_tcg && cgstats.negative_curvature && accepted;
    hist.ritz_escape_candidate(end) = used_ritz_escape;
    hist.ritz_escape_accepted(end) = used_ritz_escape && accepted;

    if Delta <= opts.delta_min && ~accepted
        sol.message = 'trust-region radius reached minimum';
        break;
    end
end
runtime = toc(timer);

Wmat = real_vector_to_complex_matrix(z,N,U);
wcell = matrix_to_wcell(cfg,Wmat);
metrics = evaluate_beamforming_solution(cfg,channels,smat,wcell);
sol.wcell = wcell;
sol.Wmat = Wmat;
sol.z = z;
sol.metrics = metrics;
sol.runtime = runtime;
sol.lambda = lambda;
sol.lambda_min = lambda_min;
sol.negative_curvature = negative_curvature;
sol.stationarity = stationarity;
sol.primal_residual = primal;
sol.delta = Delta;
sol.work = work;
sol.options = opts;
hist.work = work;
hist.options = opts;
sol.hist = hist;

    function [f,grad_f,c,Gz,crlbq,gout] = derivatives(zq)
        work.derivative_calls = work.derivative_calls+1;
        Wq = real_vector_to_complex_matrix(zq,N,U);
        [f,g_complex,gout] = compute_crlb_objective_grad(cfg,smat,Wq);
        crlbq = gout.crlb;
        [c,Gcomplex] = compute_constraints_grad(cfg,channels,smat,Wq);
        grad_f = complex_matrix_to_real_vector(g_complex);
        Gz = zeros(dim,numel(c));
        for jj = 1:numel(c)
            Gz(:,jj) = complex_matrix_to_real_vector(Gcomplex(:,:,jj));
        end
    end

    function active = select_active_constraints(cq)
        active = false(size(cq));
        active(1:U) = cq(1:U) >= -opts.sinr_active_tol;
        active(U+1:end) = cq(U+1:end) >= -opts.power_active_tol;
    end

    function Q = active_normal_basis(Gq,active)
        basis = Gq(:,active);
        if opts.project_phase_directions
            Wq = real_vector_to_complex_matrix(z,N,U);
            phase_basis = zeros(dim,U);
            for uu = 1:U
                phase_direction = zeros(N,U);
                phase_direction(:,uu) = 1i*Wq(:,uu);
                phase_basis(:,uu) = ...
                    complex_matrix_to_real_vector(phase_direction);
            end
            basis = [basis phase_basis];
        end
        if ~isempty(basis)
            Q = orth(basis,opts.rank_tol);
        else
            Q = zeros(dim,0);
        end
    end

    function lambdaq = estimate_multipliers(Gq,gq,active)
        lambdaq = zeros(size(Gq,2),1);
        if any(active)
            work.multiplier_calls = work.multiplier_calls+1;
            Aactive = Gq(:,active);
            switch lower(string(opts.multiplier_mode))
                case "nnls"
                    work.nnls_calls = work.nnls_calls+1;
                    lambdaq(active) = lsqnonneg(Aactive,-gq);
                case "projected_ls"
                    gram = Aactive'*Aactive;
                    regularization = opts.multiplier_regularization* ...
                        max(1,trace(gram)/max(size(gram,1),1));
                    lambdaq(active) = max(0, ...
                        -(gram+regularization*eye(size(gram)))\ ...
                        (Aactive'*gq));
                case "zero"
                    % Diagnostic only: ignores active constraint curvature.
                otherwise
                    error('Unknown multiplier_mode: %s',opts.multiplier_mode);
            end
        end
    end

    function Hv = lagrangian_hess_vec(zq,v,lambdaq)
        work.hvp_calls = work.hvp_calls+1;
        work.fd_hvp_calls = work.fd_hvp_calls+1;
        vnorm = norm(v);
        if vnorm <= realmin
            Hv = zeros(size(v));
            return;
        end
        displacement = opts.fd_rel_step*max(1,norm(zq));
        if opts.fd_trust_scaling
            displacement = min(displacement, ...
                opts.fd_trust_fraction*max(Delta,opts.delta_min));
        end
        displacement = max(displacement,opts.fd_min_displacement);
        h = displacement/vnorm;
        [~,gp,~,Gp] = derivatives(zq+h*v);
        [~,gm,~,Gm] = derivatives(zq-h*v);
        Hv = ((gp+Gp*lambdaq)-(gm+Gm*lambdaq))/(2*h);
        Hv = real(Hv(:));
    end

    function Hv = exact_lagrangian_hess_vec_real(zq,v,lambdaq)
        work.hvp_calls = work.hvp_calls+1;
        work.exact_hvp_calls = work.exact_hvp_calls+1;
        Wq = real_vector_to_complex_matrix(zq,N,U);
        Vq = real_vector_to_complex_matrix(real(v(:)),N,U);
        Hcomplex = exact_lagrangian_hess_vec( ...
            cfg,channels,smat,Wq,Vq,lambdaq);
        Hv = complex_matrix_to_real_vector(Hcomplex);
        Hv = real(Hv(:));
    end

    function [J,Bpsd] = build_generalized_gauss_newton_model(zq,gout)
        % The CRLB objective is a scalar function of the three quadratic
        % sufficient statistics (g_a,g_b,g_ab). Its generalized
        % Gauss-Newton model keeps J' * B_+ * J and omits both the
        % second derivatives of those statistics and active-constraint
        % curvature. B_+ is the PSD projection of the 3-by-3 outer Hessian.
        work.gn_model_builds = work.gn_model_builds+1;
        Wq = real_vector_to_complex_matrix(zq,N,U);
        beta2 = cfg.radio.effective_bandwidth_sq(:);
        matrices = {smat.Aa,smat.Ab,smat.Aab};
        J = zeros(dim,3);
        for stat_idx = 1:3
            gradient_stat = zeros(N,U);
            Astat = matrices{stat_idx};
            for user_idx = 1:U
                gradient_stat(:,user_idx) = 2*beta2(user_idx)* ...
                    Astat*Wq(:,user_idx);
            end
            J(:,stat_idx) = complex_matrix_to_real_vector(gradient_stat);
        end

        a = max(real(gout.ga),realmin);
        b = max(real(gout.gb),realmin);
        cstat = real(gout.gab);
        s = max(a+b,realmin);
        d = max(a*b-cstat^2,realmin);
        outer_gradient = [1/s-b/d; 1/s-a/d; 2*cstat/d];
        Bphi = [ ...
            -1/s^2+b^2/d^2, -1/s^2+cstat^2/d^2, -2*b*cstat/d^2; ...
            -1/s^2+cstat^2/d^2, -1/s^2+a^2/d^2, -2*a*cstat/d^2; ...
            -2*b*cstat/d^2, -2*a*cstat/d^2, 2/d+4*cstat^2/d^2];
        if strcmpi(gout.objective_type,'crlb_ratio')
            Bouter = gout.ratio*(Bphi+outer_gradient*outer_gradient');
        else
            Bouter = Bphi;
        end
        Bouter = real((Bouter+Bouter')/2);
        [V,Douter] = eig(Bouter,'vector');
        scale = max([1;abs(Douter(:))]);
        clipped = max(real(Douter),opts.gn_eigen_floor*scale);
        Bpsd = real(V*diag(clipped)*V');
        Bpsd = (Bpsd+Bpsd')/2;
    end

    function Hv = generalized_gauss_newton_hess_vec(v,J,Bpsd)
        work.hvp_calls = work.hvp_calls+1;
        Hv = J*(Bpsd*(J'*real(v(:)))) + opts.gn_damping*real(v(:));
        Hv = real(Hv(:));
    end

    function znew = retract_active_power(zq,step,cq)
        Wq = real_vector_to_complex_matrix(zq+step,N,U);
        Nt = cfg.sys.Nt;
        for mm = 1:cfg.sys.M
            rows = (mm-1)*Nt+(1:Nt);
            block = Wq(rows,:);
            power = real(sum(abs(block(:)).^2));
            is_active_power = cq(U+mm) >= -opts.power_active_tol;
            if is_active_power && power > 0
                block = block*sqrt(cfg.power.Pmax/power);
            elseif power > cfg.power.Pmax
                block = block*sqrt(cfg.power.Pmax/power);
            end
            Wq(rows,:) = block;
        end
        znew = complex_matrix_to_real_vector(Wq);
    end

    function zq = correct_sinr_retraction(zq,cbase,trust_radius)
        % Newton-like normal correction for the intersection of the active
        % power spheres and near-active SINR surfaces. This is part of the
        % retraction, not a separate optimization or line-search step.
        for correction_iter = 1:opts.sinr_retraction_max_iter
            [~,~,cq,Gq] = derivatives(zq);
            idx = find(cq(1:U) > -opts.sinr_retraction_margin);
            if isempty(idx)
                return;
            end

            active_power = cbase(U+1:end) >= -opts.power_active_tol;
            if any(active_power)
                Qpower = orth(Gq(:,U+find(active_power)),opts.rank_tol);
                project_power_tangent = @(v) v-Qpower*(Qpower'*v);
            else
                project_power_tangent = @(v) v;
            end
            A = project_power_tangent(Gq(:,idx));
            residual = cq(idx)+opts.sinr_retraction_margin;
            normal_matrix = A'*A + opts.sinr_retraction_regularization*eye(numel(idx));
            correction = -A*(normal_matrix\residual);
            correction = project_power_tangent(correction);
            correction_cap = opts.sinr_retraction_max_fraction*trust_radius;
            if norm(correction) > correction_cap
                correction = correction*(correction_cap/norm(correction));
            end
            zq = retract_active_power(zq,correction,cbase);
            work.normal_correction_steps = work.normal_correction_steps+1;
        end
    end

    function [zq,fq,cq] = correct_sinr_retraction_evaluated( ...
            zq,cbase,trust_radius)
        % Shares the final trial evaluation with the acceptance test. When
        % no SINR correction is needed, this removes one full derivative
        % evaluation compared with the reference implementation.
        [fq,~,cq,Gq] = derivatives(zq);
        if ~opts.enable_sinr_correction
            return;
        end
        for correction_iter = 1:opts.sinr_retraction_max_iter
            idx = find(cq(1:U) > -opts.sinr_retraction_margin);
            if isempty(idx)
                return;
            end
            active_power = cbase(U+1:end) >= -opts.power_active_tol;
            if any(active_power)
                Qpower = orth(Gq(:,U+find(active_power)),opts.rank_tol);
                project_power_tangent = @(v) v-Qpower*(Qpower'*v);
            else
                project_power_tangent = @(v) v;
            end
            A = project_power_tangent(Gq(:,idx));
            residual = cq(idx)+opts.sinr_retraction_margin;
            normal_matrix = A'*A + ...
                opts.sinr_retraction_regularization*eye(numel(idx));
            correction = -A*(normal_matrix\residual);
            correction = project_power_tangent(correction);
            correction_cap = ...
                opts.sinr_retraction_max_fraction*trust_radius;
            if norm(correction) > correction_cap
                correction = correction*(correction_cap/norm(correction));
            end
            zq = retract_active_power(zq,correction,cbase);
            work.normal_correction_steps = work.normal_correction_steps+1;
            [fq,~,cq,Gq] = derivatives(zq);
        end
    end
end

function hist = init_history()
hist.iter=[]; hist.objective=[]; hist.primal_residual=[];
hist.crlb=[]; hist.time_s=[];
hist.stationarity=[]; hist.lambda_min=[]; hist.negative_curvature=[];
hist.delta=[]; hist.snapshot_iter=[]; hist.z_snapshots=[];
hist.active_count=[]; hist.accepted=[]; hist.ratio=[]; hist.step_norm=[];
hist.step_type={}; hist.cg_iter=[]; hist.cg_negative_curvature=[];
hist.tcg_called=[]; hist.tcg_reason={}; hist.tcg_final_res_norm=[];
hist.tcg_min_curvature_est=[];
hist.tcg_nc_accepted=[]; hist.lanczos_checked=[];
hist.ritz_negative_detected=[]; hist.ritz_escape_candidate=[];
hist.ritz_escape_accepted=[];
hist.max_constraint=[];
hist.tangent_grad_norm=[];
hist.active_membership=[];
hist.active_change_count=[];
end

function hist = store_history(hist,iter,f,crlb,time_s,c,stationarity,primal,lambda_min, ...
    negative_curvature,Delta,z,active,accepted,ratio,step_norm,step_type)
hist.iter(end+1)=iter; hist.objective(end+1)=f;
hist.crlb(end+1)=crlb; hist.time_s(end+1)=time_s;
hist.primal_residual(end+1)=primal; hist.stationarity(end+1)=stationarity;
hist.lambda_min(end+1)=lambda_min;
hist.negative_curvature(end+1)=negative_curvature;
hist.delta(end+1)=Delta; hist.snapshot_iter(end+1)=iter;
if isempty(hist.z_snapshots) || size(hist.z_snapshots,1)==numel(z)
    hist.z_snapshots(:,end+1)=z;
end
hist.active_count(end+1)=nnz(active);
hist.active_membership(:,end+1)=logical(active(:));
if size(hist.active_membership,2)==1
    hist.active_change_count(end+1)=0;
else
    hist.active_change_count(end+1)=nnz(xor( ...
        hist.active_membership(:,end),hist.active_membership(:,end-1)));
end
hist.accepted(end+1)=accepted; hist.ratio(end+1)=ratio;
hist.step_norm(end+1)=step_norm; hist.step_type{end+1}=step_type;
hist.cg_iter(end+1)=0; hist.cg_negative_curvature(end+1)=false;
hist.tcg_called(end+1)=false; hist.tcg_reason{end+1}='not_called';
hist.tcg_final_res_norm(end+1)=NaN;
hist.tcg_min_curvature_est(end+1)=NaN;
hist.tcg_nc_accepted(end+1)=false;
hist.lanczos_checked(end+1)=false;
hist.ritz_negative_detected(end+1)=false;
hist.ritz_escape_candidate(end+1)=false;
hist.ritz_escape_accepted(end+1)=false;
hist.max_constraint(end+1)=max(c);
hist.tangent_grad_norm(end+1)=NaN;
end

function stats = empty_cgstats(reason)
stats.iter=0; stats.reason=reason; stats.negative_curvature=false;
stats.near_zero_curvature=false; stats.hit_boundary=true;
stats.min_curvature_est=NaN; stats.final_res_norm=NaN;
stats.model_value=NaN; stats.model_decrease=NaN; stats.step_norm=0;
end

function opts = fill_defaults(opts)
if ~isfield(opts,'max_iter'); opts.max_iter=120; end
if ~isfield(opts,'delta0'); opts.delta0=5e-2; end
if ~isfield(opts,'delta_max'); opts.delta_max=5e-1; end
if ~isfield(opts,'delta_min'); opts.delta_min=1e-9; end
if ~isfield(opts,'eta_accept'); opts.eta_accept=0.05; end
if ~isfield(opts,'eta_shrink'); opts.eta_shrink=0.25; end
if ~isfield(opts,'eta_expand'); opts.eta_expand=0.75; end
if ~isfield(opts,'gamma_shrink'); opts.gamma_shrink=0.5; end
if ~isfield(opts,'gamma_expand'); opts.gamma_expand=1.8; end
if ~isfield(opts,'boundary_tol'); opts.boundary_tol=0.9; end
if ~isfield(opts,'pred_tol'); opts.pred_tol=1e-12; end
if ~isfield(opts,'min_actual_decrease'); opts.min_actual_decrease=1e-12; end
if ~isfield(opts,'sinr_active_tol'); opts.sinr_active_tol=1e-3; end
if ~isfield(opts,'power_active_tol'); opts.power_active_tol=1e-5; end
if ~isfield(opts,'sinr_retraction_margin'); opts.sinr_retraction_margin=1e-5; end
if ~isfield(opts,'sinr_retraction_max_iter'); opts.sinr_retraction_max_iter=4; end
if ~isfield(opts,'sinr_retraction_regularization'); opts.sinr_retraction_regularization=1e-12; end
if ~isfield(opts,'sinr_retraction_max_fraction'); opts.sinr_retraction_max_fraction=1.0; end
if ~isfield(opts,'rank_tol'); opts.rank_tol=1e-10; end
if ~isfield(opts,'fd_rel_step'); opts.fd_rel_step=1e-5; end
if ~isfield(opts,'fd_trust_scaling'); opts.fd_trust_scaling=true; end
if ~isfield(opts,'fd_trust_fraction'); opts.fd_trust_fraction=0.1; end
if ~isfield(opts,'fd_min_displacement'); opts.fd_min_displacement=1e-8; end
if ~isfield(opts,'cg_max_iter'); opts.cg_max_iter=30; end
if ~isfield(opts,'cg_tol'); opts.cg_tol=1e-3; end
if ~isfield(opts,'cg_nc_tol'); opts.cg_nc_tol=1e-8; end
if ~isfield(opts,'cg_flat_tol'); opts.cg_flat_tol=1e-12; end
if ~isfield(opts,'cg_forcing_power'); opts.cg_forcing_power=[]; end
if ~isfield(opts,'cg_forcing_constant'); opts.cg_forcing_constant=1; end
if ~isfield(opts,'cg_absolute_floor'); opts.cg_absolute_floor=0; end
if ~isfield(opts,'project_phase_directions')
    opts.project_phase_directions=false;
end
if ~isfield(opts,'grad_floor'); opts.grad_floor=1e-10; end
if ~isfield(opts,'lanczos_max_iter'); opts.lanczos_max_iter=20; end
if ~isfield(opts,'lanczos_tol'); opts.lanczos_tol=1e-8; end
if ~isfield(opts,'curvature_check_every'); opts.curvature_check_every=10; end
if ~isfield(opts,'curvature_check_stationarity'); opts.curvature_check_stationarity=2e-3; end
if ~isfield(opts,'feas_tol'); opts.feas_tol=1e-6; end
if ~isfield(opts,'stationarity_tol'); opts.stationarity_tol=1e-3; end
if ~isfield(opts,'curvature_tol'); opts.curvature_tol=1e-5; end
if ~isfield(opts,'seed'); opts.seed=2026; end
if ~isfield(opts,'curvature_policy'); opts.curvature_policy='periodic'; end
if ~isfield(opts,'reuse_tcg_model'); opts.reuse_tcg_model=false; end
if ~isfield(opts,'fuse_trial_evaluation'); opts.fuse_trial_evaluation=false; end
if ~isfield(opts,'enable_sinr_correction'); opts.enable_sinr_correction=true; end
if ~isfield(opts,'multiplier_mode'); opts.multiplier_mode='nnls'; end
if ~isfield(opts,'hessian_model'); opts.hessian_model='exact_lagrangian'; end
if ~isfield(opts,'enable_tcg_negative_curvature')
    opts.enable_tcg_negative_curvature=true;
end
if ~isfield(opts,'enable_ritz_escape'); opts.enable_ritz_escape=false; end
if ~isfield(opts,'gn_eigen_floor'); opts.gn_eigen_floor=0; end
if ~isfield(opts,'gn_damping'); opts.gn_damping=1e-10; end
if ~isfield(opts,'multiplier_regularization')
    opts.multiplier_regularization=1e-10;
end
end
