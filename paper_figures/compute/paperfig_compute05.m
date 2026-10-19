function archive=paperfig_compute05(options)
%PAPERFIG_COMPUTE05 Controlled local Newton test for the conditional rate result.
% The reference starts from the CURRENT solve_costar endpoint. The local
% experiment uses the same geometry, exact HVP, tCG and level-set normal
% correction, with quadratic forcing explicitly imposed. It is not a claim
% that the default fixed-relative tCG tolerance yields a quadratic trajectory.
if nargin<1; options=struct(); end
s=paperfig_setup(); o=paperfig_options(options);
if ~isfield(options,'seeds'); o.seeds=o.seeds(1); end
% The default experiment measures a controlled residual map; regularity
% diagnostics are recorded independently of the observed contraction order.
if ~isfield(o,'local_spec'); o.local_spec=s.default; end
if ~isfield(o,'require_local_assumptions'); o.require_local_assumptions=false; end
if ~isfield(o,'local_errors'); o.local_errors=logspace(log10(5e-6),log10(2.5e-5),11); end
if ~isfield(o,'local_face_constant'); o.local_face_constant=0.01; end
if ~isfield(o,'local_face_ratio_cap'); o.local_face_ratio_cap=1e3; end
folder=fullfile(s.results_dir,'paper_fig05'); if ~isfolder(folder); mkdir(folder); end
archive=fullfile(folder,'data.mat'); signature=paperfig_signature(); audit_signature=paperfig_audit_signature(5);
d=struct('figure_id',5,'signature',signature,'audit_signature',audit_signature,'options',o,'seeds',o.seeds, ...
    'trials',{cell(1,numel(o.seeds))},'done',false(1,numel(o.seeds)),'complete',false, ...
    'experiment','Empirical controlled horizontal Newton residual map with quadratic forcing; regularity diagnostics recorded separately');
if o.resume && isfile(archive)
    old=load(archive);
    if isfield(old,'signature') && strcmp(old.signature,signature) && ...
            isfield(old,'audit_signature') && strcmp(old.audit_signature,audit_signature) && isequaln(old.options,o); d=old; end
end
for t=1:numel(o.seeds)
    if d.done(t); continue; end
    key='local_empirical'; if isequaln(o.local_spec,s.default); key='default'; end
    [p,run]=paperfig_case(key,o.local_spec,o.seeds(t),'costar',o);
    d.trials{t}=local_trial(p,run,o); d.done(t)=true; simtwc_save_atomic(archive,d);
end
x=[]; y=[]; rows=cell(0,15);
for t=1:numel(d.trials)
    a=d.trials{t}; use=a.fit_mask;
    x=[x a.r0(use)]; y=[y a.r1(use)]; %#ok<AGROW>
    for j=1:numel(a.r0)
        rows(end+1,:)={a.seed,a.errors(j),a.r0(j),a.r1(j), ...
            a.r1(j)/a.r0(j)^2,a.accepted(j),a.feasible(j),a.interior(j), ...
            a.stable(j),a.linear_residual(j),a.normal_ratio(j),a.forcing_ratio(j), ...
            a.face_residual(j),a.face_ratio(j),use(j)}; %#ok<AGROW>
    end
end
assert(numel(x)>=3,'Fewer than three admissible local Newton samples; all records were saved.');
X=[ones(numel(x),1),log10(x(:))]; b=X\log10(y(:)); residual=log10(y(:))-X*b;
covariance=(residual'*residual/max(1,numel(x)-2))*pinv(X'*X);
d.pooled_slope=b(2); d.pooled_slope_ci95=b(2)+[-1 1]*1.96*sqrt(max(0,covariance(2,2)));
d.r0=x; d.r1=y; d.quadratic_cv=std(y./x.^2)/mean(y./x.^2);
d.complete=all(d.done); d.completed_at=char(datetime('now')); simtwc_save_atomic(archive,d);
writetable(cell2table(rows,'VariableNames',{'seed','relative_perturbation','r_k','r_k_plus_1', ...
    'quadratic_ratio','accepted','feasible','interior','stable_working_set', ...
    'linear_residual','normal_ratio','forcing_ratio','face_residual','face_ratio', ...
    'used_in_fit'}),fullfile(folder,'local_steps.csv'));
writetable(table(d.pooled_slope,d.pooled_slope_ci95(1),d.pooled_slope_ci95(2), ...
    d.quadratic_cv,numel(x),'VariableNames',{'order','CI95_low','CI95_high','quadratic_CV','samples'}), ...
    fullfile(folder,'fit_summary.csv'));
paperfig_reference_conditions(d,folder);
fprintf('[Paper Fig5] Controlled local order %.4f, CV %.3g, %d qualified samples\n',d.pooled_slope,d.quadratic_cv,numel(x));
end

function a=local_trial(p,run,o)
z=run.solution.z; active=run.solution.active; U=p.cfg.sys.U;
p.paperfig_working_thresholds=run.solution.working_thresholds;
% Polish only inside the declared conditional local Newton experiment.
z=level_restore(p,z,active,zeros(nnz(active),1));
polish=cell(1,15);
for k=1:15
    st=state(p,z,active); if norm(st.gt)<1e-12; break; end
    q=one_step(p,z,active,0.5); polish{k}=q;
    if ~q.accepted || q.r1>=q.r0; break; end
    z=q.z1;
end
zstar=level_restore(p,z,active,zeros(nnz(active),1)); st=state(p,zstar,active);
h=run.history.rtr; candidates=[]; distances=[];
for k=1:size(h.z_snapshots,2)
    dist=norm(h.z_snapshots(:,k)-zstar)/max(1,norm(zstar));
    if dist>2*max(o.local_errors); candidates(end+1)=k; distances(end+1)=dist; end %#ok<AGROW>
end
if isempty(candidates)
    stream=RandStream('mt19937ar','Seed',p.seed+18000); direction=randn(stream,numel(z),1);
else
    [~,j]=min(distances); direction=h.z_snapshots(:,candidates(j))-zstar;
end
direction=st.project(direction);
% Restrict the perturbation to identifiable horizontal curvature directions.
direction=st.project(hvp(p,zstar,direction,st.lambda)); direction=direction/max(norm(direction),realmin);
n=numel(o.local_errors);
reference_check=reference_audit(p,zstar,active,st);
disp(reference_check);
local_assumptions=reference_check.licq && reference_check.positive_multipliers && ...
    reference_check.horizontal_positive_definite && reference_check.tangent_residual<=1e-9 && ...
    reference_check.face_residual<=1e-11;
reference_check.all_local_assumptions=local_assumptions;
if o.require_local_assumptions
    assert(local_assumptions,'Local reference fails the requested regularity/accuracy tests.');
end
a=struct('status','ok','seed',p.seed,'errors',o.local_errors,'reference_z',zstar, ...
    'reference_check',reference_check, ...
    'reference_tangent_residual',norm(st.gt),'reference_solver_success',run.solution.success, ...
    'reference_algorithm',run.solution.algorithm_version,'polish',{polish}, ...
    'r0',nan(1,n),'r1',nan(1,n),'accepted',false(1,n),'feasible',false(1,n), ...
    'interior',false(1,n),'stable',false(1,n),'linear_residual',nan(1,n), ...
    'normal_ratio',nan(1,n),'forcing_ratio',nan(1,n),'face_residual',nan(1,n), ...
    'face_ratio',nan(1,n),'fit_mask',false(1,n),'steps',{cell(1,n)});
for j=1:n
    perturbation=o.local_errors(j)*max(1,norm(zstar));
    target=-o.local_face_constant*perturbation^2*ones(nnz(active),1);
    z0=level_restore(p,zstar+perturbation*direction,active,target);
    q=one_step(p,z0,active,0.5); a.steps{j}=q;
    a.r0(j)=q.r0; a.r1(j)=q.r1; a.accepted(j)=q.accepted; a.feasible(j)=q.feasible;
    a.interior(j)=q.interior; a.stable(j)=q.stable;
    a.linear_residual(j)=q.linear_residual; a.normal_ratio(j)=q.normal_ratio;
    a.forcing_ratio(j)=q.linear_residual/max(q.r0^2,realmin);
    a.face_residual(j)=q.face_residual;
    a.face_ratio(j)=q.face_residual/max(q.r0^2,realmin);
    a.fit_mask(j)=q.accepted && q.feasible && q.interior && q.stable && ...
        q.r0>0 && q.r1>0 && q.normal_ratio<=2000 && ...
        q.linear_residual<=max(0.2*q.r0^2,2e-14) && ...
        a.face_ratio(j)<=o.local_face_ratio_cap;
    fprintf('[Paper Fig5] seed%d sample%d r0=%.3g r1=%.3g used=%d\n', ...
        p.seed,j,q.r0,q.r1,a.fit_mask(j));
end
a.active_count=nnz(active); a.user_count=U;
end

function a=reference_audit(p,z,active,st)
A=st.G(:,active); scales=max(vecnorm(A),realmin);
singular=svd(A./scales,'econ'); rankA=nnz(singular>1e-11);
% An offline dense horizontal spectrum checks the SOSC used by the local test.
[Q,~]=qr(st.Q); Z=Q(:,size(st.Q,2)+1:end); HZ=zeros(size(Z));
for j=1:size(Z,2); HZ(:,j)=hvp(p,z,Z(:,j),st.lambda); end
Hh=Z'*HZ; Hh=(Hh+Hh')/2;
if isempty(Hh); eigenvalues=inf; else; eigenvalues=real(eig(Hh,'vector')); end
minimum=min(eigenvalues); scaleH=max([1;abs(eigenvalues)]);
rg=norm(st.g+st.G*st.lambda)/max(1,norm(st.g));
rc=max([0;st.c]); rl=norm(st.lambda.*st.c,inf)/max(1,norm(st.lambda,inf));
a=struct('licq',rankA==nnz(active),'positive_multipliers',all(st.lambda(active)>0), ...
    'minimum_active_multiplier',min([inf;st.lambda(active)]), ...
    'horizontal_lambda_min',minimum,'horizontal_positive_definite',minimum>1e-10*scaleH, ...
    'horizontal_spectral_scale',scaleH,'face_residual',norm(st.c(active)), ...
    'tangent_residual',norm(st.gt),'active_rank',rankA,'horizontal_dimension',size(Z,2), ...
    'full_nnls_kkt',max([rc,rg,rl]),'primal_residual',rc, ...
    'normalized_lagrangian_gradient',rg,'complementarity',rl);
end

function q=one_step(p,z,active,Delta)
st=state(p,z,active); H=@(v)st.project(hvp(p,z,st.project(v),st.lambda));
cg=struct('max_iter',240,'forcing_power',2,'forcing_constant',0.1, ...
    'absolute_floor',1e-14,'nc_tol',0,'flat_tol',0,'detect_negative_curvature',true);
[s,info]=trust_region_tcg(st.gt,H,Delta,cg); s=st.project(s);
% Preserve the nearby active level set with an inward cubic displacement.
target=st.c(active)-norm(s)^3;
z1=level_restore(p,z+s,active,target); next=state(p,z1,active);
pred=-(st.g'*s+0.5*s'*H(s)); ared=st.f-next.f;
feasible=max([0;next.c])<=1e-10;
normratio=norm(z1-z-s)/max(norm(s)^2,realmin);
threshold=[p.paperfig_working_thresholds(1)*ones(p.cfg.sys.U,1); ...
    p.paperfig_working_thresholds(2)*ones(p.cfg.sys.M,1)];
stable=isequal(st.c>=-threshold,active) && isequal(next.c>=-threshold,active);
accepted=feasible && pred>0 && ared>=0.05*pred && normratio<=2000;
q=struct('z0',z,'z1',z1,'r0',norm(st.gt),'r1',norm(next.gt), ...
    'accepted',accepted,'feasible',feasible,'interior',norm(s)<0.99*Delta, ...
    'stable',stable,'linear_residual',info.final_res_norm,'normal_ratio',normratio, ...
    'face_residual',norm(st.c(active)), ...
    'predicted',pred,'actual',ared,'ratio',ared/pred,'cg',info,'step_norm',norm(s));
end

function st=state(p,z,active)
N=p.cfg.sys.M*p.cfg.sys.Nt; U=p.cfg.sys.U;
W=real_vector_to_complex_matrix(z,N,U); [f,gc]=compute_crlb_objective_grad(p.cfg,p.smat,W);
g=complex_matrix_to_real_vector(gc); [c,G]=costar_sosp_constraints(p.cfg,p.channels,W);
A=G(:,active); lambda=zeros(size(c));
if ~isempty(A)
    scales=max(vecnorm(A),realmin); mu=lsqnonneg(A./scales,-g); lambda(active)=mu./scales.';
end
phase=zeros(numel(z),U);
for u=1:U
    V=zeros(N,U); V(:,u)=1i*W(:,u); phase(:,u)=complex_matrix_to_real_vector(V);
end
B=[A phase]; B=B./max(vecnorm(B),realmin);
Q=orth(B,1e-11); project=@(v)v-Q*(Q'*v);
st=struct('f',f,'g',g,'c',c,'G',G,'lambda',lambda,'Q',Q,'project',project,'gt',project(g));
end

function z=level_restore(p,z,active,target)
for k=1:20
    W=real_vector_to_complex_matrix(z,p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
    [c,G]=costar_sosp_constraints(p.cfg,p.channels,W); r=c(active)-target;
    if norm(r,inf)<=2e-13; break; end
    correction=lsqminnorm(G(:,active)',-r,1e-11); z=z+correction;
end
end

function v=hvp(p,z,d,lambda)
N=p.cfg.sys.M*p.cfg.sys.Nt; U=p.cfg.sys.U;
v=complex_matrix_to_real_vector(exact_lagrangian_hess_vec(p.cfg,p.channels,p.smat, ...
    real_vector_to_complex_matrix(z,N,U),real_vector_to_complex_matrix(d,N,U),lambda));
end
