function result=paperfig_compute_efficiency(options)
%PAPERFIG_COMPUTE_EFFICIENCY Paired timings at returned COSTAR endpoints.
% Three fixed independent realizations per system size; each is represented
% by its returned iterate, whether or not the solver certified SOSP. Each time is
% the median of three interleaved warm calls. Curves show the median across
% realizations. Problem generation, optimization and audits are not timed.
if nargin<1; options=struct(); end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root),'-begin');
o=defaults(options);
signature=paperfig_signature();
benchmark_signature=paperfig_efficiency_signature();
thread_count=maxNumCompThreads;
assert(thread_count==1, ...
    'Use a fresh MATLAB -singleCompThread process for paired timing.');
folder=fullfile(root,'results','paper_fig06');
jobs=fullfile(folder,'jobs'); if ~exist(jobs,'dir'); mkdir(jobs); end
s=simtwc_config(); base=s.default;
points=struct('panel',{},'x',{},'M',{},'U',{});
for M=o.M_values; points(end+1)=struct('panel','AP','x',M,'M',M,'U',o.fixed_U); end %#ok<AGROW>
for U=o.U_values; points(end+1)=struct('panel','User','x',U,'M',o.fixed_M,'U',U); end %#ok<AGROW>
rows=struct([]); completed=0;
for j=1:numel(points)
    point=points(j);
    for seed=o.seeds
        spec=base; spec.M=point.M; spec.U=point.U; spec.Nt=o.Nt; spec.Nr=o.Nr;
        tag=sprintf('M%d_U%d_Nt%d_seed%d_R%d',spec.M,spec.U,spec.Nt,seed,o.repetitions);
        path=fullfile(jobs,[tag '.mat']);
        cached=false;
        if o.resume && exist(path,'file')
            a=load(path,'job'); job=a.job;
            cached=strcmp(job.schema,o.schema) && isfield(job,'signature') && ...
                strcmp(job.signature,signature) && isequal(job.solver_options,o.solver_options) && ...
                isequal(job.spec,spec) && job.summary.repetitions==o.repetitions && ...
                isfield(job,'benchmark_source_signature') && ...
                strcmp(job.benchmark_source_signature,benchmark_signature);
            if cached && ~isfield(job.summary,'snapshot_source')
                % A successful solver validates its final-point certificate.
                % The saved iteration must also equal the final loop counter.
                cached=job.summary.solver_success && ...
                    job.summary.snapshot_iter==job.solver_work.rtr_iterations;
                if cached
                    job.summary.snapshot_source='returned_endpoint';
                    job.summary.actual_terminal_spsc=true;
                    job.summary.terminal_iter=job.solver_work.rtr_iterations;
                    job.policy_migration=['Retained measured calls: saved successful solver exit ' ...
                        'and snapshot_iter==saved solver_work.rtr_iterations establish that ' ...
                        'the saved last SPSC call was at the returned endpoint.'];
                    save(path,'job','-v7.3');
                end
            end
        end
        if ~cached
            job=run_job(spec,seed,o);
            job.signature=signature; job.solver_options=o.solver_options;
            job.benchmark_source_signature=benchmark_signature;
            save(path,'job','-v7.3');
        end
        row=job.summary; row.panel=point.panel; row.x=point.x;
        if isempty(rows); rows=row; else; rows(end+1)=row; end %#ok<AGROW>
        completed=completed+1;
        fprintf('Fig6 %d/%d %s=%d seed%d: SPSC %.6gs, dense %.6gs, speedup %.2fx, %s, solver_pass=%d\n', ...
            completed,numel(points)*numel(o.seeds),point.panel,point.x,seed, ...
            row.spsc_median_s,row.dense_median_s,row.speedup,row.spsc_status,row.solver_success);
    end
end
result=struct('schema',o.schema,'options',o,'rows',rows, ...
    'created',datestr(now,30),'matlab_version',version,'computer',computer, ...
    'max_computational_threads',thread_count,'methodology',methodology(), ...
    'signature',signature,'point_trial_rows',numel(rows), ...
    'benchmark_source_signature',benchmark_signature, ...
    'unique_trajectory_jobs',size(unique([[rows.M].',[rows.U].',[rows.seed].'],'rows'),1));
save(fullfile(folder,'data.mat'),'result','-v7.3');
writetable(struct2table(rows),fullfile(folder,'paired_trials.csv'));
write_summary(result,folder);
end

function job=run_job(spec,seed,o)
source=prepare_endpoint(spec,seed,o);
p=source.problem; sol=source.solution;
z=sol.z; N=spec.M*spec.Nt; U=spec.U;
W=real_vector_to_complex_matrix(z,N,U);
[~,gc]=compute_crlb_objective_grad(p.cfg,p.smat,W);
[c,G]=costar_sosp_constraints(p.cfg,p.channels,W);
geo=costar_paper_geometry(complex_matrix_to_real_vector(gc),c,G,U,sol.endpoint_geometry_nu,sol.options);
assert(geo.success,'Cannot reconstruct the final working geometry.');
assert(isequal(geo.active,sol.active) && ...
    norm(geo.lambda-sol.lambda)<=1e-10*max(1,norm(sol.lambda)), ...
    'Reconstructed endpoint multipliers do not match the solver.');
lambda=geo.lambda; Q=geo.Q; epsilon_H=sol.options.epsilon_H;
structured=@()costar_sosp_structured_psd(p.cfg,p.channels,p.smat,W,lambda,epsilon_H,Q);
warm=structured();
assert(isfinite(warm.shift) && (warm.conditional_psd || strcmp(warm.status,'negative_curvature')), ...
    'SPSC did not resolve the saved call.');
dense=@()paperfig_efficiency_dense_evd(p,W,lambda,Q,warm.shift);
reference=dense();
audit=audit_pair(p,W,lambda,Q,warm,reference,epsilon_H,seed);
spsc_times=zeros(o.repetitions,1); dense_times=spsc_times;
for rep=1:o.repetitions
    % Alternate order to distribute warm-cache and slow-drift effects.
    if mod(rep,2)==1
        t=tic; si=structured(); spsc_times(rep)=toc(t);
        t=tic; di=dense(); dense_times(rep)=toc(t);
    else
        t=tic; di=dense(); dense_times(rep)=toc(t);
        t=tic; si=structured(); spsc_times(rep)=toc(t);
    end
    assert(si.tangent_negative_inertia==di.negative_inertia);
    assert(si.conditional_psd==di.certified && si.shift==warm.shift);
end
summary=struct('M',spec.M,'U',U,'Nt',spec.Nt,'Nr',spec.Nr,'seed',seed, ...
    'dimension',2*N*U,'normal_rank',size(Q,2),'snapshot_iter',sol.terminal_iterations, ...
    'actual_spsc_calls',source.actual_spsc_calls,'solver_success',sol.success, ...
    'snapshot_kkt',geo.kkt,'solver_runtime_s',sol.runtime, ...
    'epsilon_H',epsilon_H,'sigma',warm.shift, ...
    'spsc_status',warm.status,'spsc_mode',warm.mode, ...
    'ambient_negative_inertia',warm.ambient_negative_inertia, ...
    'tangent_negative_inertia',warm.tangent_negative_inertia, ...
    'lambda_min_tangent',reference.lambda_min,'spsc_lower_bound',warm.eigenvalue_lower_bound, ...
    'spsc_rayleigh',warm.rayleigh,'hvp_relative_error',audit.hvp_relative_error, ...
    'inertia_match',audit.inertia_match,'direction_pass',audit.direction_pass, ...
    'spsc_median_s',median(spsc_times),'dense_median_s',median(dense_times), ...
    'speedup',median(dense_times)/median(spsc_times),'repetitions',o.repetitions, ...
    'snapshot_source','returned_endpoint', ...
    'actual_terminal_spsc',source.actual_terminal_spsc, ...
    'terminal_iter',sol.terminal_iterations);
snapshot=struct('z',z,'lambda',lambda,'Q',Q,'active',geo.active, ...
    'nu',sol.endpoint_geometry_nu,'iteration',sol.terminal_iterations, ...
    'original_status',source.original_status);
job=struct('schema',o.schema,'spec',spec,'problem',p,'snapshot',snapshot, ...
    'summary',summary,'spsc',warm,'dense',reference,'audit',audit, ...
    'spsc_times_s',spsc_times,'dense_times_s',dense_times, ...
    'solver_message',sol.message,'solver_work',sol.work);
end

function source=prepare_endpoint(spec,seed,o)
% Preserve costly optimizer outputs before curvature timing/audits begin.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
folder=fullfile(root,'results','paper_fig06','endpoints');
if ~isfolder(folder); mkdir(folder); end
path=fullfile(folder,sprintf('M%d_U%d_Nt%d_seed%d.mat',spec.M,spec.U,spec.Nt,seed));
signature=paperfig_signature();
benchmark_signature=paperfig_efficiency_signature();
if o.resume && isfile(path)
    a=load(path,'source'); source=a.source;
    if strcmp(source.signature,signature) && isequal(source.spec,spec) && ...
            isequal(source.solver_options,o.solver_options) && ...
            isfield(source,'benchmark_source_signature') && ...
            strcmp(source.benchmark_source_signature,benchmark_signature); return; end
end
p=simtwc_make_problem(spec,seed);
assert(p.initialization.success,'Initialization failed: %s',p.initialization.message);
[sol,hist]=solve_costar(p.cfg,p.channels,p.smat,p.initialization.wcell,o.solver_options);
trace=hist.rtr;
assert(isequal(sol.z,trace.z_snapshots(:,end)) && trace.iter(end)==sol.terminal_iterations);
source=struct('problem',p,'solution',rmfield(sol,'hist'), ...
    'actual_spsc_calls',nnz(~strcmp(trace.meo_status,'not_called')), ...
    'actual_terminal_spsc',~strcmp(trace.meo_status{end},'not_called'), ...
    'original_status',trace.meo_status{end},'signature',signature, ...
    'spec',spec,'solver_options',o.solver_options, ...
    'benchmark_source_signature',benchmark_signature);
save(path,'source','-v7.3');
end

function audit=audit_pair(p,W,lambda,Q,spsc,dense,epsilon_H,seed)
H=paperfig_efficiency_dense_hessian(p,W,lambda); d=size(H,1);
old=rng; cleanup=onCleanup(@()rng(old)); rng(seed+900000,'twister'); %#ok<NASGU>
relative=0;
for j=1:3
    v=randn(d,1); v=v/norm(v);
    Hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec(p.cfg,p.channels,p.smat,W, ...
        real_vector_to_complex_matrix(v,size(W,1),size(W,2)),lambda));
    relative=max(relative,norm(H*v-Hv)/max(1,norm(Hv)));
end
assert(relative<5e-10,'Independent dense Hessian and exact HVP disagree.');
inertia_match=spsc.tangent_negative_inertia==dense.negative_inertia;
assert(inertia_match && spsc.conditional_psd==dense.certified,'Tangent inertia mismatch.');
direction_pass=true; tangent_error=NaN; rayleigh=NaN; unit_error=NaN;
if ~isempty(spsc.negative_direction)
    v=spsc.negative_direction; unit_error=abs(norm(v)-1);
    tangent_error=norm(Q'*v); rayleigh=v'*H*v;
    direction_pass=unit_error<1e-9 && tangent_error<1e-8 && rayleigh<-epsilon_H/2 ...
        && abs(rayleigh-spsc.rayleigh)<1e-9*max(1,abs(rayleigh));
    assert(direction_pass,'Recovered SPSC direction failed independent checks.');
end
if spsc.conditional_psd
    assert(dense.lambda_min>=spsc.eigenvalue_lower_bound-1e-9*max(1,norm(H,1)));
end
audit=struct('hvp_relative_error',relative,'inertia_match',inertia_match, ...
    'direction_pass',direction_pass,'tangent_error',tangent_error, ...
    'unit_error',unit_error,'rayleigh',rayleigh,'vectors_audited',3);
end

function o=defaults(options)
o=struct('schema','structured_efficiency_actual_spsc_v1', ...
    'M_values',[8 12 16 24 32],'U_values',[4 8 12 16], ...
    'fixed_U',8,'fixed_M',16,'Nt',4,'Nr',4,'seeds',2026:2028, ...
    'repetitions',3,'resume',true,'solver_options',struct());
names=fieldnames(options);
for j=1:numel(names)
    if isfield(o,names{j}); o.(names{j})=options.(names{j});
    elseif ~ismember(names{j},{'compute','plot'}); error('Unknown efficiency option: %s',names{j}); end
end
assert(numel(unique(o.seeds))==numel(o.seeds) && o.repetitions>=1);
end

function text=methodology()
text=['Returned endpoint of each fresh COSTAR trajectory, with the same rule for successful and ' ...
    'iteration-budget exits; solver_success and actual_terminal_spsc are recorded for every case. ' ...
    'This benchmark evaluates curvature computation at fixed feasible iterates, independently of whether ' ...
    'the optimization has reached its stationarity tolerances. Both methods receive the same ' ...
    'W, nonnegative NNLS lambda, orthonormal working-normal Q, and scalar sigma selected by SPSC. ' ...
    'SPSC recomputes all block factorizations and its full-space/tangent tests, including direction recovery ' ...
    'when needed. Dense tangent EVD independently assembles the analytic ambient Hessian, completes Q by ' ...
    'full QR, projects to its orthogonal complement, and computes the complete eigensystem. Both methods ' ...
    'are warmed once, then measured in interleaved order in a single-computational-thread MATLAB process. ' ...
    'Per-realization median of repeated calls; curve median across declared seeds. Problem construction, ' ...
    'COSTAR optimization, geometry/NNLS construction, independent HVP audits and all graphics are excluded. ' ...
    'The saved lower bound is not an eigenvalue estimate. No older experiment timings are imported.'];
end

function write_summary(result,folder)
rows=result.rows; data=struct([]); panels={'AP','User'};
for p=1:2
    use=strcmp({rows.panel},panels{p}); xx=unique([rows(use).x]);
    for x=xx
        idx=use & [rows.x]==x; r=rows(idx);
        a=struct('panel',panels{p},'x',x,'M',r(1).M,'U',r(1).U, ...
            'trials',numel(r),'spsc_median_s',median([r.spsc_median_s]), ...
            'dense_median_s',median([r.dense_median_s]), ...
            'paired_speedup_median',median([r.speedup]), ...
            'ambient_certificates',nnz(strcmp({r.spsc_mode},'ambient')), ...
            'tangent_certificates',nnz(strcmp({r.spsc_mode},'tangent')), ...
            'negative_directions',nnz(strcmp({r.spsc_mode},'negative_curvature')), ...
            'all_audits_pass',all([r.inertia_match]) && all([r.direction_pass]));
        if isempty(data); data=a; else; data(end+1)=a; end %#ok<AGROW>
    end
end
writetable(struct2table(data),fullfile(folder,'curve_summary.csv'));
fid=fopen(fullfile(folder,'TIMING_METHOD.md'),'w','n','UTF-8'); cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Computational efficiency of structured curvature evaluation\n\n%s\n\n',result.methodology);
fprintf(fid,'MATLAB: %s. Computational threads: %d.\n\n',result.matlab_version,result.max_computational_threads);
fprintf(fid,'Seeds: %s. Warm timed repeats: %d.\n\n',mat2str(result.options.seeds),result.options.repetitions);
fprintf(fid,'Point-trial rows: %d. Unique trajectory jobs: %d (the common M/U intersection is reused across panels).\n\n', ...
    result.point_trial_rows,result.unique_trajectory_jobs);
fprintf(fid,['Each trace contributes its returned endpoint. Kernel branch counts appear in ' ...
    '`curve_summary.csv`; all checkpoints and raw repeated timings are saved under `jobs/`. ' ...
    'An ambient certificate bypasses tangent solves by design; its timing is the actual SPSC workflow.\n']);
end
