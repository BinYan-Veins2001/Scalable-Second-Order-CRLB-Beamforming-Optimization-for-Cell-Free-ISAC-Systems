function archive = simtwc_compute_fig20(options)
%SIMTWC_COMPUTE_FIG20 Conditional Theorem 3 numerical validation.
if nargin<1; options=struct(); end
s=simtwc_setup(); o=simtwc_theorem3_options(options);
suffix=''; if ~isempty(o.output_tag); suffix=['_' o.output_tag]; end
archive=fullfile(s.data_dir,['fig20_data' suffix '.mat']);
signature=struct('version','theorem3_tangent_audit_v1','options',o, ...
    'profile',[s.smoke s.quick],'sources',source_signature(s.root));
data=struct('signature',signature,'options',o,'trials', ...
    {cell(numel(o.regimes),numel(o.seeds))}, ...
    'done',false(numel(o.regimes),numel(o.seeds)), ...
    'complete',false,'trajectory',[],'oracle',[], ...
    'notes',{{'Conditional exact-arithmetic theory; numerical evidence only.', ...
    'Dense spectra and cone sampling are offline diagnostics excluded from solver runtime.', ...
    'Finite restoration does not prove entry into the exact identified-set regime.'}});
if o.resume && exist(archive,'file')
    previous=load(archive);
    if isfield(previous,'signature') && isequaln(previous.signature,signature)
        data=previous;
    else
        fprintf('[Fig20] Settings/source changed; recomputing this experiment.\n');
    end
end
for ir=1:numel(o.regimes)
    for is=1:numel(o.seeds)
        if data.done(ir,is); continue; end
        spec=o.regimes(ir); seed=o.seeds(is);
        fprintf('[Fig20] %s, seed %d (%d/%d)\n',spec.name,seed,is,numel(o.seeds));
        trial=struct('seed',seed,'regime',spec.name,'status','not_run', ...
            'message','','initialization_runtime_s',NaN,'optimization_runtime_s',NaN, ...
            'entry_feasible',false,'solver_conditional_sosp',false, ...
            'numerical_joint_pass',false,'numerical_geometry_ok',false, ...
            'audit',[],'solution',[],'history',[]);
        try
            p=simtwc_make_problem(spec,seed);
            trial.initialization_runtime_s=p.initialization.runtime_s;
            trial.entry_feasible=p.initialization.success;
            if ~p.initialization.success
                trial.status='entry_failed'; trial.message=p.initialization.message;
            else
                term=p.cfg.costar.sosp;
                term.max_iter=o.terminal_max_iter;
                term.epsilon_KKT=o.epsilon_KKT; term.epsilon_H=o.epsilon_H;
                term.failure_probability=o.failure_probability;
                term.curvature_space='active_tangent';
                term.use_structured_psd=o.use_structured_psd;
                term.meo_record_history=true; term.meo_overlap_safety_factor=1;
                term.seed=seed+51000;
                [sol,hist]=solve_costar(p.cfg,p.channels,p.smat, ...
                    p.initialization.wcell,struct('terminal_options',term));
                trial.optimization_runtime_s=sol.runtime;
                p.cfg.costar.sosp=sol.terminal_options;
                audit=simtwc_theorem3_audit(p,sol.z,o,true);
                trial.audit=audit; trial.history=hist;
                if isfield(sol,'hist'); sol=rmfield(sol,'hist'); end
                trial.solution=sol; trial.status='completed'; trial.message=sol.message;
                trial.solver_conditional_sosp=sol.conditional_sosp_reached;
                trial.numerical_joint_pass=audit.primal<=sol.terminal_options.feas_tol && ...
                    audit.joint_ratio<=1;
                trial.numerical_geometry_ok=audit.licq && audit.positive_multipliers && ...
                    audit.face_residual<=sol.terminal_options.face_tol && ...
                    audit.nnls_ls_defect<=1e-8;
                fprintf('  KKT=%.3g, dense lambda_T=%.3g, conditional=%d, q=%.3g\n', ...
                    audit.kkt,audit.dense_lambda_min,sol.conditional_sosp_reached,audit.joint_ratio);
            end
        catch exception
            trial.status='error'; trial.message=exception.message;
            warning('Fig20:TrialFailed','%s, seed %d: %s',spec.name,seed,exception.message);
        end
        data.trials{ir,is}=trial; data.done(ir,is)=true;
        simtwc_save_atomic(archive,data);
    end
end

% A single declared realization is used for the trajectory, not a best-seed selection.
representative=data.trials{1,1};
if ~strcmp(representative.status,'completed')
    export_tables(data,s,suffix); simtwc_save_atomic(archive,data);
    error('Fig20:RepresentativeFailed', ...
        'The declared representative trial failed: %s. All trial records were saved.', ...
        representative.message);
end
p=simtwc_make_problem(o.regimes(1),o.seeds(1));
p.cfg.costar.sosp=representative.solution.terminal_options;
if isempty(data.trajectory)
    data.trajectory=trajectory_audit(p,representative.history.rtr,o);
    simtwc_save_atomic(archive,data);
end
if isempty(data.oracle)
    [endpoint,op]=simtwc_theorem3_audit(p,representative.solution.z,o,true);
    nt=size(op.Ht,1);
    oracle_opts=struct('epsilon_H',o.epsilon_H, ...
        'failure_probability',o.failure_probability,'norm_bound',norm(op.Ht,1), ...
        'effective_dim',nt,'max_iter',o.oracle_trace_max_iter, ...
        'record_history',true,'overlap_safety_factor',1,'seed',o.base_seed+170003);
    fprintf('[Fig20] Independent endpoint oracle, tangent dimension %d\n',nt);
    timer=tic;
    info=costar_sosp_meo(@(v) op.Ht*v,nt,@(v) v,oracle_opts);
    data.oracle=struct('trace',info,'dense_lambda_min',endpoint.dense_lambda_min, ...
        'diagnostic_runtime_s',toc(timer),'random_trials',[],'negative_case_iter',NaN);
    % Dense reduced actions accelerate repeated OFFLINE oracle audits only.
    data.oracle.random_trials=repeat_oracle(op.Ht,endpoint.dense_lambda_min,o,0);
    indices=find(data.trajectory.lambda_min<-o.epsilon_H);
    if ~isempty(indices)
        index=indices(1); k=data.trajectory.spectral_indices(index);
        [negative,nop]=simtwc_theorem3_audit(p, ...
            representative.history.rtr.z_snapshots(:,k),o,true);
        data.oracle.negative_case_iter=representative.history.rtr.iter(k);
        data.oracle.random_trials=[data.oracle.random_trials; ...
            repeat_oracle(nop.Ht,negative.dense_lambda_min,o,1)];
    end
    simtwc_save_atomic(archive,data);
end
data.complete=all(data.done(:)); data.completed_at=char(datetime('now'));
export_tables(data,s,suffix);
simtwc_save_atomic(archive,data);
fprintf('[Fig20] MAT archive: %s\n',archive);
end

function t=trajectory_audit(p,h,o)
n=numel(h.iter); t=struct('iter',h.iter,'kkt',nan(1,n), ...
    'primal',nan(1,n),'spectral_indices',[], ...
    'lambda_min',[],'terminal_start_iter',h.terminal_start_iter);
indices=unique([1:o.curvature_stride:n,n,h.terminal_start_iter+1, ...
    max(1,h.terminal_start_iter)]);
indices=indices(indices>=1 & indices<=n); t.spectral_indices=indices;
for k=1:n
    spectral=any(indices==k);
    a=simtwc_theorem3_audit(p,h.z_snapshots(:,k),o,spectral);
    t.kkt(k)=a.kkt; t.primal(k)=a.primal;
    if spectral
        t.lambda_min(end+1)=a.dense_lambda_min;
        fprintf('[Fig20] Spectrum at iteration %d/%d: %.4g\n',k-1,n-1,a.dense_lambda_min);
    end
end
end

function rows=repeat_oracle(Ht,truth,o,case_id)
rows=zeros(o.random_trials,11); n=size(Ht,1);
for j=1:o.random_trials
    opts=struct('epsilon_H',o.epsilon_H,'norm_bound',norm(Ht,1), ...
        'failure_probability',o.failure_probability,'max_iter',o.random_audit_budget, ...
        'effective_dim',n,'seed',o.base_seed+700001+case_id*10000+j, ...
        'overlap_safety_factor',1);
    info=costar_sosp_meo(@(v) Ht*v,n,@(v) v,opts);
    margin=1e-9*max(1,norm(Ht,1));
    valid=isfinite(info.eigenvalue_lower_bound);
    covered=valid && info.eigenvalue_lower_bound<=truth+margin;
    negative=strcmp(info.status,'negative_curvature');
    false_pass=info.conditional_psd && truth<-o.epsilon_H-margin;
    rows(j,:)=[case_id,j,truth,info.eigenvalue_lower_bound,info.rayleigh, ...
        info.iterations,info.hvp_calls,valid,covered,negative,false_pass];
end
end

function export_tables(data,s,suffix)
rows=cell(0,34);
for ir=1:size(data.trials,1)
    for is=1:size(data.trials,2)
        t=data.trials{ir,is}; if isempty(t); continue; end
        values=nan(1,22); scope=''; oracle_status='';
        if ~isempty(t.audit)
            a=t.audit; w=t.solution.work; wt=t.solution.work_terminal;
            scope=t.solution.certificate.scope;
            oracle_status=t.solution.certificate.oracle.status;
            values=[a.crlb,a.kkt,a.dense_lambda_min,a.joint_ratio, ...
                a.sigma_min_raw,a.lambda_min_active,a.face_residual,a.nnls_ls_defect, ...
                a.epsilon2_requested,a.epsilon2_residual,a.cone_sample_violation, ...
                a.cone_sample_count,a.hvp_relative_error, ...
                w.tcg_negative_curvature_detections,w.tcg_negative_curvature_accepted, ...
                wt.ritz_negative_curvature_detections,wt.ritz_escape_accepted,wt.meo_calls, ...
                a.dense_ambient_minimum,a.ambient_joint_pass, ...
                t.solution.work_fast.ritz_negative_curvature_detections,a.nnls_exitflag];
        end
        rows(end+1,:)=[{t.regime,t.seed,t.status,t.entry_feasible, ...
            t.solver_conditional_sosp,t.numerical_joint_pass,t.numerical_geometry_ok, ...
            t.initialization_runtime_s,t.optimization_runtime_s,t.message,scope,oracle_status}, ...
            num2cell(values)]; %#ok<AGROW>
    end
end
names={'Regime','Seed','Status','EntryFeasible','ConditionalStop','NumericalJointPass', ...
    'NumericalGeometryOK','InitializationRuntime_s','OptimizationRuntime_s', ...
    'TerminationMessage','CurvatureScope','OracleStatus', ...
    'CRLB','KKT','DenseLambdaT','JointRatio','SigmaMinRaw', ...
    'MinActiveMultiplier','FaceResidual','NNLSLSDefect','Epsilon2Requested', ...
    'Epsilon2Residual','SampleConeViolation','ConeSampleCount','HVPRelativeError', ...
    'TCGNCDetected','TCGNCAccepted','OracleNCDetected','OracleNCAccepted','OracleCalls', ...
    'DenseAmbientMinimum','AmbientJointPass','FastRitzNCDetected','NNLSExitFlag'};
writetable(cell2table(rows,'VariableNames',names), ...
    fullfile(s.figure_dir,['fig20_theorem3_trials' suffix '.csv']));
summary=cell(size(data.trials,1),10);
for ir=1:size(data.trials,1)
    trials=data.trials(ir,:); attempted=~cellfun(@isempty,trials); trials=trials(attempted);
    n=numel(trials); entries=0; stops=0; joints=0; geometry=0;
    runtime=[]; kkt=[]; curvature=[];
    for is=1:n
        t=trials{is}; entries=entries+t.entry_feasible;
        stops=stops+t.solver_conditional_sosp; joints=joints+t.numerical_joint_pass;
        geometry=geometry+t.numerical_geometry_ok;
        if ~isempty(t.audit)
            runtime(end+1)=t.optimization_runtime_s; kkt(end+1)=t.audit.kkt; %#ok<AGROW>
            curvature(end+1)=t.audit.dense_lambda_min; %#ok<AGROW>
        end
    end
    summary(ir,:)={data.options.regimes(ir).name,n,entries,stops,joints,geometry, ...
        100*stops/max(1,n),median(runtime),max([NaN kkt],[],'omitnan'), ...
        min([NaN curvature],[],'omitnan')};
end
writetable(cell2table(summary,'VariableNames',{'Regime','Trials','FeasibleEntries', ...
    'ConditionalStops','NumericalJointPasses','NumericalGeometryOK', ...
    'ConditionalStopRate_pct','MedianOptimizationRuntime_s','WorstKKT','MinDenseLambdaT'}), ...
    fullfile(s.figure_dir,['fig20_theorem3_summary' suffix '.csv']));
if ~isempty(data.oracle)
    names={'CaseID','Repeat','DenseLambdaT','LowerBound','Rayleigh','SearchDimension', ...
        'HVPCalls','FiniteBound','BoundCovered','NCDetected','FalsePass'};
    writetable(array2table(data.oracle.random_trials,'VariableNames',names), ...
        fullfile(s.figure_dir,['fig20_theorem3_oracle' suffix '.csv']));
    audit=data.oracle.random_trials; cases=unique(audit(:,1)); orows=zeros(numel(cases),8);
    for j=1:numel(cases)
        group=audit(audit(:,1)==cases(j),:); valid=group(:,8)>0;
        count=nnz(valid); failures=sum(group(valid,9)==0);
        upper=NaN;
        if count>0
            if failures==count; upper=1;
            else; upper=betaincinv(0.95,failures+1,count-failures); end
        end
        orows(j,:)=[cases(j),size(group,1),count,failures,upper, ...
            sum(group(:,10)),sum(group(:,11)),data.options.failure_probability];
    end
    writetable(array2table(orows,'VariableNames',{'CaseID','Trials','FiniteBounds', ...
        'ObservedEnclosureFailures','FailureProbability95Upper', ...
        'NCDetections','FalsePasses','RequestedPerCallDelta'}), ...
        fullfile(s.figure_dir,['fig20_theorem3_oracle_summary' suffix '.csv']));
end
end

function signature=source_signature(root)
folders={'optimization','config','baselines','channels','scenario','crlb','utils', ...
    fullfile('sim_twc','common'),fullfile('sim_twc','compute')};
digest=java.security.MessageDigest.getInstance('SHA-256');
for j=1:numel(folders)
    files=dir(fullfile(root,folders{j},'*.m'));
    for k=1:numel(files)
        path=fullfile(files(k).folder,files(k).name);
        digest.update(unicode2native([folders{j} '/' files(k).name fileread(path)],'UTF-8'));
    end
end
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
