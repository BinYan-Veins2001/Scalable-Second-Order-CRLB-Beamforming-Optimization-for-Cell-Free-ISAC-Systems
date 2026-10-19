function archive=paperfig_compute04(options)
%PAPERFIG_COMPUTE04 Dense independent spectra plus exact full-stop records.
if nargin<1; options=struct(); end
s=paperfig_setup(); o=paperfig_options(options); regimes=s.fig06.regimes;
if isfield(o,'regimes'); regimes=o.regimes; end
folder=fullfile(s.results_dir,'paper_fig04'); if ~isfolder(folder); mkdir(folder); end
archive=fullfile(folder,'data.mat'); signature=paperfig_signature(); audit_signature=paperfig_audit_signature(4);
d=struct('figure_id',4,'signature',signature,'audit_signature',audit_signature,'options',o,'regimes',regimes, ...
    'seeds',o.seeds,'trials',{cell(numel(regimes),numel(o.seeds))}, ...
    'done',false(numel(regimes),numel(o.seeds)),'complete',false);
if o.resume && isfile(archive)
    old=load(archive);
    if isfield(old,'signature') && strcmp(old.signature,signature) && ...
            isfield(old,'audit_signature') && strcmp(old.audit_signature,audit_signature) && isequaln(old.options,o); d=old; end
end
keys={'default','high_sinr','more_aps'};
for ir=1:numel(regimes)
    if ir>numel(keys); keys{ir}=sprintf('regime_%d',ir); end
    spec=rmfield(regimes(ir),'name');
    for t=1:numel(o.seeds)
        if d.done(ir,t); continue; end
        fprintf('[Paper Fig4] %s, seed %d\n',regimes(ir).name,o.seeds(t));
        [p,run]=paperfig_case(keys{ir},spec,o.seeds(t),'costar',o);
        audit=paperfig_dense_audit(p,run.solution,o);
        d.trials{ir,t}=struct('regime',regimes(ir).name,'seed',o.seeds(t), ...
            'status','completed','audit',audit,'solution',run.solution, ...
            'runtime_s',run.runtime_median_s,'history',run.history);
        d.done(ir,t)=true; simtwc_save_atomic(archive,d);
        fprintf('  KKT=%.3g eigT=%.3g q_stop=%.3g fullstop=%d epsilon2=%.3g\n', ...
            audit.kkt,audit.dense_lambda_min,audit.joint_ratio,run.solution.success,audit.epsilon2);
    end
end
d.complete=all(d.done(:)); d.completed_at=char(datetime('now')); simtwc_save_atomic(archive,d);
rows=cell(0,18);
for ir=1:numel(regimes)
    for t=1:numel(o.seeds)
        q=d.trials{ir,t}; a=q.audit;
        rows(end+1,:)={q.regime,q.seed,a.crlb,a.primal,a.kkt,a.dense_lambda_min, ...
            a.joint_ratio,a.epsilon2,a.joint_ratio_certificate,a.joint_ratio_full, ...
            a.solver_success,a.fullspace_certificate,a.curvature_lower_bound, ...
            a.sigma_min,a.minimum_active_multiplier,a.hvp_relative_error,q.runtime_s,q.solution.message}; %#ok<AGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'regime','seed','CRLB_m2','primal', ...
    'r_KKT','dense_lambda_min_T','q_stop_dense','epsilon2_guard','q_certificate', ...
    'q_full','solver_full_stop','fullspace_certificate','ell','sigma_min', ...
    'min_active_multiplier','hvp_error','runtime_s','message'}),fullfile(folder,'joint_stationarity_trials.csv'));
end
