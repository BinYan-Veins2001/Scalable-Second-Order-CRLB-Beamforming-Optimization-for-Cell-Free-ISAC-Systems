function archive=paperfig_compute03(options)
%PAPERFIG_COMPUTE03 Common finite-Lanczos audit of seven fresh solver endpoints.
if nargin<1; options=struct(); end
s=paperfig_setup(); o=paperfig_options(options);
folder=fullfile(s.results_dir,'paper_fig03'); if ~isfolder(folder); mkdir(folder); end
archive=fullfile(folder,'data.mat'); signature=paperfig_signature(); audit_signature=paperfig_audit_signature(3); methods=s.methods;
if isfield(o,'methods'); methods=o.methods; end
d=struct('figure_id',3,'signature',signature,'audit_signature',audit_signature,'options',o, ...
    'methods',{methods},'seeds',o.seeds,'runs',{cell(numel(methods),numel(o.seeds))}, ...
    'audits',{cell(numel(methods),numel(o.seeds))},'done',false(numel(methods),numel(o.seeds)), ...
    'complete',false,'runtime_isolated',false, ...
    'runtime_role','Endpoint diagnostics only; concurrent runs, excluded from runtime claims.');
if o.resume && isfile(archive)
    old=load(archive);
    if isfield(old,'signature') && strcmp(old.signature,signature) && ...
            isfield(old,'audit_signature') && strcmp(old.audit_signature,audit_signature) && isequaln(old.options,o); d=old; end
end
for t=1:numel(o.seeds)
    for a=1:numel(methods)
        if d.done(a,t); continue; end
        [p,run]=paperfig_case('default',s.default,o.seeds(t),methods{a},o);
        d.runs{a,t}=run;
        if strcmp(run.status,'ok')
            d.audits{a,t}=paperfig_horizontal(p,run.Wmat,o);
        else
            d.audits{a,t}=struct('lambda_min_horizontal',NaN,'kkt',NaN,'lanczos_iter',0);
        end
        d.done(a,t)=true; simtwc_save_atomic(archive,d);
    end
end
d.complete=all(d.done(:)); d.completed_at=char(datetime('now'));
simtwc_save_atomic(archive,d); rows=cell(0,9);
for a=1:numel(methods)
    for t=1:numel(o.seeds)
        r=d.runs{a,t}; q=d.audits{a,t};
        rows(end+1,:)={methods{a},o.seeds(t),r.status,r.feasible, ...
            r.metrics.crlb,q.kkt,q.lambda_min_horizontal,q.lanczos_iter,r.runtime_median_s}; %#ok<AGROW>
    end
end
writetable(cell2table(rows,'VariableNames',{'method','seed','status','feasible', ...
    'CRLB_m2','r_KKT','lambda_min_H_estimate','Lanczos_iterations','runtime_s'}), ...
    fullfile(folder,'endpoint_trials.csv'));
end
