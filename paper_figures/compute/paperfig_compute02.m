function archive=paperfig_compute02(options)
%PAPERFIG_COMPUTE02 Updated COSTAR CRLB, complete KKT, and horizontal history.
if nargin<1; options=struct(); end
s=paperfig_setup(); o=paperfig_options(options); seed=o.seeds(1);
folder=fullfile(s.results_dir,'paper_fig02'); if ~isfolder(folder); mkdir(folder); end
archive=fullfile(folder,'data.mat'); signature=paperfig_signature(); audit_signature=paperfig_audit_signature(2);
[p,run]=paperfig_case('default',s.default,seed,'costar',o);
h=run.history.rtr; n=numel(h.iter);
d=struct('figure_id',2,'signature',signature,'audit_signature',audit_signature,'options',o,'seed',seed, ...
    'run',run,'iter',h.iter,'crlb',h.crlb,'kkt',h.kkt, ...
    'tangent_grad_norm',h.tangent_grad_norm,'horizontal',nan(1,n), ...
    'done',false(1,n),'complete',false);
if o.resume && isfile(archive)
    old=load(archive);
    if isfield(old,'signature') && strcmp(old.signature,signature) && ...
            isfield(old,'audit_signature') && strcmp(old.audit_signature,audit_signature) && ...
            isequaln(old.options,o) && isequaln(old.iter,h.iter); d=old; end
end
for k=find(~d.done)
    [geo,~,~,~,~]=paperfig_geometry(p,h.z_snapshots(:,k),h.nu(k),run.solution.options);
    W=real_vector_to_complex_matrix(h.z_snapshots(:,k),p.cfg.sys.M*p.cfg.sys.Nt,p.cfg.sys.U);
    a=paperfig_horizontal(p,W,o,geo);
    d.horizontal(k)=a.lambda_min_horizontal; d.done(k)=true;
    if mod(k,10)==0 || k==n
        fprintf('[Paper Fig2] Horizontal audit %d/%d\n',k,n); simtwc_save_atomic(archive,d);
    end
end
d.complete=all(d.done); d.completed_at=char(datetime('now'));
simtwc_save_atomic(archive,d);
writetable(table(d.iter(:),d.crlb(:),d.kkt(:),d.tangent_grad_norm(:),d.horizontal(:), ...
    'VariableNames',{'iteration','CRLB_m2','r_KKT','tangent_gradient','lambda_min_H_estimate'}), ...
    fullfile(folder,'trajectory.csv'));
end
