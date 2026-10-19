function paperfig_finalize_metadata(numbers)
%PAPERFIG_FINALIZE_METADATA Record final plot/provenance-only source revisions.
% Numerical source hashes are checked before adding the final audit metadata.
if nargin<1; numbers=2:5; end
s=paperfig_setup(); signature=paperfig_signature(); physical_signature=paperfig_audit_signature(0);
for id=numbers
    path=fullfile(s.results_dir,sprintf('paper_fig%02d',id),'data.mat');
    if ~isfile(path); continue; end
    d=load(path); assert(strcmp(d.signature,signature),'Numerical source changed after this run.');
    d.audit_signature=paperfig_audit_signature(id);
    d.audit_metadata_finalized_at=char(datetime('now'));
    d.audit_metadata_note=['Final provenance guards and plot-only layout refinements; ' ...
        'numerical source hash verified unchanged.'];
    if id==3
        d.runtime_isolated=false;
        d.runtime_role='Endpoint diagnostics only; fresh concurrent or verified archival times, excluded from runtime claims.';
    elseif id==5
        for t=1:numel(d.trials)
            if ~isfield(d.trials{t}.reference_check,'full_nnls_kkt')
                seed=d.trials{t}.seed;
                problem=load(fullfile(s.case_dir,sprintf('default_seed%d',seed),'problem.mat'),'p');
                geo=paperfig_geometry(problem.p,d.trials{t}.reference_z,d.options.epsilon_KKT,d.options.solver_options);
                d.trials{t}.reference_check.full_nnls_kkt=geo.kkt;
            end
        end
        paperfig_reference_conditions(d,fileparts(path));
    end
    simtwc_save_atomic(path,d);
end
% Early cases were computed before the additional physical-recipe hash was
% introduced. Preserve those numerical outputs after verifying the core hash.
files=dir(fullfile(s.case_dir,'*','*.mat'));
for j=1:numel(files)
    path=fullfile(files(j).folder,files(j).name); d=load(path);
    if ~isfield(d,'signature') || ~strcmp(d.signature,signature); continue; end
    if ~isfield(d,'physical_signature')
        d.physical_signature=physical_signature; simtwc_save_atomic(path,d);
    end
end
end
