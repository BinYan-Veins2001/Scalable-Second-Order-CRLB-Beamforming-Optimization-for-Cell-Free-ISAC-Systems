function report=paperfig_import_verified_baselines()
%PAPERFIG_IMPORT_VERIFIED_BASELINES Preserve unchanged endpoints with provenance.
% This explicit one-time import never imports COSTAR or old curvature audits.
s=paperfig_setup(); folder=fullfile(s.results_dir,'paper_fig03','baseline_archive');
identity=jsondecode(fileread(fullfile(folder,'source_identity.json')));
assert(identity.all_baseline_sources_match && identity.baseline_physical_config_prefix_matches);
report=readtable(fullfile(s.results_dir,'paper_fig03','archive_reuse_verification.csv'));
assert(height(report)==60 && all(report.eligible_initial_model));
fid=fopen(identity.archive,'rb'); bytes=fread(fid,inf,'*uint8'); fclose(fid);
digest=java.security.MessageDigest.getInstance('SHA-256'); digest.update(bytes);
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
assert(strcmp(hash,identity.archive_sha256),'Historical baseline archive changed.');
old=load(identity.archive); indices=find(~strcmp(old.methods,'costar'));
methods=old.methods(indices); seeds=old.seeds; runs=old.runs(indices,:);
source=struct('path',identity.archive,'sha256',hash,'spec',old.spec, ...
    'original_signature',old.signature,'baseline_sync_source',old.baseline_sync_source, ...
    'baseline_sync_gamma_dB',old.baseline_sync_gamma_dB);
save(fullfile(folder,'raw_baselines.mat'),'methods','seeds','runs','source','-v7.3');
signature=paperfig_signature(); physical_signature=paperfig_audit_signature(0);
solver_options=struct(); provenance_rows=cell(0,5);
for t=1:numel(s.seeds)
    seed=s.seeds(t); jt=find(old.seeds==seed,1);
    problem_spec=s.default; spec=problem_spec; p=simtwc_make_problem(spec,seed);
    casefolder=fullfile(s.case_dir,sprintf('default_seed%d',seed));
    if ~isfolder(casefolder); mkdir(casefolder); end
    save(fullfile(casefolder,'problem.mat'),'p','spec','signature','physical_signature','-v7.3');
    for a=indices
        method=old.methods{a}; path=fullfile(casefolder,[method '.mat']);
        fresh=false;
        if isfile(path)
            current=load(path);
            fresh=isfield(current,'signature') && strcmp(current.signature,signature) && ...
                isfield(current,'problem_spec') && isequaln(current.problem_spec,problem_spec) && ...
                isfield(current,'run') && strcmp(current.run.status,'ok');
        end
        if fresh
            run=current.run; origin='fresh current-folder solve';
        else
            run=old.runs{a,jt}; origin='verified unchanged baseline archive';
        end
        run.provenance=struct('origin',origin,'archive_source',identity.archive, ...
            'archive_sha256',hash,'initial_model_verified',true, ...
            'runtime_role','Diagnostic endpoint record only; not used for runtime claims.');
        if ~fresh; run.provenance.curvature_reused=false; end
        save(path,'run','signature','solver_options','problem_spec','seed','physical_signature','-v7.3');
        provenance_rows(end+1,:)={method,seed,origin,hash,false}; %#ok<AGROW>
    end
end
writetable(cell2table(provenance_rows,'VariableNames', ...
    {'method','seed','endpoint_origin','archive_sha256','old_curvature_reused'}), ...
    fullfile(folder,'endpoint_provenance.csv'));
fprintf('[Paper Fig3] Verified archival baselines imported; new COSTAR endpoints and fresh baseline solves retained.\n');
end
