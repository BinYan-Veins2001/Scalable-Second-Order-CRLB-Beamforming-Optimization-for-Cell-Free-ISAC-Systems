function report=paperfig_efficiency_finalize_provenance()
%PAPERFIG_EFFICIENCY_FINALIZE_PROVENANCE Backfill first-run provenance only.
% Run once after reviewing the first recorded experiment. No optimization,
% Hessian assembly, eigensolve or timed calls are performed here.
% The first run preceded the figure-specific fingerprint addition; its raw
% measurements and all independent audit results are retained exactly.
s=paperfig_setup(); folder=fullfile(s.results_dir,'paper_fig06');
path=fullfile(folder,'data.mat'); a=load(path,'result'); result=a.result;
numerical=paperfig_signature();
assert(strcmp(result.signature,numerical),'Numerical optimizer sources changed.');
assert(all([result.rows.inertia_match]) && all([result.rows.direction_pass]) && ...
    max([result.rows.hvp_relative_error])<5e-10,'First-run audits did not pass.');
[fingerprint,manifest]=paperfig_efficiency_signature();
note=['Initial experiment provenance backfill: timed dense-Hessian/EVD and SPSC kernel ' ...
    'implementations and physical problem-generation files were verified unchanged. ' ...
    'The compute wrapper was subsequently amended to use all returned endpoints, ' ...
    'preserve unsuccessful optimizer exits, and validate complete cache fingerprints. ' ...
    'Older successful samples were retained only after the saved final-iteration proof. ' ...
    'No timed calls or raw timing arrays were replaced during this provenance step.'];
jobs=dir(fullfile(folder,'jobs','*.mat'));
for j=1:numel(jobs)
    jobpath=fullfile(jobs(j).folder,jobs(j).name); b=load(jobpath,'job'); job=b.job;
    assert(strcmp(job.signature,numerical) && job.audit.inertia_match && job.audit.direction_pass);
    assert(strcmp(job.summary.snapshot_source,'returned_endpoint'));
    before_s=job.spsc_times_s; before_d=job.dense_times_s;
    job.benchmark_source_signature=fingerprint; job.source_manifest=manifest;
    job.provenance_note=note; save(jobpath,'job','-v7.3');
    c=load(jobpath,'job');
    assert(isequal(before_s,c.job.spsc_times_s) && isequal(before_d,c.job.dense_times_s));
end
endpoints=dir(fullfile(folder,'endpoints','*.mat'));
for j=1:numel(endpoints)
    pointpath=fullfile(endpoints(j).folder,endpoints(j).name);
    b=load(pointpath,'source'); source=b.source;
    assert(strcmp(source.signature,numerical));
    source.benchmark_source_signature=fingerprint;
    source.provenance_note=note; save(pointpath,'source','-v7.3');
end
result.benchmark_source_signature=fingerprint; result.source_manifest=manifest;
result.provenance_note=note;
save(path,'result','-v7.3');
[~,unique_indices]=unique([[result.rows.M].',[result.rows.U].',[result.rows.seed].'],'rows');
unique_rows=result.rows(unique_indices);
report=struct('passed',true,'jobs',numel(jobs),'endpoints',numel(endpoints), ...
    'benchmark_source_signature',fingerprint,'all_raw_timing_arrays_preserved',true, ...
    'optimizer_certified_endpoints',nnz([unique_rows.solver_success]), ...
    'optimizer_budget_exits',nnz(~[unique_rows.solver_success]), ...
    'ambient_kernel_certificates',nnz(strcmp({unique_rows.spsc_mode},'ambient')), ...
    'tangent_kernel_certificates',nnz(strcmp({unique_rows.spsc_mode},'tangent')), ...
    'kernel_negative_directions',nnz(strcmp({unique_rows.spsc_mode},'negative_curvature')));
fid=fopen(fullfile(folder,'source_manifest.json'),'w','n','UTF-8');
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(struct('report',report,'manifest',manifest,'note',note)));
md=fopen(fullfile(folder,'PROVENANCE.md'),'w','n','UTF-8');
mdclean=onCleanup(@()fclose(md)); %#ok<NASGU>
fprintf(md,'# Fig. 6 provenance and branch coverage\n\n%s\n\n',note);
fprintf(md,'Unique returned endpoints: %d; optimizer-certified endpoints: %d; iteration-budget exits: %d.\n\n', ...
    numel(unique_rows),report.optimizer_certified_endpoints,report.optimizer_budget_exits);
fprintf(md,['Timed SPSC branches: %d ambient certificates, %d tangent certificates, %d negative-direction recoveries. ' ...
    'The plotted curves therefore measure the ambient certification path of SPSC. ' ...
    'Tangent certification and negative-direction recovery are separately verified by ' ...
    '`tests/test_paperfig_efficiency.m`; their runtime is not represented by these curves.\n\n'], ...
    report.ambient_kernel_certificates,report.tangent_kernel_certificates,report.kernel_negative_directions);
fprintf(md,'Benchmark source SHA-256: `%s`. Per-file hashes are saved in `source_manifest.json`.\n',fingerprint);
disp(report);
end
