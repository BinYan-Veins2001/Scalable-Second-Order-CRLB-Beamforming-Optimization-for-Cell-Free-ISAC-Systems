function value=paperfig_audit_signature(id)
%PAPERFIG_AUDIT_SIGNATURE Curvature, plotting, and physical recipe provenance.
% Kept separate from the core numerical hash used by timed figure 6 jobs.
s=paperfig_setup(); paths={fullfile('sim_twc','common','simtwc_config.m'), ...
    fullfile('sim_twc','common','simtwc_make_problem.m')};
if id>0
    names={'setup','options','case','geometry','horizontal','dense_audit','signature','audit_signature','export','manifest'};
    for j=1:numel(names)
        paths{end+1}=fullfile('paper_figures','common',['paperfig_' names{j} '.m']); %#ok<AGROW>
    end
    files=dir(fullfile(s.root,'sim_twc','common','*.m'));
    for j=1:numel(files); paths{end+1}=fullfile('sim_twc','common',files(j).name); end %#ok<AGROW>
    paths{end+1}=fullfile('paper_figures','compute',sprintf('paperfig_compute%02d.m',id));
    paths{end+1}=fullfile('paper_figures','plot',sprintf('paperfig_plot%02d.m',id));
    if id==3
        paths{end+1}=fullfile('paper_figures','common','paperfig_verify_baseline_archive.m');
        paths{end+1}=fullfile('paper_figures','common','paperfig_import_verified_baselines.m');
        paths{end+1}=fullfile('results','paper_fig03','baseline_archive','source_identity.json');
    elseif id==5
        paths{end+1}=fullfile('paper_figures','common','paperfig_reference_conditions.m');
    end
end
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(paperfig_signature(),'UTF-8'));
for j=1:numel(paths)
    path=fullfile(s.root,paths{j});
    if isfile(path); digest.update(unicode2native([paths{j} fileread(path)],'UTF-8')); end
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
