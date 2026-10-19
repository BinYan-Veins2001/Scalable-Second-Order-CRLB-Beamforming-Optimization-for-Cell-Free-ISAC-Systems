function jobs=prepare_alignment_validation()
%PREPARE_ALIGNMENT_VALIDATION Three frozen, default-size paired inputs.
root=fileparts(fileparts(mfilename('fullpath')));
old_root=fullfile(fileparts(root),'COSTAR_TWC_SOSP_SharedZF_20261003');
assert(isfile(fullfile(old_root,'optimization','solve_costar.m')));
vroot=fullfile(root,'benchmark','validation');
for name={'problems','jobs','raw','logs'}
    if ~isfolder(fullfile(vroot,name{1})); mkdir(fullfile(vroot,name{1})); end
end
jobs=cell(6,1); count=0;
for seed=2026:2028
    problem_file=fullfile(vroot,'problems',sprintf('default_seed%d.mat',seed));
    assert(isfile(problem_file),'Frozen physical input is missing.');
    if mod(seed,2)==0; labels={'new','old'}; else; labels={'old','new'}; end
    for j=1:2
        label=labels{j}; project_root=root;
        if strcmp(label,'old'); project_root=old_root; end
        stem=sprintf('%s_default_seed%d',label,seed);
        job=struct('project_root',project_root,'benchmark_root',fullfile(root,'benchmark'), ...
            'problem_file',problem_file,'label',label,'case_name','default', ...
            'seed',seed,'compute_threads',1,'solver_options',struct(), ...
            'dense_max_dim',2048,'output_file',fullfile(vroot,'raw',[stem '.mat']));
        job_file=fullfile(vroot,'jobs',[stem '.mat']);
        save(job_file,'job'); count=count+1; jobs{count}=job_file;
    end
end
save(fullfile(vroot,'jobs.mat'),'jobs');
fid=fopen(fullfile(vroot,'jobs.txt'),'w','n','UTF-8');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jobs{:});
fprintf('Prepared %d jobs; no optimization was run.\n',numel(jobs));
end
