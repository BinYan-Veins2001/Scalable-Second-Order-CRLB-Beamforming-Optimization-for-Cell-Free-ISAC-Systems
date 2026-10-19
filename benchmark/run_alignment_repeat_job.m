function result=run_alignment_repeat_job(job_file)
%RUN_ALIGNMENT_REPEAT_JOB Decode one frozen repetition job in a fresh process.
if endsWith(job_file,'.json')
    job=jsondecode(fileread(job_file));
    [folder,name]=fileparts(job_file);
    job_file=fullfile(folder,[name '.mat']);
    save(job_file,'job');
end
result=run_paired_job(job_file);
end
