function report=summarize_alignment_validation()
%SUMMARIZE_ALIGNMENT_VALIDATION Compare performance and independent audits.
root=fileparts(fileparts(mfilename('fullpath')));
vroot=fullfile(root,'benchmark','validation');
loaded=load(fullfile(vroot,'jobs.mat'),'jobs'); rows=cell(numel(loaded.jobs),1);
for j=1:numel(loaded.jobs)
    input=load(loaded.jobs{j},'job');
    data=load(input.job.output_file,'result','sol');
    rows{j}=data.result.row;
    assert(rows{j}.executed && rows{j}.joint_pass, ...
        'Common endpoint audit failed for %s.',input.job.output_file);
    if strcmp(rows{j}.label,'new')
        assert(data.sol.success && data.sol.certificate.oracle_is_at_final_point);
        assert(rows{j}.paper_joint_pass && rows{j}.paper_primal<=1e-10);
        assert(data.sol.critical_cone_tolerance<=data.sol.options.epsilon_H);
        assert(data.sol.kkt_residual<=data.sol.options.epsilon_KKT);
        assert(data.sol.work.fd_hvp_calls==0 && data.sol.work.meo_calls==0);
    end
end
T=struct2table(vertcat(rows{:}));
writetable(T,fullfile(vroot,'trial_results.csv'));
old=T(strcmp(T.label,'old'),:); new=T(strcmp(T.label,'new'),:);
old=sortrows(old,'seed'); new=sortrows(new,'seed');
assert(isequal(old.seed,new.seed));
relative_crlb=abs(new.crlb-old.crlb)./old.crlb;
report=struct('passed',true,'paired_trials',height(old), ...
    'default_size',[16,4,4,8],'threads',1, ...
    'maximum_relative_crlb_difference',max(relative_crlb), ...
    'mean_old_runtime_s',mean(old.optimization_s), ...
    'mean_new_runtime_s',mean(new.optimization_s), ...
    'mean_runtime_ratio',mean(new.optimization_s)/mean(old.optimization_s), ...
    'old_joint_passes',nnz(old.joint_pass),'new_joint_passes',nnz(new.joint_pass), ...
    'new_own_multiplier_joint_passes',nnz(new.paper_joint_pass), ...
    'maximum_new_own_kkt',max(new.solver_kkt), ...
    'minimum_new_own_tangent_eigenvalue',min(new.paper_lambda_min_tangent), ...
    'maximum_new_cone_tolerance',max(new.solver_critical_cone_tolerance), ...
    'timing_note','Fresh single-thread processes; execute jobs sequentially with other heavy jobs paused.');
assert(report.maximum_relative_crlb_difference<1e-6, ...
    'New solver no longer reproduces old CRLB.');
save(fullfile(root,'audit','alignment_pairs.mat'),'report','T');
fid=fopen(fullfile(root,'audit','alignment_pairs.json'),'w','n','UTF-8');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
disp(report);
end
