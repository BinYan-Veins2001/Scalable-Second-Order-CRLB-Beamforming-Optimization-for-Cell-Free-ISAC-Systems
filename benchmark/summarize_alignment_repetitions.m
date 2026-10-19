function report=summarize_alignment_repetitions()
%SUMMARIZE_ALIGNMENT_REPETITIONS Use all runs, with median per paired input.
root=fileparts(fileparts(mfilename('fullpath')));
vroot=fullfile(root,'benchmark','validation');
files=strip(splitlines(fileread(fullfile(vroot,'jobs_all_repetitions.txt'))));
files=files(strlength(files)>0); rows=cell(numel(files),1);
references=struct(); consistency=struct();
for j=1:numel(files)
    f=char(files(j));
    if endsWith(f,'.json'); job=jsondecode(fileread(f));
    else; input=load(f,'job'); job=input.job; end
    data=load(job.output_file,'result','sol','hist'); row=data.result.row;
    if isfield(job,'repetition'); repetition=job.repetition; else; repetition=1; end
    row.repetition=repetition; rows{j}=row;
    assert(row.executed && row.joint_pass,'Endpoint audit failed.');
    if strcmp(row.label,'new')
        assert(data.sol.success && data.sol.certificate.oracle_is_at_final_point);
        assert(row.paper_joint_pass && row.paper_primal<=1e-10);
        assert(data.sol.critical_cone_tolerance<=data.sol.options.epsilon_H);
        assert(data.sol.kkt_residual<=data.sol.options.epsilon_KKT);
        assert(data.sol.work.fd_hvp_calls==0 && data.sol.work.meo_calls==0);
    end
    key=sprintf('%s_%d',row.label,row.seed);
    if repetition==1
        references.(key)=data;
        consistency.(key)=struct('same_work_counts',true, ...
            'same_endpoint',true,'same_numerical_trace',true, ...
            'maximum_endpoint_difference',0,'compared_repeats',0);
    else
        base=references.(key); c=consistency.(key);
        c.same_work_counts=c.same_work_counts && isequaln(base.sol.work,data.sol.work);
        c.same_endpoint=c.same_endpoint && isequaln(base.sol.Wmat,data.sol.Wmat);
        c.maximum_endpoint_difference=max(c.maximum_endpoint_difference, ...
            norm(base.sol.Wmat-data.sol.Wmat,'fro'));
        fields={'objective','crlb','kkt','stationarity','delta','accepted', ...
            'step_type','cg_iter','z_snapshots','lambda_min'};
        for k=1:numel(fields)
            name=fields{k};
            if isfield(base.hist.rtr,name) && isfield(data.hist.rtr,name)
                c.same_numerical_trace=c.same_numerical_trace && ...
                    isequaln(base.hist.rtr.(name),data.hist.rtr.(name));
            end
        end
        c.compared_repeats=c.compared_repeats+1; consistency.(key)=c;
    end
end
T=struct2table(vertcat(rows{:}));
writetable(T,fullfile(vroot,'trial_results_all_repetitions.csv'));
seeds=unique(T.seed); medians=zeros(numel(seeds),5);
for j=1:numel(seeds)
    old=T(T.seed==seeds(j) & strcmp(T.label,'old'),:);
    new=T(T.seed==seeds(j) & strcmp(T.label,'new'),:);
    assert(height(old)==3 && height(new)==3,'Exactly three repetitions required.');
    mo=median(old.optimization_s); mn=median(new.optimization_s);
    medians(j,:)=[seeds(j),mo,mn,mn/mo, ...
        abs(new.crlb(1)-old.crlb(1))/old.crlb(1)];
end
P=array2table(medians,'VariableNames',{'seed','old_median_s', ...
    'new_median_s','paired_runtime_ratio','relative_crlb_difference'});
writetable(P,fullfile(vroot,'paired_repetition_medians.csv'));
report=struct('passed',true,'paired_inputs',height(P),'repetitions_per_solver',3, ...
    'fresh_processes',height(T),'threads',1, ...
    'mean_old_median_s',mean(P.old_median_s), ...
    'mean_new_median_s',mean(P.new_median_s), ...
    'ratio_of_mean_medians',mean(P.new_median_s)/mean(P.old_median_s), ...
    'mean_paired_runtime_ratio',mean(P.paired_runtime_ratio), ...
    'maximum_relative_crlb_difference',max(P.relative_crlb_difference), ...
    'old_joint_passes',nnz(T.joint_pass & strcmp(T.label,'old')), ...
    'new_joint_passes',nnz(T.joint_pass & strcmp(T.label,'new')), ...
    'new_own_multiplier_joint_passes',nnz(T.paper_joint_pass & strcmp(T.label,'new')), ...
    'consistency',consistency, ...
    'aggregation','Median of all three fresh-process repeats per solver/seed; mean across three paired seeds. Initial runs retained.');
assert(report.maximum_relative_crlb_difference<1e-6);
names=fieldnames(consistency);
for j=1:numel(names)
    c=consistency.(names{j});
    assert(c.compared_repeats==2 && c.same_work_counts && c.same_endpoint && c.same_numerical_trace, ...
        'A repeated run changed its numerical work or trajectory.');
end
save(fullfile(root,'audit','alignment_repetitions.mat'),'report','T','P');
fid=fopen(fullfile(root,'audit','alignment_repetitions.json'),'w','n','UTF-8');
assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true)); disp(P); disp(report);
end
