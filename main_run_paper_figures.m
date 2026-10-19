function report=main_run_paper_figures(options)
%MAIN_RUN_PAPER_FIGURES Compute manuscript figures 2--6; figure 1 is opt-in.
% Example: main_run_paper_figures(struct('figures',[2 4 5]));
% Per-figure options: options.fig04=struct('seeds',2026:2035).
if nargin<1; options=struct(); end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
if ~isfield(options,'figures'); options.figures=2:6; end
numbers=options.figures(:).';
assert(all(ismember(numbers,1:6)),'Figure numbers must be between 1 and 6.');
if any(numbers==1)
    assert(isfield(options,'run_dimension_scaling') && options.run_dimension_scaling, ...
        'Figure 1 requires explicit run_dimension_scaling=true.');
end
if ~isfield(options,'continue_on_error'); options.continue_on_error=true; end
entries=paperfig_manifest();
report=repmat(struct('number',0,'entry','','status','','archive','','message',''),0,1);
for number=numbers
    item=entries(number); key=sprintf('fig%02d',number); current=struct();
    if isfield(options,key); current=options.(key); end
    row=struct('number',number,'entry',item.entry,'status','running', ...
        'archive','','message','');
    fprintf('[Paper figure %d] %s\n',number,item.name);
    try
        value=feval(item.entry,current);
        if ischar(value) || isstring(value)
            row.archive=char(value);
        else
            row.archive=fullfile(root,'results',sprintf('paper_fig%02d',number),'data.mat');
        end
        row.status='completed';
    catch exception
        row.status='error'; row.message=exception.message;
        fprintf(2,'%s\n',getReport(exception,'extended','hyperlinks','off'));
        if ~options.continue_on_error; rethrow(exception); end
    end
    report(end+1)=row; %#ok<AGROW>
    folder=fullfile(root,'results'); if ~exist(folder,'dir'); mkdir(folder); end
    save(fullfile(folder,'paper_run_report.mat'),'report','options');
    writetable(struct2table(report,'AsArray',true),fullfile(folder,'paper_run_report.csv'));
end
end
