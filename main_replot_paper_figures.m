function report=main_replot_paper_figures(numbers)
%MAIN_REPLOT_PAPER_FIGURES Replot stored new-version data without simulation.
if nargin<1; numbers=2:6; end
root=fileparts(mfilename('fullpath')); addpath(genpath(root),'-begin');
entries=paperfig_manifest();
report=repmat(struct('number',0,'archive','','status','','message',''),0,1);
for number=numbers(:).'
    assert(ismember(number,1:6),'Unknown manuscript figure.');
    archive=fullfile(root,'results',sprintf('paper_fig%02d',number),'data.mat');
    row=struct('number',number,'archive',archive,'status','missing','message','');
    if exist(archive,'file')
        try
            feval(entries(number).plotter,archive); row.status='completed';
        catch exception
            row.status='error'; row.message=exception.message;
            fprintf(2,'Figure %d: %s\n',number,exception.message);
        end
    else
        row.message='No data archive; simulation was not launched.';
        fprintf('[Paper figure %d] No stored data; skipping.\n',number);
    end
    report(end+1)=row; %#ok<AGROW>
end
end
