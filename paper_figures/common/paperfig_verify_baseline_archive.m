function report=paperfig_verify_baseline_archive()
%PAPERFIG_VERIFY_BASELINE_ARCHIVE Read-only initialization and metric checks.
s=paperfig_setup(); source=fullfile(fileparts(s.root), ...
    'COSTAR_TWC_Simulations_ExactHVP_PTR_20260907','data','fig03_data.mat');
old=load(source); rows=cell(0,8);
for t=1:numel(s.seeds)
    seed=s.seeds(t); j=find(old.seeds==seed,1); assert(~isempty(j));
    p=simtwc_make_problem(s.default,seed); z0=complex_matrix_to_real_vector(p.initialization.Wmat);
    for method=s.methods
        key=method{1}; if strcmp(key,'costar'); continue; end
        a=find(strcmp(old.methods,key),1); run=old.runs{a,j};
        initial_error=NaN; initial_found=false;
        if any(strcmp(key,s.iterative_methods))
            if isfield(run.solution,'initial_Wmat')
                zi=complex_matrix_to_real_vector(run.solution.initial_Wmat); initial_found=true;
            elseif isfield(run.history,'z_snapshots') && ~isempty(run.history.z_snapshots)
                zi=run.history.z_snapshots(:,1); initial_found=true;
            end
            if initial_found; initial_error=norm(zi-z0)/max(1,norm(z0)); end
        else
            initial_found=true; initial_error=0;
        end
        metrics=evaluate_beamforming_solution(p.cfg,p.channels,p.smat,matrix_to_wcell(p.cfg,run.Wmat));
        crlb_error=abs(metrics.crlb-run.metrics.crlb)/max(abs(run.metrics.crlb),realmin);
        sinr_error=norm(metrics.sinr-run.metrics.sinr)/max(1,norm(run.metrics.sinr));
        pass=initial_found && initial_error<=1e-12 && crlb_error<=1e-12 && sinr_error<=1e-12;
        rows(end+1,:)={key,seed,initial_found,initial_error,crlb_error,sinr_error,pass,source}; %#ok<AGROW>
    end
end
report=cell2table(rows,'VariableNames',{'method','seed','initial_found','initial_error', ...
    'CRLB_error','SINR_error','eligible_initial_model','source'});
folder=fullfile(s.results_dir,'paper_fig03'); if ~isfolder(folder); mkdir(folder); end
writetable(report,fullfile(folder,'archive_reuse_verification.csv')); disp(report);
end
