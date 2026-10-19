function archive = simtwc_compute_fig08()
%SIMTWC_COMPUTE_FIG08 Nested compact and large-area AP scaling.

s=simtwc_setup(); spec=s.fig08;
validate_sdr_limit(spec);
compact_methods=s.methods;
large_methods=setdiff(s.methods,{'sdr_sca'},'stable');
method_meta=struct('compact',{compact_methods},'large',{large_methods});
data=simtwc_initialize_data(8,spec,method_meta); path=simtwc_archive_path(8);
data.execution_policy=struct( ...
    'sdr_sca_max_M',spec.sdr_sca_max_M, ...
    'compact_sdr_candidate_M',spec.compact_sdr_M, ...
    'rule','run SDR-SCA only at candidate M with M <= sdr_sca_max_M');
fprintf(['[Fig08] SDR-SCA compact-grid limit: M <= %g; ' ...
    'candidate M = %s.\n'],spec.sdr_sca_max_M,mat2str(spec.compact_sdr_M));

nA=numel(compact_methods); nM=numel(spec.compact.M_values);
nT=numel(spec.seeds_compact);
if ~isfield(data,'compact_runs') || ...
        ~same_size3(data.compact_runs,nA,nM,nT)
    data.compact_runs=cell(nA,nM,nT);
end
for t=1:nT
    mother_spec=spec.compact; mother_spec.M=max(spec.compact.M_values);
    mother=simtwc_make_problem(mother_spec,spec.seeds_compact(t));
    for m=1:nM
        M=spec.compact.M_values(m); p=simtwc_nested_problem(mother,M);
        for a=1:nA
            if strcmp(compact_methods{a},'sdr_sca') && ...
                    (~ismember(M,spec.compact_sdr_M) || ...
                    M>spec.sdr_sca_max_M)
                % Preserve an already completed SDR result, but never launch
                % a new solve outside the current cap.  Failed/empty entries
                % become explicit skips and therefore count as completed.
                if ~successful(data.compact_runs{a,m,t})
                    if M>spec.sdr_sca_max_M
                        reason=sprintf( ...
                            'SDR-SCA skipped: M=%d exceeds configured maximum M=%g', ...
                            M,spec.sdr_sca_max_M);
                    else
                        reason='SDR-SCA skipped: outside sparse evaluation grid';
                    end
                    data.compact_runs{a,m,t}=skipped(reason);
                end
                continue;
            end
            if done(data.compact_runs{a,m,t}); continue; end
            fprintf('[Fig08 compact] M=%d, seed=%d, %s\n', ...
                M,spec.seeds_compact(t),compact_methods{a});
            data.compact_runs{a,m,t}=simtwc_try_run( ...
                compact_methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end

nA=numel(large_methods); nM=numel(spec.large.M_values);
nT=numel(spec.seeds_large);
if ~isfield(data,'large_runs') || ~same_size3(data.large_runs,nA,nM,nT)
    data.large_runs=cell(nA,nM,nT);
end
for t=1:nT
    mother_spec=spec.large; mother_spec.M=max(spec.large.M_values);
    mother=simtwc_make_problem(mother_spec,spec.seeds_large(t));
    for m=1:nM
        M=spec.large.M_values(m); p=simtwc_nested_problem(mother,M);
        for a=1:nA
            if done(data.large_runs{a,m,t}); continue; end
            fprintf('[Fig08 large] M=%d, seed=%d, %s\n', ...
                M,spec.seeds_large(t),large_methods{a});
            data.large_runs{a,m,t}=simtwc_try_run( ...
                large_methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end
data.compact_methods=compact_methods; data.large_methods=large_methods;
[data.compact_final_crlb,data.compact_runtime_s,data.compact_feasible]= ...
    extract_results(data.compact_runs);
[data.large_final_crlb,data.large_runtime_s,data.large_feasible]= ...
    extract_results(data.large_runs);
data.complete=all(cellfun(@done,data.compact_runs(:))) && ...
    all(cellfun(@done,data.large_runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function x=skipped(reason)
x=struct('status','skipped','message',reason,'feasible',false, ...
    'runtime_median_s',NaN,'metrics',struct('crlb',NaN));
end
function yes=done(x)
yes=isstruct(x)&&isfield(x,'status')&&any(strcmp(x.status,{'ok','skipped'}));
end
function yes=successful(x)
yes=isstruct(x)&&isfield(x,'status')&&strcmp(x.status,'ok');
end

function validate_sdr_limit(spec)
if ~isfield(spec,'sdr_sca_max_M') || ...
        ~isscalar(spec.sdr_sca_max_M) || ...
        ~isreal(spec.sdr_sca_max_M) || ...
        isnan(spec.sdr_sca_max_M) || spec.sdr_sca_max_M<0
    error(['s.fig08.sdr_sca_max_M must be a nonnegative real scalar ' ...
        '(or Inf).']);
end
end

function yes=same_size3(x,n1,n2,n3)
yes=size(x,1)==n1 && size(x,2)==n2 && size(x,3)==n3;
end

function [C,T,F]=extract_results(runs)
sz=size(runs); C=nan(sz); T=nan(sz); F=false(sz);
for k=1:numel(runs)
    q=runs{k};
    if ~isstruct(q) || ~isfield(q,'status') || ~strcmp(q.status,'ok'); continue; end
    C(k)=q.metrics.crlb; T(k)=q.runtime_median_s; F(k)=q.feasible;
end
end
