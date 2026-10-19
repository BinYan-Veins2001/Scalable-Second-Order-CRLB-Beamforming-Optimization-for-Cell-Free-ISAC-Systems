function archive = simtwc_compute_fig15()
%SIMTWC_COMPUTE_FIG15 AP, transmit-antenna, and user-count scaling.
%
% The experiment is fully controlled by s.fig15 in simtwc_config.m.  The
% three sweeps use paired Monte Carlo seeds.  SDR-SCA is resource-capped only
% in the AP sweep; every configured method runs throughout the Nt and U grids.

s=simtwc_setup(); spec=s.fig15; methods=cellstr(string(spec.methods));
validate_spec(spec);
data=simtwc_initialize_data(15,spec,methods); path=simtwc_archive_path(15);
nA=numel(methods); nT=numel(spec.seeds);

%% Column 1: number of APs, with a nested maximum-M realization per seed.
M_values=spec.ap.M_values;
nM=numel(M_values);
if ~isfield(data,'ap_runs') || ~same_size3(data.ap_runs,nA,nM,nT)
    data.ap_runs=cell(nA,nM,nT);
end
data.execution_policy=struct('sdr_sca_max_M',spec.sdr_sca_max_M, ...
    'rule','AP sweep runs SDR-SCA only for M <= the inclusive cap');
fprintf('[Fig15] AP-sweep SDR-SCA limit: M <= %g.\n',spec.sdr_sca_max_M);
for t=1:nT
    mother_spec=common_fields(spec);
    mother_spec.M=max(M_values); mother_spec.U=spec.ap.U;
    mother_spec.Nt=spec.ap.Nt; mother_spec.Nr=spec.ap.Nr;
    mother=simtwc_make_problem(mother_spec,spec.seeds(t));
    for m=1:nM
        M=M_values(m); p=simtwc_nested_problem(mother,M);
        for a=1:nA
            is_sdr=strcmp(methods{a},'sdr_sca');
            if is_sdr && M>spec.sdr_sca_max_M
                % Preserve a previously completed solve in the MAT archive,
                % but mark the current point as disallowed for plotting.
                if ~ok(data.ap_runs{a,m,t})
                    data.ap_runs{a,m,t}=skipped(sprintf( ...
                        'SDR-SCA skipped: M=%d exceeds configured cap M=%g', ...
                        M,spec.sdr_sca_max_M));
                end
                continue;
            end
            % A formerly skipped point is recomputed automatically if the
            % user later raises the cap; only successful results are reused.
            if ok(data.ap_runs{a,m,t}); continue; end
            fprintf('[Fig15 AP] M=%d, seed=%d, %s\n', ...
                M,spec.seeds(t),methods{a});
            data.ap_runs{a,m,t}=simtwc_try_run(methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end

%% Column 2: number of transmit antennas per AP.
Nt_values=spec.antenna.Nt_values;
nN=numel(Nt_values);
if ~isfield(data,'antenna_runs') || ~same_size3(data.antenna_runs,nA,nN,nT)
    data.antenna_runs=cell(nA,nN,nT);
end
for n=1:nN
    pspec=common_fields(spec);
    pspec.M=spec.antenna.M; pspec.U=spec.antenna.U;
    pspec.Nt=Nt_values(n);
    if spec.antenna.couple_receive_antennas
        pspec.Nr=Nt_values(n);
    else
        pspec.Nr=spec.antenna.Nr;
    end
    for t=1:nT
        p=simtwc_make_problem(pspec,spec.seeds(t));
        for a=1:nA
            if ok(data.antenna_runs{a,n,t}); continue; end
            fprintf('[Fig15 Nt] Nt=%d, Nr=%d, seed=%d, %s\n', ...
                pspec.Nt,pspec.Nr,spec.seeds(t),methods{a});
            data.antenna_runs{a,n,t}=simtwc_try_run( ...
                methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end

%% Column 3: number of users.
U_values=spec.user.U_values;
nU=numel(U_values);
if ~isfield(data,'user_runs') || ~same_size3(data.user_runs,nA,nU,nT)
    data.user_runs=cell(nA,nU,nT);
end
for u=1:nU
    pspec=common_fields(spec);
    pspec.M=spec.user.M; pspec.U=U_values(u);
    pspec.Nt=spec.user.Nt; pspec.Nr=spec.user.Nr;
    for t=1:nT
        p=simtwc_make_problem(pspec,spec.seeds(t));
        for a=1:nA
            if ok(data.user_runs{a,u,t}); continue; end
            fprintf('[Fig15 users] U=%d, seed=%d, %s\n', ...
                pspec.U,spec.seeds(t),methods{a});
            data.user_runs{a,u,t}=simtwc_try_run(methods{a},p,false,struct());
            simtwc_save_atomic(path,data);
        end
    end
end

data.methods=methods; data.seeds=spec.seeds;
data.M_values=M_values; data.Nt_values=Nt_values; data.U_values=U_values;
data.ap_sdr_allowed=M_values<=spec.sdr_sca_max_M;
[data.ap_crlb,data.ap_runtime_s,data.ap_feasible,data.ap_kkt]= ...
    extract(data.ap_runs);
[data.antenna_crlb,data.antenna_runtime_s, ...
    data.antenna_feasible,data.antenna_kkt]=extract(data.antenna_runs);
[data.user_crlb,data.user_runtime_s,data.user_feasible,data.user_kkt]= ...
    extract(data.user_runs);
data.complete=all(cellfun(@terminal,data.ap_runs(:))) && ...
    all(cellfun(@terminal,data.antenna_runs(:))) && ...
    all(cellfun(@terminal,data.user_runs(:)));
data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end

function p=common_fields(spec)
p=struct('gamma_dB',spec.gamma_dB, ...
    'min_ap_distance',spec.min_ap_distance);
end

function [C,T,F,K]=extract(runs)
sz=size(runs); C=nan(sz); T=nan(sz); F=false(sz); K=nan(sz);
for k=1:numel(runs)
    q=runs{k}; if ~ok(q); continue; end
    C(k)=q.metrics.crlb; T(k)=q.runtime_median_s;
    F(k)=q.feasible; K(k)=q.kkt.r_kkt;
end
end

function q=skipped(message)
q=struct('status','skipped','message',message,'feasible',false, ...
    'runtime_median_s',NaN,'metrics',struct('crlb',NaN), ...
    'kkt',struct('r_kkt',NaN));
end

function validate_spec(spec)
if ~isscalar(spec.sdr_sca_max_M) || ~isreal(spec.sdr_sca_max_M) || ...
        isnan(spec.sdr_sca_max_M) || spec.sdr_sca_max_M<0
    error('s.fig15.sdr_sca_max_M must be a nonnegative scalar or Inf.');
end
assert(isvector(spec.ap.M_values) && all(diff(spec.ap.M_values)>0), ...
    's.fig15.ap.M_values must be a strictly increasing vector.');
assert(isvector(spec.antenna.Nt_values) && ...
    all(diff(spec.antenna.Nt_values)>0), ...
    's.fig15.antenna.Nt_values must be a strictly increasing vector.');
assert(isvector(spec.user.U_values) && all(diff(spec.user.U_values)>0), ...
    's.fig15.user.U_values must be a strictly increasing vector.');
end

function yes=same_size3(x,n1,n2,n3)
yes=size(x,1)==n1 && size(x,2)==n2 && size(x,3)==n3;
end

function yes=terminal(q)
yes=isstruct(q)&&isfield(q,'status')&&any(strcmp(q.status,{'ok','skipped'}));
end

function yes=ok(q)
yes=isstruct(q)&&isfield(q,'status')&&strcmp(q.status,'ok');
end
