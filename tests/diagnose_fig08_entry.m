function report=diagnose_fig08_entry()
%DIAGNOSE_FIG08_ENTRY Distinguish ZF construction failure from SINR infeasibility.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
previous=getenv('COSTAR_SMOKE'); quick=getenv('COSTAR_QUICK');
cleanup=onCleanup(@() restore(previous,quick)); %#ok<NASGU>
setenv('COSTAR_SMOKE','1'); setenv('COSTAR_QUICK','');
s=simtwc_config(); spec=s.fig08.large; spec.M=max(spec.M_values);
mother=simtwc_make_problem(spec,2001); report=struct();
for M=[4 6]
    p=simtwc_nested_problem(mother,M); entry=p.initialization;
    cfg=p.cfg; U=cfg.sys.U; N=M*cfg.sys.Nt;
    H=p.channels.h_stack';
    A=H*sqrt(cfg.power.Pmax/cfg.comm.noise_power);
    upper=zeros(U,1);
    for u=1:U
        gain=0;
        for m=1:M; gain=gain+norm(p.channels.h{m,u}); end
        upper(u)=cfg.power.Pmax*gain^2/cfg.comm.noise_power;
    end
    cap=1-cfg.costar.global_zf_init.power_margin_fraction;
    required=sqrt(cfg.comm.Gamma*(1+cfg.costar.global_zf_init.sinr_margin_fraction));
    cvx_begin quiet
        cvx_solver sdpt3
        variable Xzf(N,U) complex
        variable tzf nonnegative
        minimize(tzf)
        subject to
            A*Xzf == required*eye(U);
            for m=1:M
                rows=(m-1)*cfg.sys.Nt+(1:cfg.sys.Nt);
                norm(reshape(Xzf(rows,:),[],1))<=tzf;
            end
    cvx_end
    zf_status=cvx_status; zf_peak=tzf^2;
    cvx_begin quiet
        cvx_solver sdpt3
        variable Xsinr(N,U) complex
        variable tsinr nonnegative
        minimize(tsinr)
        subject to
            for u=1:U
                imag(A(u,:)*Xsinr(:,u))==0;
                real(A(u,:)*Xsinr(:,u))>=sqrt(cfg.comm.Gamma)* ...
                    norm([A(u,:)*Xsinr(:,setdiff(1:U,u)),1]);
            end
            for m=1:M
                rows=(m-1)*cfg.sys.Nt+(1:cfg.sys.Nt);
                norm(reshape(Xsinr(rows,:),[],1))<=tsinr;
            end
    cvx_end
    sinr_status=cvx_status; sinr_peak=tsinr^2;
    delivered=A*Xsinr;
    achieved=zeros(U,1);
    for u=1:U
        achieved(u)=abs(delivered(u,u))^2/ ...
            (1+sum(abs(delivered(u,setdiff(1:U,u))).^2));
    end
    q=struct('M',M,'seed',2001,'initialization',entry, ...
        'single_user_snr_upper_dB',10*log10(upper), ...
        'zf_min_peak_power_ratio',zf_peak,'zf_cvx_status',zf_status, ...
        'zf_power_cap_ratio',cap,'sinr_min_peak_power_ratio',sinr_peak, ...
        'sinr_cvx_status',sinr_status,'sinr_candidate_dB',10*log10(achieved), ...
        'positions',p.scenario);
    report.(sprintf('M%d',M))=q;
    info=entry.info;
    fprintf('\nM=%d U=%d Nt=%d seed=2001, entry=%d, rank=%d, cond=%.5g\n', ...
        M,U,cfg.sys.Nt,entry.success,info.aggregate_rank,info.condition_number);
    fprintf('floor AP powers / Pmax: %s\n',mat2str(info.floor_ap_power'/cfg.power.Pmax,7));
    fprintf('projection cycles=%d converged=%d normalized response residual=%.6g\n', ...
        info.projection_info.cycles,info.projection_info.converged,info.response_residual);
    fprintf('actual SINR(dB): %s; AP power/Pmax: %s; CRLB %.6g\n', ...
        mat2str(10*log10(info.metrics.sinr(:))',7), ...
        mat2str(info.metrics.ap_power(:)'/cfg.power.Pmax,7),info.metrics.crlb);
    fprintf('single-user SNR upper(dB): %s\n',mat2str(q.single_user_snr_upper_dB',7));
    fprintf('ZF minimum peak/Pmax=%.9g [cap %.9g], CVX %s\n',zf_peak,cap,zf_status);
    fprintf('general SINR minimum peak/Pmax=%.9g, CVX %s, achieved SINR %s\n', ...
        sinr_peak,sinr_status,mat2str(q.sinr_candidate_dB',7));
end
% Initialization-only audit of the actual publication large-area settings.
setenv('COSTAR_SMOKE',''); s=simtwc_config(); spec=s.fig08.large;
spec.M=max(spec.M_values); mother=simtwc_make_problem(spec,2001);
p=simtwc_nested_problem(mother,min(spec.M_values));
report.publication_entry=p.initialization;
fprintf('\nPUBLICATION M=%d U=%d Nt=%d: shared entry success=%d, SINRmin=%.7g dB\n', ...
    p.cfg.sys.M,p.cfg.sys.U,p.cfg.sys.Nt,p.initialization.success, ...
    10*log10(min(p.initialization.metrics.sinr)));
save(fullfile(root,'tests','fig08_entry_diagnosis.mat'),'report','-v7.3');
end

function restore(previous,quick)
setenv('COSTAR_SMOKE',previous); setenv('COSTAR_QUICK',quick);
end
