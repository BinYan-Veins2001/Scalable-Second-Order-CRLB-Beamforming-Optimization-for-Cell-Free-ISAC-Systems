function p = simtwc_nested_problem(mother,M)
%SIMTWC_NESTED_PROBLEM Restrict a maximum-M realization to its first M APs.

cfg=mother.cfg; cfg.sys.M=M;
scenario=mother.scenario;
scenario.ap_pos=scenario.ap_pos(1:M,:); scenario.M=M;
channels=mother.channels;
channels.h=channels.h(1:M,:);
channels.h_stack=channels.h_stack(1:M*cfg.sys.Nt,:);
for name={'large_scale','dist','theta'}
    channels.(name{1})=channels.(name{1})(1:M,:);
end
smat=build_sensing_matrices(cfg,scenario);
p=struct('cfg',cfg,'scenario',scenario,'channels',channels,'smat',smat, ...
    'seed',mother.seed,'spec',mother.spec);
p.spec.M=M;
p.initialization=initialize_shared_stacked_zf(cfg,channels,smat);
end
