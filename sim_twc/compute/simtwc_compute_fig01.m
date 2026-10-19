function archive = simtwc_compute_fig01()
%SIMTWC_COMPUTE_FIG01 Topology and SDR-SCA/COSTAR position-CRLB maps.

s=simtwc_setup(); spec=s.fig01; keys=spec.map_methods;
data=simtwc_initialize_data(1,spec,keys);
path=simtwc_archive_path(1);
if data.complete; archive=path; return; end
p=simtwc_make_problem(spec,spec.seed);
data.scenario=p.scenario;
data.cfg=p.cfg;
data.runs=cell(1,numel(keys));
for k=1:numel(keys)
    fprintf('[Fig01] %s\n',keys{k});
    data.runs{k}=simtwc_try_run(keys{k},p,false,struct());
    simtwc_save_atomic(path,data);
end
x=spec.grid_axis; y=spec.grid_axis;
maps=nan(numel(y),numel(x),numel(keys));
for iy=1:numel(y)
    for ix=1:numel(x)
        scenario_q=p.scenario;
        scenario_q.target_pos=[x(ix),y(iy)];
        smat_q=build_sensing_matrices(p.cfg,scenario_q);
        for k=1:numel(keys)
            if strcmp(data.runs{k}.status,'ok')
                maps(iy,ix,k)=compute_crlb_from_w( ...
                    p.cfg,smat_q,data.runs{k}.wcell);
            end
        end
    end
end
data.grid_x=x; data.grid_y=y; data.crlb_maps=maps;
data.map_methods=keys;
data.complete=true; data.completed_at=char(datetime('now'));
simtwc_save_atomic(path,data); archive=path;
end
