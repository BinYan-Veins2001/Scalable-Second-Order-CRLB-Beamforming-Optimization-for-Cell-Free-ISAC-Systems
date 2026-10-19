function data = simtwc_initialize_data(id,spec,methods)
%SIMTWC_INITIALIZE_DATA Resume only schema-compatible checkpoints.
%
% Figures 8 and 15 use an SDR-SCA AP cap as an execution-resource safeguard
% rather than a physical experiment parameter. Excluding it from the
% compatibility signature lets the user lower or raise the cap without
% discarding an already large AP-scaling checkpoint.

s=simtwc_config();
path=simtwc_archive_path(id);
signature_spec=spec;
if any(id==[8 15]) && isfield(signature_spec,'sdr_sca_max_M')
    signature_spec=rmfield(signature_spec,'sdr_sca_max_M');
end
signature=jsonencode(struct('spec',signature_spec,'methods',{methods}, ...
    'costar_options',s.costar_options,'initialization',s.initialization_policy, ...
    'runtime',s.runtime_definition,'defaults',default_config(), ...
    'kkt_rule',[s.kkt_sinr_active_tol,s.kkt_power_active_tol]));
if s.resume && isfile(path)
    old=load(path);
    if isfield(old,'schema_version') && strcmp(old.schema_version,s.schema_version) && ...
            isfield(old,'signature') && strcmp(old.signature,signature)
        data=old;
        % Keep the archived numerical results while recording the currently
        % requested resource policy in the next checkpoint save.
        data.spec=spec;
        return;
    end
end
data=struct('schema_version',s.schema_version,'code_version',s.code_version, ...
    'figure_id',id,'signature',signature,'spec',spec,'methods',{methods}, ...
    'quick',s.quick,'complete',false, ...
    'created_at',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
end
