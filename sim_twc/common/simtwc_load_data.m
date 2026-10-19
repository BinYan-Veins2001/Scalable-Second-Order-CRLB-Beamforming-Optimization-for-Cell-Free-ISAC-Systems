function data = simtwc_load_data(id)
%SIMTWC_LOAD_DATA Load a figure checkpoint, computing it when absent.

s=simtwc_config(); path=simtwc_archive_path(id);
compatible=false;
if isfile(path)
    data=load(path);
    compatible=isfield(data,'schema_version') && strcmp(data.schema_version,s.schema_version);
    name=sprintf('fig%02d',id);
    if compatible && isfield(s,name) && isfield(data,'methods')
        expected=simtwc_initialize_data(id,s.(name),data.methods);
        compatible=isfield(data,'signature') && strcmp(data.signature,expected.signature);
    end
end
if ~compatible
    if isappdata(groot,'COSTAR_PLOT_ONLY') && getappdata(groot,'COSTAR_PLOT_ONLY')
        error('simtwc:MissingArchive','Plot-only mode requires a compatible archive: %s.',path);
    end
    simtwc_compute_figure(id);
end
data=load(path);
end
