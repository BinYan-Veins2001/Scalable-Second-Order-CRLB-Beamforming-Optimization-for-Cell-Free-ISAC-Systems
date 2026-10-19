function path = simtwc_archive_path(id)
%SIMTWC_ARCHIVE_PATH Canonical MAT checkpoint path for one figure.

s=simtwc_config();
suffix=''; if s.quick; suffix='_quick'; end
path=fullfile(s.data_dir,sprintf('fig%02d_data%s.mat',id,suffix));
end
