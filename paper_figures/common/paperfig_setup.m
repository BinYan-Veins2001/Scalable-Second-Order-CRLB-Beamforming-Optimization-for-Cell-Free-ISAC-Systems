function s=paperfig_setup()
%PAPERFIG_SETUP Current-paper output paths, with all solver calls local.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(root),'-begin');
assert(~strcmpi(strtrim(getenv('COSTAR_SMOKE')),'1') && ...
    ~strcmpi(strtrim(getenv('COSTAR_QUICK')),'1'), ...
    'PaperFigures:LegacyProfile', ...
    'Clear COSTAR_SMOKE/COSTAR_QUICK for paper figures; use explicit options.seeds for a shorter run.');
s=simtwc_config(); s.root=root;
s.figure_dir=fullfile(root,'figures');
s.results_dir=fullfile(root,'results');
s.result_dir=s.results_dir;
s.case_dir=fullfile(s.results_dir,'shared_cases');
s.schema_version='paper_figures_safeguarded_spsc_20261019_v1';
s.seeds=2026:2035;
s.epsilon_KKT=1e-5; s.epsilon_H=1e-5;
for p={s.figure_dir,s.results_dir,s.case_dir}
    if ~isfolder(p{1}); mkdir(p{1}); end
end
resolved=which('solve_costar');
assert(startsWith(strrep(resolved,'\','/'),strrep(root,'\','/')), ...
    'PaperFigures:WrongSolver','solve_costar must resolve inside this package.');
end
