% SETUP_PROJECT Add the replication code to the MATLAB path.

project_root = fileparts(mfilename('fullpath'));
code_dir = fullfile(project_root, 'code');
helpers_dir = fullfile(code_dir, 'helpers');
generated_dir = fullfile(project_root, 'results', 'generated');

addpath(code_dir, helpers_dir);

if exist(generated_dir, 'dir') ~= 7
    mkdir(generated_dir);
end

cd(project_root);
fprintf('Replication project ready: %s\n', project_root);
