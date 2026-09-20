function Run_All_Figures
%RUN_ALL_FIGURES Reproduce numerical Figures 2--9 in manuscript order.
%
% Figure 1 is the system-model schematic and is not MATLAB-generated.
% Figure 2 generates/overwrites Master_Numerical_Data.mat. Later scripts
% load and append to that shared file. Each script exports Figure_0N.pdf
% and Figure_0N.png into this repository folder.

rootDir = fileparts(mfilename('fullpath'));
scripts = { ...
    'Code_Figure_02.m', ...
    'Code_Figure_03.m', ...
    'Code_Figure_04.m', ...
    'Code_Figure_05.m', ...
    'Code_Figure_06.m', ...
    'Code_Figure_07.m', ...
    'Code_Figure_08.m', ...
    'Code_Figure_09.m'};

fprintf('\nReproducing Figures 2--9 from: %s\n', rootDir);
for k = 1:numel(scripts)
    scriptPath = fullfile(rootDir, scripts{k});
    if ~isfile(scriptPath)
        error('Missing required script: %s', scriptPath);
    end
    fprintf('\n[%d/%d] Running %s ...\n', k, numel(scripts), scripts{k});
    escapedPath = strrep(scriptPath, , '');
    evalin('base', sprintf('run(''%s'');', escapedPath));
    close all force;
end
fprintf('
Completed. Figures 2--9 were regenerated.
');
end
