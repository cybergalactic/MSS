function mssCraft()
% mssCraft lists vessel-structure MAT-files and CRAFT models in MSS.
%
%   mssCraft     Print the MSS craft in the Command Window.
%
% These MAT-files contain the MSS vessel structure produced by hydrodynamic
% solvers and their import/export tools.
%
% Searches the HYDRO vessel catalogues in this MSS installation, including
% subfolders. Other MAT-files outside these catalogues, such as model-specific
% RAOs are not listed. The files are listed without loading their contents.
%
% Author:    Thor I. Fossen
% Date:      2026-10-28
% Revisions:

mssRoot = fileparts(mfilename('fullpath'));
catalogues = {'vessels_wamit', 'vessels_capytaine'};
sources = {'WAMIT', 'Capytaine'};
vessels = struct('source', {}, 'name', {}, 'file', {});

for i = 1:numel(catalogues)
    folder = fullfile(mssRoot, 'HYDRO', catalogues{i});
    vessels = [vessels scanCatalogue(folder, sources{i})]; %#ok<AGROW>
end

fprintf('\nMarine craft in MSS with vessel-structure MAT-files (%d files)\n', numel(vessels));
fprintf('Locations relative to: %s\n\n', mssRoot);
fprintf('%-15s %-22s %s\n', 'Source', 'MAT-file', 'Location');
for i = 1:numel(vessels)
    relativeFile = vessels(i).file(length(mssRoot) + 2:end);
    fprintf('%-15s %-22s %s\n', vessels(i).source, ...
        [vessels(i).name '.mat'], relativeFile);
end

craftRoot = fullfile(mssRoot, 'CRAFT');
models = struct('type', {}, 'name', {}, 'file', {}, 'simulation', {});
craftTypes = {'AUV', 'USV', 'SHIP'};
for i = 1:numel(craftTypes)
    folder = fullfile(craftRoot, craftTypes{i});
    files = dir(fullfile(folder, 'models', '*.m'));
    [~, order] = sort({files.name});
    for j = order
        models(end + 1) = craftModel(fullfile(folder, 'models', ...
            files(j).name), craftTypes{i}, folder); %#ok<AGROW>
    end
end

% These dynamics models sit outside the main models directories.
models(end + 1) = craftModel(fullfile(craftRoot, 'SHIP', 'models', ...
    'semisubModels', 'rig.m'), 'Floating', fullfile(craftRoot, 'SHIP'));
models(end + 1) = craftModel(fullfile(craftRoot, 'hydroVessel.m'), ...
    'General', craftRoot);

fprintf('\nCRAFT dynamics models (%d models)\n\n', numel(models));
fprintf('%-10s %-18s %-20s %s\n', 'Type', 'Model', 'Simulation', 'Location');
for i = 1:numel(models)
    relativeFile = models(i).file(length(mssRoot) + 2:end);
    fprintf('%-10s %-18s %-20s %s\n', models(i).type, ...
        [models(i).name '.m'], models(i).simulation, relativeFile);
end

end

function model = craftModel(file, craftType, simulationFolder)
[~, name] = fileparts(file);
simulationName = name;
switch lower(name)
    case 'lcontainer'
        simulationName = 'container';
    case 'otter_legacy'
        simulationName = 'otter';
    case 'rig'
        simulationName = 'semisub';
end
scripts = dir(fullfile(simulationFolder, 'SIM*.m'));
match = find(strcmpi({scripts.name}, ['SIM' simulationName '.m']), 1);
simulation = '-';
if ~isempty(match)
    simulation = scripts(match).name;
end
model = struct('type', craftType, 'name', name, 'file', file, ...
    'simulation', simulation);
end

function vessels = scanCatalogue(folder, source)
% Use ordinary dir calls for MATLAB and GNU Octave compatibility.
vessels = struct('source', {}, 'name', {}, 'file', {});
if exist(folder, 'dir') ~= 7
    return
end
entries = dir(folder);
[~, order] = sort({entries.name});
entries = entries(order);
for i = 1:numel(entries)
    name = entries(i).name;
    if isempty(name) || name(1) == '.'
        continue
    end
    file = fullfile(folder, name);
    if entries(i).isdir
        vessels = [vessels scanCatalogue(file, source)]; %#ok<AGROW>
    else
        [~, stem, extension] = fileparts(name);
        if strcmpi(extension, '.mat') && isempty(regexpi(stem, 'ABC$', 'once'))
            vessels(end + 1) = struct('source', source, ...
                'name', stem, 'file', file); %#ok<AGROW>
        end
    end
end
end
