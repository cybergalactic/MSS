function mssVessels()
% mssVessels lists marine craft in MSS with a MAT-file containing the vessel structure.
%
%   mssVessels      Print source, filename, and location in the Command Window.
%
% These MAT-files contain the MSS vessel structure produced by hydrodynamic
% solvers and their import/export tools, such as ShipX/VERES, WAMIT, and
% MSS-Capytaine.
%
% Searches the HYDRO vessel catalogues in this MSS installation, including
% subfolders. Companion *ABC.mat files are omitted. Other MAT-files outside
% these catalogues, such as model-specific RAOs and FDI data, are not listed.
% The files are listed without loading their contents. This function has no
% output and does not create or change variables in the caller's workspace.

mssRoot = fileparts(mfilename('fullpath'));
catalogues = {'vessels_shipx', 'vessels_wamit', 'vessels_capytaine'};
sources = {'ShipX / VERES', 'WAMIT', 'Capytaine'};
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
