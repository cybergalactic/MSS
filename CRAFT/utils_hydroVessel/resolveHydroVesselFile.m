function filePath = resolveHydroVesselFile(matFile)
% resolveHydroVesselFile locates an MSS hydrodynamic vessel MAT-file.
% filePath = resolveHydroVesselFile(matFile) first searches the HYDRO vessel
% catalogues belonging to the same MSS installation as this function. It then
% falls back to the MATLAB or GNU Octave search path. This makes the vessel GUI
% work immediately after new catalogue directories are added, without requiring
% the user to rebuild a previously saved MSS path.
%
% Input:
%   matFile: Vessel data filename, for example 'testShip.mat'
%
% Output:
%   filePath: Absolute path to the vessel data file
%
% Author:    Thor I. Fossen
% Date:      2026-10-04

if ~ischar(matFile) || isempty(matFile)
    error('The vessel data filename must be a nonempty character vector.');
end

utilsPath = fileparts(mfilename('fullpath'));
mssRoot = fileparts(fileparts(utilsPath));
hydroPath = fullfile(mssRoot, 'HYDRO');
catalogues = {'vessels_shipx', 'vessels_wamit', 'vessels_capytaine'};

for i = 1:numel(catalogues)
    cataloguePath = fullfile(hydroPath, catalogues{i});

    % Also support a MAT-file placed directly in a catalogue directory.
    candidate = fullfile(cataloguePath, matFile);
    if exist(candidate, 'file') == 2
        filePath = candidate;
        return
    end

    entries = dir(cataloguePath);
    for j = 1:numel(entries)
        if entries(j).isdir && ...
                ~strcmp(entries(j).name, '.') && ...
                ~strcmp(entries(j).name, '..')
            candidate = fullfile(cataloguePath, entries(j).name, matFile);
            if exist(candidate, 'file') == 2
                filePath = candidate;
                return
            end
        end
    end
end

filePath = which(matFile);
if isempty(filePath)
    error('Unable to find vessel data file "%s".', matFile);
end

end
