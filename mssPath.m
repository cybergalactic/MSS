% MSS Path Update Script
% This script checks for the presence of the 'MSS' directory on MATLAB's 
% path. If found, it removes all existing instances of this path and its 
% sub-directories and then reinstalls them to ensure the path includes all 
% current sub-directories. This is useful for situations where the directory 
% structure of 'MSS' may change due to updates or modifications in the 
% file system.
%
% Usage:
% - Run this script to refresh the MATLAB path entries related to the
%   'MSS' directory.
% - It handles multiple instances of 'MSS' directories by using the
%   first found.
% - It removes obsolete MSS entries, including folders that no longer exist.
% - It excludes hidden subfolders, such as .git, from the MATLAB path.
%
% Author:    Thor I. Fossen
% Date:      2024-04-26
% Revisions:
%   2026-10-25: Remove stale path entries and exclude hidden subfolders.

if exist('OCTAVE_VERSION', 'builtin')
    disp(['Add the MSS path with subfolders by starting' ...
        ' "octave --gui", then click, edit - set path']);
    return;
end

% Locate directories named 'MSS' using the MATLAB what function
mssInfo = what('MSS');

% Evaluate the presence and number of 'MSS' directories found
if isempty(mssInfo)
    error(['MSS directory not found. Add the MSS root folder to the' ...
        ' MATLAB path, then run mssPath again.']);
elseif length(mssInfo) > 1
    % Handle multiple MSS directories found by notifying the user
    % and using the first one
    warning(['%d MSS directories found. Using the first one' ...
        ' found at: %s'], length(mssInfo), char(mssInfo(1).path));
    % Use the first directory found, ensure it's char
    basePath = char(mssInfo(1).path);  
else
    % Confirm when one MSS path is found
    basePath = char(mssInfo.path);  % Ensure the path is char
    fprintf('MSS directory found at: %s\n', basePath);
end

% Remove the saved MSS entries, rather than scanning the current directory
% tree. A scan cannot find entries for folders that have already disappeared.
currentFolders = strsplit(path, pathsep);
if ispc
    isMssFolder = strcmpi(currentFolders, basePath) | ...
        strncmpi(currentFolders, [basePath filesep], length(basePath) + 1);
else
    isMssFolder = strcmp(currentFolders, basePath) | ...
        strncmp(currentFolders, [basePath filesep], length(basePath) + 1);
end
path(strjoin(currentFolders(~isMssFolder), pathsep));

% Retain genpath's handling of private, class, and namespace folders, but
% exclude hidden subfolders and their descendants (for example .git/refs).
newFolders = strsplit(genpath(basePath), pathsep);
newFolders = newFolders(~cellfun('isempty', newFolders));
relativeFolders = cellfun(@(folder) folder(length(basePath) + 1:end), ...
    newFolders, 'UniformOutput', false);
isHiddenFolder = ~cellfun('isempty', ...
    regexp(relativeFolders, '(^|[\\/])\.', 'once'));
addpath(strjoin(newFolders(~isHiddenFolder), pathsep));

% Save the updated path
savepath;
fprintf('MATLAB path updated and saved successfully in MATLAB.\n');
